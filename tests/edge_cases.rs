//! The places where geodesy maths usually breaks: the exact poles, the antimeridian, the Earth's
//! centre, extreme heights, exact antipodes, and NaN or infinite input. Every public conversion
//! must give a finite answer for valid input, pass invalid input through as NaN or `None`, and
//! never panic.

use geo3d::*;

/// Positions at the awkward places: exact poles (also at a non-zero longitude), both sides of the
/// antimeridian, deep underground, far out in space, and ordinary points for contrast.
fn awkward() -> Vec<Geodetic> {
    vec![
        Geodetic::new(90.0, 0.0, 0.0),
        Geodetic::new(-90.0, 0.0, 100.0),
        Geodetic::new(90.0, 123.0, 5_000.0),
        Geodetic::new(0.0, 180.0, 0.0),
        Geodetic::new(0.0, -180.0, 0.0),
        Geodetic::new(10.0, 179.999_999_9, 10.0),
        Geodetic::new(-10.0, -179.999_999_9, 10.0),
        Geodetic::new(45.0, 45.0, 1.0e8),
        Geodetic::new(0.0, 0.0, -6.0e6),
        Geodetic::new(0.0, 0.0, 0.0),
        Geodetic::new(-33.87, 151.21, 20.0),
    ]
}

fn finite(g: Geodetic) -> bool {
    g.lat_deg.is_finite() && g.lon_deg.is_finite() && g.height_m.is_finite()
}

fn finite3(v: [f64; 3]) -> bool {
    v.iter().all(|x| x.is_finite())
}

/// How far apart two geodetic positions are, through space.
fn apart(a: Geodetic, b: Geodetic) -> f64 {
    distance3d_m(a, b)
}

#[test]
fn conversions_round_trip_at_every_awkward_place() {
    for g in awkward() {
        let e = geodetic_to_ecef(g);
        assert!(finite3(e.to_array()), "{g:?}");
        let back = ecef_to_geodetic(e);
        assert!(finite(back), "{g:?} -> {back:?}");
        // Position error, scaled to the size of the coordinates.
        assert!(
            apart(g, back) < 1e-6 * (1.0 + g.height_m.abs() / 6.4e6),
            "{g:?} -> {back:?}"
        );
    }
    // The Earth's centre: latitude 90°, height −b.
    let c = ecef_to_geodetic(Ecef::new(0.0, 0.0, 0.0));
    assert_eq!((c.lat_deg, c.lon_deg, c.height_m), (90.0, 0.0, -WGS84.b()));
}

#[test]
fn local_frames_and_pointing_between_every_pair() {
    let pts = awkward();
    for &obs in &pts {
        let frame = LocalFrame::new(obs);
        for &target in &pts {
            let enu = frame.enu_of(target);
            let ned = frame.ned_of(target);
            assert!(finite3(enu.into()) && finite3(ned.into()), "{obs:?} -> {target:?}");
            let aer = frame.look_angles(target);
            assert!(
                aer.azimuth_deg.is_finite() && aer.elevation_deg.is_finite() && aer.range_m.is_finite(),
                "{obs:?} -> {target:?}: {aer:?}"
            );
            assert!((0.0..360.0).contains(&aer.azimuth_deg), "{aer:?}");
            assert!((-90.0..=90.0).contains(&aer.elevation_deg), "{aer:?}");
            // Pointing there and back lands on the target.
            let there = frame.destination(aer);
            let tol = 1e-6 * (1.0 + aer.range_m / 6.4e6);
            assert!(apart(there, target) < tol, "{obs:?} -> {target:?}: {there:?}");
            assert!(apart(frame.geodetic_of(enu), target) < tol);
            assert!(apart(frame.geodetic_of_ned(ned), target) < tol);
        }
        // Pointing at the Earth's centre: straight down.
        let down = frame.look_angles_ecef(Ecef::new(0.0, 0.0, 0.0));
        assert!(
            down.elevation_deg < -80.0 && down.range_m.is_finite(),
            "{obs:?}: {down:?}"
        );
    }
}

#[test]
fn distances_and_bearings_between_every_pair() {
    let pts = awkward();
    for &a in &pts {
        for &b in &pts {
            assert!(distance3d_m(a, b).is_finite());
            let gc = WGS84.great_circle_distance_m(a, b);
            let brg = WGS84.initial_bearing_deg(a, b);
            assert!(gc.is_finite() && gc >= 0.0, "{a:?} -> {b:?}: {gc}");
            assert!((0.0..360.0).contains(&brg), "{a:?} -> {b:?}: {brg}");
            assert!(finite(WGS84.great_circle_destination(a, brg, gc)), "{a:?} -> {b:?}");
            assert!(finite(WGS84.great_circle_midpoint(a, b)), "{a:?} -> {b:?}");
            let g = WGS84.geodesic_inverse(a, b).expect("every finite pair has a geodesic");
            assert!(g.distance_m.is_finite() && g.distance_m >= 0.0, "{a:?} -> {b:?}: {g:?}");
            assert!((0.0..360.0).contains(&g.initial_bearing_deg) && (0.0..360.0).contains(&g.final_bearing_deg));
            let (end, azi) = WGS84.geodesic_direct(a, g.initial_bearing_deg, g.distance_m);
            assert!(finite(end) && azi.is_finite(), "{a:?} -> {b:?}");
            // Travelling the geodesic lands on b (on the surface: heights are ignored).
            assert!(
                apart(end.with_height(0.0), b.with_height(0.0)) < 1e-5,
                "{a:?} -> {b:?}: landed at {end:?}"
            );
            assert!(finite(WGS84.interpolate(a, b, 0.5)), "{a:?} -> {b:?}");
            // A triangle through both, which may wrap a pole or cross the antimeridian.
            let p = WGS84.polygon_area(&[a, b, Geodetic::new(-33.87, 151.21, 0.0)]).unwrap();
            assert!(
                p.perimeter_m.is_finite() && p.area_m2.abs() <= 2.6e14,
                "{a:?} -> {b:?}: {p:?}"
            );
        }
    }
    // Exact antipodes, the case Vincenty's method could not solve.
    let half = WGS84
        .geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 180.0, 0.0))
        .unwrap();
    assert!((half.distance_m - 20_003_931.458_6).abs() < 1e-3);
}

#[test]
fn rays_volumes_and_motion_never_panic() {
    let pts = awkward();
    for &obs in &pts {
        for el in [-90.0, -45.0, -0.001, 0.0, 0.001, 45.0, 90.0] {
            for az in [0.0, 90.0, 180.0, 359.999] {
                if let Some(g) = ground_strike(obs, az, el) {
                    assert!(finite(g), "{obs:?} az {az} el {el}: {g:?}");
                }
            }
        }
        let e = geodetic_to_ecef(obs);
        for &other in &pts {
            let o = geodetic_to_ecef(other);
            let _ = WGS84.line_of_sight(e, o);
            let _ = WGS84.ray_intersect(e, o.vec() - e.vec());
            let cpa = closest_approach(e, Vec3::new(1.0, 2.0, 3.0), o, Vec3::new(-3.0, 0.0, 1.0));
            assert!(cpa.time_s.is_finite() && cpa.distance_m.is_finite());
            let _ = Sphere::new(obs, 1_000.0).contains(other);
            let _ = Cone::new(obs, 10.0, -90.0, 30.0, 1.0e7).contains(other);
            let _ = Cylinder::new(obs, 5_000.0, AltitudeBand::below(1.0e4)).contains(other);
        }
        let frame = LocalFrame::new(obs);
        let v = frame.velocity_of_course(Course::new(45.0, 100.0, -5.0));
        let c = frame.course_of_velocity(v);
        assert!(
            (c.course_deg - 45.0).abs() < 1e-9 && (c.ground_speed_mps - 100.0).abs() < 1e-9,
            "{obs:?}: {c:?}"
        );
    }
}

#[test]
fn nan_and_infinity_in_give_nan_infinity_or_none_out() {
    let bad = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
    let ok = Geodetic::new(10.0, 20.0, 30.0);
    for &x in &bad {
        for g in [
            Geodetic::new(x, 0.0, 0.0),
            Geodetic::new(0.0, x, 0.0),
            Geodetic::new(0.0, 0.0, x),
        ] {
            assert!(!finite3(geodetic_to_ecef(g).to_array()), "{g:?}");
            // Heights are ignored along the surface, so an unknown (NaN) height still has a
            // distance; a bad latitude or longitude has none.
            let surface = WGS84.geodesic_inverse(ok, g);
            if g.height_m.is_finite() {
                assert!(surface.is_none(), "{g:?}");
            } else {
                assert!(surface.is_some_and(|s| s.distance_m.is_finite()), "{g:?}");
            }
            assert_eq!(
                WGS84.polygon_area(&[ok, g, Geodetic::new(0.0, 0.0, 0.0)]).is_none(),
                surface.is_none()
            );
            assert!(
                !finite(WGS84.geodesic_direct(g, 10.0, 1_000.0).0) || g.height_m.is_infinite(),
                "{g:?}"
            );
            assert!(ground_strike(g, 0.0, -10.0).is_none(), "{g:?}");
            let aer = LocalFrame::new(ok).look_angles(g);
            assert!(!(aer.azimuth_deg.is_finite() && aer.elevation_deg.is_finite() && aer.range_m.is_finite()));
        }
        for e in [
            Ecef::new(x, 0.0, 0.0),
            Ecef::new(1.0e6, x, 0.0),
            Ecef::new(1.0e6, 0.0, x),
        ] {
            assert!(!finite(ecef_to_geodetic(e)), "{e:?}");
        }
        assert!(LocalFrame::new(ok).ground_strike(x, -10.0).is_none());
        assert!(LocalFrame::new(ok).ground_strike(0.0, x).is_none());
    }
}
