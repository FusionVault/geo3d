//! Cross-validation of the geodesics against geographiclib-rs (GeographicLib's algorithms in Rust,
//! itself validated against GeographicLib's GeodTest set): random pairs worldwide, nearly antipodal
//! pairs (where Vincenty's method fails), and the special lines along the equator and meridians.

use geo3d::{Ellipsoid, Geodetic, GRS80, WGS84};
use geographiclib_rs::{DirectGeodesic, Geodesic, InverseGeodesic};

/// A small deterministic generator (xorshift64*), so the test set is the same on every run.
struct Rng(u64);
impl Rng {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// Angles compare on the circle.
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    d.min(360.0 - d)
}

fn pairs() -> Vec<(Geodetic, Geodetic)> {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut v = Vec::new();
    for _ in 0..20_000 {
        let a = Geodetic::new(rng.range(-90.0, 90.0), rng.range(-180.0, 180.0), 0.0);
        let b = Geodetic::new(rng.range(-90.0, 90.0), rng.range(-180.0, 180.0), 0.0);
        v.push((a, b));
    }
    // Nearly antipodal: the second point within half a degree of the first one's antipode.
    for _ in 0..5_000 {
        let a = Geodetic::new(rng.range(-90.0, 90.0), rng.range(-180.0, 180.0), 0.0);
        let b = Geodetic::new(
            (-a.lat_deg + rng.range(-0.5, 0.5)).clamp(-90.0, 90.0),
            a.lon_deg + 180.0 + rng.range(-0.5, 0.5),
            0.0,
        )
        .normalized();
        v.push((a, b));
    }
    // Equatorial, meridional and polar lines.
    for lon in [0.0, 1.0, 90.0, 179.0, 179.5, 179.9, 180.0] {
        v.push((Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, lon, 0.0)));
        v.push((Geodetic::new(10.0, 0.0, 0.0), Geodetic::new(-10.0, lon, 0.0)));
        v.push((Geodetic::new(90.0, 0.0, 0.0), Geodetic::new(-30.0, lon, 0.0)));
    }
    v
}

fn check(ell: Ellipsoid) {
    let reference = Geodesic::new(ell.a(), ell.f());
    let (mut worst_s, mut worst_az) = (0.0f64, 0.0f64);
    for (a, b) in pairs() {
        let ours = ell.geodesic_inverse(a, b).expect("every finite pair has a geodesic");
        let (s12, azi1, azi2, _a12): (f64, f64, f64, f64) =
            reference.inverse(a.lat_deg, a.lon_deg, b.lat_deg, b.lon_deg);
        let ds = (ours.distance_m - s12).abs();
        worst_s = worst_s.max(ds);
        assert!(ds < 1e-7, "{a:?} -> {b:?}: s12 {} vs {s12}", ours.distance_m);
        // Bearings are ill-conditioned only for nearly antipodal or coincident points.
        if s12 > 1.0 && s12 < 0.9 * std::f64::consts::PI * ell.b() {
            let daz = angle_diff(ours.initial_bearing_deg, azi1).max(angle_diff(ours.final_bearing_deg, azi2));
            worst_az = worst_az.max(daz);
            assert!(daz < 1e-9, "{a:?} -> {b:?}: azimuths off by {daz}°");
        }
        // The direct problem from our answer lands where geographiclib-rs's does.
        let (end, azi) = ell.geodesic_direct(a, ours.initial_bearing_deg, ours.distance_m);
        let (lat2, lon2, azi2r): (f64, f64, f64) =
            reference.direct(a.lat_deg, a.lon_deg, ours.initial_bearing_deg, ours.distance_m);
        let miss = ell
            .to_ecef(end)
            .distance_to(ell.to_ecef(Geodetic::new(lat2, lon2, 0.0)));
        assert!(
            miss < 1e-7,
            "{a:?} az {} s {}: direct differs by {miss} m",
            ours.initial_bearing_deg,
            ours.distance_m
        );
        assert!(
            angle_diff(azi, azi2r) < 1e-9 || s12 < 1.0,
            "{a:?}: final azimuth {azi} vs {azi2r}"
        );
    }
    eprintln!("{ell:?}: worst distance difference {worst_s:e} m, worst azimuth difference {worst_az:e}°");
}

#[test]
fn inverse_and_direct_match_geographiclib_on_wgs84() {
    check(WGS84);
}

#[test]
fn inverse_and_direct_match_geographiclib_on_other_ellipsoids() {
    check(GRS80);
    check(Ellipsoid::new(6_378_137.0, 1.0 / 150.0)); // strongly oblate
    check(Ellipsoid::new(6_378_137.0, -1.0 / 150.0)); // prolate
}
