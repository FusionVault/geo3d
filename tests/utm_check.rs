//! UTM/UPS against GeographicLib's own GeoConvert (version 2.5.2, run at nanometre precision on
//! geographiclib.sourceforge.io, 2026-10-14). GeographicLib uses Krüger's series to sixth order
//! (about 5 nm); geo3d to eighth, so the two agree to a few nanometres.

// The reference values are kept exactly as GeoConvert printed them.
#![allow(clippy::excessive_precision)]

use geo3d::{Geodetic, UtmUps, WGS84};

/// (latitude, longitude, requested zone (-1 = standard), zone, north, easting, northing)
#[rustfmt::skip]
const FORWARD: [(f64, f64, i32, u8, bool, f64, f64); 29] = [
    (33.3, 44.4, -1, 38, true, 444140.544918426, 3684706.355549777),
    (-33.87, 151.21, -1, 56, false, 334435.706141829, 6250816.397804656),
    (40.748433, -73.985656, -1, 18, true, 585632.081341873, 4511326.153663830),
    (51.5, -0.12, -1, 30, true, 699889.806985106, 5709362.292819306),
    (0.0, 0.0, -1, 31, true, 166021.443080540, 0.000000000),
    (60.0, 5.0, -1, 32, true, 276979.926401006, 6658157.202407252),
    (63.9, 2.9, -1, 31, true, 495091.311822963, 7085874.811547142),
    (78.0, 8.9, -1, 31, true, 636716.845988947, 8665261.549757802),
    (78.0, 9.0, -1, 33, true, 360973.603634535, 8665496.995776532),
    (78.0, 20.9, -1, 33, true, 636716.845988947, 8665261.549757802),
    (78.0, 21.0, -1, 35, true, 360973.603634535, 8665496.995776532),
    (78.0, 41.9, -1, 37, true, 567282.208954672, 8660035.492844149),
    (83.9, -179.9, -1, 1, true, 465609.168730046, 9317795.753332415),
    (-79.9, 123.0, -1, 51, false, 500000.000000000, 1129575.677877436),
    (84.0, 10.0, -1, 0, true, 2115776.050744780, 1343401.388263992),
    (-80.0, 10.0, -1, 32, false, 519384.803295973, 1118247.585192557),
    (-80.0001, 10.0, -1, 0, false, 2193259.994693194, 3096031.894362924),
    (89.999, 45.0, -1, 0, true, 2000078.505692944, 1999921.494307056),
    (90.0, 0.0, -1, 0, true, 2000000.000000000, 2000000.000000000),
    (-90.0, 0.0, -1, 0, false, 2000000.000000000, 2000000.000000000),
    (45.0, 3.0, -1, 31, true, 500000.000000000, 4982950.400226552),
    (-45.0, -177.0, -1, 1, false, 500000.000000000, 5017049.599773448),
    (10.0, 179.9999999, -1, 60, true, 828928.725084910, 1106908.854143274),
    (-1e-06, 177.0, -1, 60, false, 500000.000000000, 9999999.889469954),
    (45.0, 8.0, 31, 31, true, 894076.745034487, 4995124.532465576),
    (40.0, -2.0, 31, 31, true, 73106.697732999, 4439746.917133751),
    (-20.0, -71.0, 18, 18, false, 918683.238420747, 7783514.436631494),
    (82.0, 30.0, 0, 0, true, 2444788.988891740, 1229602.872592318),
    (-78.0, -150.0, 0, 0, false, 1331509.796517138, 842141.003205625),
];

/// Requested zones GeoConvert refuses: too far from the zone for its grid limits.
#[rustfmt::skip]
const REFUSED: [(f64, f64, u8); 5] = [
    (45.0, 15.0, 31),
    (30.0, 0.0, 32),
    (87.0, 0.0, 31),
    (75.0, 100.0, 0),
    (-75.0, -100.0, 0),
];

/// (zone, north, easting, northing) → (latitude, longitude).
#[rustfmt::skip]
const REVERSE: [(u8, bool, f64, f64, f64, f64); 6] = [
    (38, true, 444140.544918426, 3684706.355549777, 33.30000000000000, 44.40000000000000),
    (56, false, 334786.0, 6252080.0, -33.85866391652991, 151.21402286089508),
    (0, true, 2000000.0, 2000000.0, 90.00000000000000, 180.00000000000000),
    (0, false, 1500000.0, 2500000.0, -83.63731756105710, -45.00000000000000),
    (31, true, 166021.443, 0.0, 0.00000000000000, -0.00000000072279),
    (60, false, 833978.557, 9999999.0, -0.00000903483105, -179.99999999927717),
];

/// A small deterministic generator (xorshift64*), so the same points every run.
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

fn apart(a: Geodetic, b: Geodetic) -> f64 {
    WGS84.to_ecef(a).distance_to(WGS84.to_ecef(b))
}

#[test]
fn forward_matches_geographiclib() {
    let mut worst = 0.0f64;
    for (lat, lon, requested, zone, north, e, n) in FORWARD {
        let g = Geodetic::new(lat, lon, 0.0);
        let u = if requested < 0 {
            assert_eq!(UtmUps::standard_zone(g), Some(zone), "{lat} {lon}");
            WGS84.to_utm_ups(g)
        } else {
            WGS84.to_utm_ups_zone(g, requested as u8)
        }
        .unwrap_or_else(|| panic!("{lat} {lon} refused"));
        assert_eq!((u.zone, u.north), (zone, north), "{lat} {lon}");
        let d = (u.easting_m - e).abs().max((u.northing_m - n).abs());
        worst = worst.max(d);
        assert!(d < 1e-8, "{lat} {lon}: {u:?} vs {e} {n}");
    }
    eprintln!("worst difference from GeographicLib: {worst:e} m");
}

#[test]
fn refuses_what_geographiclib_refuses() {
    for (lat, lon, zone) in REFUSED {
        assert!(
            WGS84.to_utm_ups_zone(Geodetic::new(lat, lon, 0.0), zone).is_none(),
            "{lat} {lon} {zone}"
        );
    }
}

#[test]
fn reverse_matches_geographiclib() {
    for (zone, north, e, n, lat, lon) in REVERSE {
        let g = WGS84.from_utm_ups(UtmUps::new(zone, north, e, n)).unwrap();
        assert!((g.lat_deg - lat).abs() < 1e-12, "{zone} {e} {n}: {g:?}");
        if lat.abs() != 90.0 {
            let d = (g.lon_deg - lon + 540.0).rem_euclid(360.0) - 180.0;
            assert!(d.abs() < 1e-12, "{zone} {e} {n}: {g:?}");
        }
    }
}

/// Grid and back lands on the same spot, everywhere, in the standard zone and in forced zones out to
/// the grid's limits.
#[test]
fn round_trips_everywhere() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let mut worst = 0.0f64;
    for i in 0..100_000 {
        let g = Geodetic::new(rng.range(-90.0, 90.0), rng.range(-180.0, 180.0), 0.0);
        let u = if i % 4 == 0 {
            // A neighbouring zone, or UPS down to 70°.
            let z = UtmUps::standard_zone(g).unwrap();
            let forced = if z == 0 { 0 } else { (z % 60) + 1 };
            match WGS84.to_utm_ups_zone(g, forced) {
                Some(u) => u,
                None => continue,
            }
        } else {
            WGS84.to_utm_ups(g).unwrap()
        };
        let back = WGS84.from_utm_ups(u).unwrap();
        let d = apart(g, back);
        worst = worst.max(d);
        assert!(d < 5e-9, "{g:?} -> {u:?} -> {back:?}: {d} m");
    }
    eprintln!("worst round trip: {worst:e} m");
}

/// Along the central meridian, UTM's northing is 0.9996 times the distance from the equator, which
/// geo3d's geodesics measure independently.
#[test]
fn central_meridian_is_k0_times_the_meridian_arc() {
    let mut worst = 0.0f64;
    for i in 0..=84 {
        let lat = i as f64 + 0.123;
        if lat >= 84.0 {
            break;
        }
        let g = Geodetic::new(lat, 3.0, 0.0);
        // Zone 31 throughout: 56° to 64° N here is the Norway exception's zone 32 by default.
        let u = WGS84.to_utm_ups_zone(g, 31).unwrap();
        let arc = WGS84
            .geodesic_inverse(Geodetic::new(0.0, 3.0, 0.0), g)
            .unwrap()
            .distance_m;
        assert_eq!(u.easting_m, 500_000.0);
        let d = (u.northing_m - 0.9996 * arc).abs();
        worst = worst.max(d);
        // A few rounding steps on numbers up to 9 300 km.
        assert!(d < 1e-15 * arc + 1e-9, "{lat}: {} vs {}", u.northing_m, 0.9996 * arc);
    }
    eprintln!("worst gap from the geodesic meridian arc: {worst:e} m");
}

/// The scale factor and convergence the projection reports match the map itself: step 1 m north and
/// 1 m east (by geodesic) and see how far, and which way, the point moves on the map. Equal scale both
/// ways at right angles also shows the map is conformal.
#[test]
fn scale_and_convergence_match_finite_differences() {
    let tm = geo3d::TransverseMercator::UTM;
    let mut rng = Rng(42);
    for _ in 0..2_000 {
        let lon0 = 3.0;
        let g = Geodetic::new(rng.range(-80.0, 84.0), lon0 + rng.range(-12.0, 12.0), 0.0);
        let p = tm.forward(lon0, g);
        let step = |bearing: f64, s: f64| {
            let (q, _) = WGS84.geodesic_direct(g, bearing, s);
            let m = tm.forward(lon0, q);
            (m.x_m, m.y_m)
        };
        let ((xn1, yn1), (xn0, yn0)) = (step(0.0, 1.0), step(0.0, -1.0));
        let ((xe1, ye1), (xe0, ye0)) = (step(90.0, 1.0), step(90.0, -1.0));
        let (dxn, dyn_, dxe, dye) = (
            (xn1 - xn0) / 2.0,
            (yn1 - yn0) / 2.0,
            (xe1 - xe0) / 2.0,
            (ye1 - ye0) / 2.0,
        );
        let (kn, ke) = (dxn.hypot(dyn_), dxe.hypot(dye));
        assert!(
            (kn - p.scale).abs() < 2e-8 && (ke - p.scale).abs() < 2e-8,
            "{g:?}: {kn} {ke} vs {}",
            p.scale
        );
        assert!((dxn * dxe + dyn_ * dye).abs() < 2e-8, "{g:?}: not at right angles");
        let gamma = -dxn.atan2(dyn_).to_degrees();
        assert!(
            (gamma - p.convergence_deg).abs() < 2e-6,
            "{g:?}: {gamma} vs {}",
            p.convergence_deg
        );
    }
}

#[test]
fn refuses_bad_input() {
    for bad in [f64::NAN, f64::INFINITY, 90.5] {
        assert!(WGS84.to_utm_ups(Geodetic::new(bad, 10.0, 0.0)).is_none());
        assert!(UtmUps::standard_zone(Geodetic::new(bad, 10.0, 0.0)).is_none());
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(WGS84.to_utm_ups(Geodetic::new(10.0, bad, 0.0)).is_none());
    }
    assert!(WGS84.to_utm_ups_zone(Geodetic::new(10.0, 10.0, 0.0), 61).is_none());
    // UPS only from 70° towards a pole; UTM only within 60° of the zone's central meridian.
    assert!(WGS84.to_utm_ups_zone(Geodetic::new(60.0, 10.0, 0.0), 0).is_none());
    assert!(WGS84.to_utm_ups_zone(Geodetic::new(0.0, 100.0, 0.0), 31).is_none());
    for u in [
        UtmUps::new(61, true, 500_000.0, 0.0),
        UtmUps::new(31, true, f64::NAN, 0.0),
        UtmUps::new(31, true, -200_000.0, 0.0),
        UtmUps::new(31, false, 500_000.0, 20_000_000.0),
        UtmUps::new(0, true, 3_000_000.0, 2_000_000.0),
    ] {
        assert!(WGS84.from_utm_ups(u).is_none(), "{u:?}");
    }
    // Central meridians.
    assert_eq!(UtmUps::new(31, true, 0.0, 0.0).central_meridian_deg(), 3.0);
    assert!(UtmUps::new(0, true, 0.0, 0.0).central_meridian_deg().is_nan());
}

/// Any ellipsoid works, and a round trip on one comes back too.
#[test]
fn other_ellipsoids() {
    for ell in [
        geo3d::GRS80,
        geo3d::INTERNATIONAL_1924,
        geo3d::Ellipsoid::new(6_378_137.0, 1.0 / 150.0),
    ] {
        let g = Geodetic::new(48.2, 16.37, 0.0);
        let u = ell.to_utm_ups(g).unwrap();
        assert_eq!(u.zone, 33);
        let back = ell.from_utm_ups(u).unwrap();
        assert!(
            (back.lat_deg - g.lat_deg).abs() < 1e-13 && (back.lon_deg - g.lon_deg).abs() < 1e-13,
            "{ell:?}: {back:?}"
        );
    }
}
