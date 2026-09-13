//! Cross-validation of the frame conversions against an independent implementation (`map_3d`) over a
//! global grid. Tolerances are set at the level the two closed forms should agree to in f64; a
//! regression in either the forward or inverse conversion, the local frame or the look angles fails
//! here even if the internal round-trip tests still pass.

use geo3d::{Geodetic, LocalFrame, WGS84};
use map_3d::Ellipsoid as M;

fn grid() -> Vec<Geodetic> {
    let mut pts = Vec::new();
    let mut lat = -88.0;
    while lat <= 88.0 {
        let mut lon = -180.0;
        while lon < 180.0 {
            for h in [-500.0, 0.0, 3_000.0, 12_000.0, 400_000.0, 20_200_000.0] {
                pts.push(Geodetic::new(lat, lon, h));
            }
            lon += 22.5;
        }
        lat += 8.0;
    }
    pts
}

#[test]
fn geodetic_to_ecef_agrees_to_a_micrometre() {
    for g in grid() {
        let ours = WGS84.to_ecef(g);
        let (x, y, z) = map_3d::geodetic2ecef(g.lat_rad(), g.lon_rad(), g.height_m, M::WGS84);
        assert!(
            (ours.x - x).abs() < 1e-6 && (ours.y - y).abs() < 1e-6 && (ours.z - z).abs() < 1e-6,
            "{g:?}: {ours:?} vs ({x},{y},{z})"
        );
    }
}

#[test]
fn ecef_to_geodetic_agrees_where_map_3d_is_exact() {
    for g in grid() {
        // map_3d's inverse is a fixed-count iteration: exact near the surface, but at GNSS-orbit
        // height it drifts by up to ~4e-5° (ours round-trips to 1e-11° there — see the ellipsoid
        // grid test), so only the terrestrial and LEO points are compared here, tightly.
        if g.height_m > 100_000.0 {
            continue;
        }
        let e = WGS84.to_ecef(g);
        let ours = WGS84.to_geodetic(e);
        let (lat, lon, h) = map_3d::ecef2geodetic(e.x, e.y, e.z, M::WGS84);
        assert!(
            (ours.lat_deg - lat.to_degrees()).abs() < 1e-8
                && (ours.lon_deg - lon.to_degrees()).abs() < 1e-8
                && (ours.height_m - h).abs() < 1e-3,
            "{g:?}: {ours:?} vs ({},{},{h})",
            lat.to_degrees(),
            lon.to_degrees()
        );
    }
}

#[test]
fn enu_and_look_angles_agree() {
    let observers = [
        Geodetic::new(0.0, 0.0, 0.0),
        Geodetic::new(51.5, -0.12, 30.0),
        Geodetic::new(-33.9, 151.2, 5.0),
        Geodetic::new(80.0, 170.0, 1000.0),
    ];
    for obs in observers {
        let frame = LocalFrame::new(obs);
        for g in grid() {
            let ours = frame.enu_of(g);
            let (e, n, u) = map_3d::geodetic2enu(
                g.lat_rad(),
                g.lon_rad(),
                g.height_m,
                obs.lat_rad(),
                obs.lon_rad(),
                obs.height_m,
                M::WGS84,
            );
            assert!(
                (ours.east - e).abs() < 1e-5 && (ours.north - n).abs() < 1e-5 && (ours.up - u).abs() < 1e-5,
                "{obs:?}->{g:?}: {ours:?} vs ({e},{n},{u})"
            );
            let aer = frame.look_angles(g);
            let (az, el, r) = map_3d::geodetic2aer(
                g.lat_rad(),
                g.lon_rad(),
                g.height_m,
                obs.lat_rad(),
                obs.lon_rad(),
                obs.height_m,
                M::WGS84,
            );
            let az = az.to_degrees().rem_euclid(360.0);
            let d = (aer.azimuth_deg - az).rem_euclid(360.0);
            let daz = d.min(360.0 - d);
            assert!(daz < 1e-9 || r < 1.0, "{obs:?}->{g:?}: az {} vs {az}", aer.azimuth_deg);
            assert!(
                (aer.elevation_deg - el.to_degrees()).abs() < 1e-9 && (aer.range_m - r).abs() < 1e-5 * r.max(1.0),
                "{obs:?}->{g:?}: {aer:?} vs ({az},{},{r})",
                el.to_degrees()
            );
        }
    }
}
