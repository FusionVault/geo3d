//! # 🧭 The NavLab examples
//!
//! FFI, the Norwegian Defence Research Establishment, publishes ten worked navigation problems
//! with the [n-vector](https://www.ffi.no/en/research/n-vector) method. They're a good tour of
//! what geo3d does, so here they are, solved with geo3d.
//!
//! Every example runs as a test. The answers are checked against the ones published by FFI and by
//! the [nvector](https://github.com/pbrod/nvector) Python package's own tests.
//!
//! Example 2 turns a vehicle's yaw, pitch and roll into a position. It needs attitude, which is on
//! the [roadmap](https://github.com/FusionVault/geo3d/issues/13).
//!
//! ## 1. Where is B, seen from A?
//!
//! A is at 1° N, 2° E, 3 m below the ellipsoid; B is at 4° N, 5° E, 6 m below. Find B's offset north,
//! east and down from A, and the direction from A to B.
//!
//! ```
//! use geo3d::{Geodetic, LocalFrame};
//!
//! let a = LocalFrame::new(Geodetic::new(1.0, 2.0, -3.0));
//! let b = Geodetic::new(4.0, 5.0, -6.0);
//!
//! let ned = a.ned_of(b);          // 331 730 m north, 332 998 m east, 17 404 m down
//! let look = a.look_angles(b);    // 45.1° from north, 2.1° below the horizon, 470 357 m away
//! assert!((ned.north - 331_730.234_780_89).abs() < 1e-6);
//! assert!((ned.east - 332_997.874_989_27).abs() < 1e-6);
//! assert!((ned.down - 17_404.271_361_94).abs() < 1e-6);
//! assert!((look.azimuth_deg - 45.109_263_24).abs() < 1e-8);
//! ```
//!
//! B is "down" even though both are near the surface: the Earth curves away over 470 km.
//!
//! ## 3. From x, y, z to latitude and longitude
//!
//! ```
//! use geo3d::{ecef_to_geodetic, Ecef};
//!
//! let b = ecef_to_geodetic(Ecef::new(6_371e3 * 0.9, -6_371e3, 6_371e3 * 1.1));
//! // 39.3787° N, 48.0128° W, 4 702 059.834 m up
//! assert!((b.lat_deg - 39.3787).abs() < 1e-4 && (b.lon_deg + 48.0128).abs() < 1e-4);
//! assert!((b.height_m - 4_702_059.834).abs() < 1e-3);
//! ```
//!
//! ## 4. From latitude and longitude to x, y, z
//!
//! ```
//! use geo3d::{geodetic_to_ecef, Geodetic};
//!
//! let p = geodetic_to_ecef(Geodetic::new(1.0, 2.0, 3.0));
//! assert!((p.x - 6_373_290.277_218_279).abs() < 1e-8);
//! assert!((p.y - 222_560.200_674_736_52).abs() < 1e-9);
//! assert!((p.z - 110_568.827_181_785_93).abs() < 1e-9);
//! ```
//!
//! These agree with the published answer to a nanometre. In fact geo3d's are slightly closer to
//! the exact values (worked out to 50 digits) in all three coordinates.
//!
//! ## 5. The distance between two points
//!
//! A is at 88° N, 0° E and B at 89° N, 170° W. FFI uses a sphere of radius 6 371 km.
//!
//! ```
//! use geo3d::{Ellipsoid, Geodetic, WGS84};
//!
//! let (a, b) = (Geodetic::new(88.0, 0.0, 0.0), Geodetic::new(89.0, -170.0, 0.0));
//! let sphere = Ellipsoid::new(6_371e3, 0.0);
//!
//! let along = sphere.great_circle_distance_m(a, b);            // 332.456 km over the surface
//! let straight = sphere.to_ecef(a).distance_to(sphere.to_ecef(b)); // 332.419 km through it
//! assert!((along - 332_456.444_11).abs() < 1e-5 && (straight - 332_418.724_86).abs() < 1e-5);
//!
//! // On the real (WGS84) shape it's 333.95 km.
//! let real = WGS84.geodesic_inverse(a, b).unwrap();
//! assert!((real.distance_m - 333_950.0).abs() < 5.0);
//! ```
//!
//! ## 6. Where was it at 16 seconds?
//!
//! An object is at 89° N, 0° E at 10 s and at 89° N, 180° E at 20 s. Where was it at 16 s?
//!
//! ```
//! use geo3d::{Geodetic, WGS84};
//!
//! let (t0, t1) = (Geodetic::new(89.0, 0.0, 0.0), Geodetic::new(89.0, 180.0, 0.0));
//! let at_16 = WGS84.great_circle_interpolate(t0, t1, (16.0 - 10.0) / (20.0 - 10.0));
//! assert!((at_16.lat_deg - 89.8).abs() < 1e-12);           // 0.2° past the pole
//! assert!((at_16.lon_deg.abs() - 180.0).abs() < 1e-9);
//! ```
//!
//! FFI's answer is 89.79998°: it blends the two n-vectors in a straight line, which is quick but
//! doesn't move at a steady speed. geo3d moves along the circle at a constant speed, which gives
//! exactly 89.8°.
//!
//! ## 7. The centre of three points
//!
//! ```
//! use geo3d::{Geodetic, WGS84};
//!
//! let points = [Geodetic::new(90.0, 0.0, 0.0), Geodetic::new(60.0, 10.0, 0.0), Geodetic::new(50.0, -20.0, 0.0)];
//! let m = WGS84.mean_position(&points).unwrap();   // 67.2362° N, 6.9175° W
//! let n = WGS84.surface_normal(m);
//! assert!((n.x - 0.384_117_170_292_6).abs() < 1e-12);
//! assert!((n.y + 0.046_602_405_485_689_447).abs() < 1e-12);
//! assert!((n.z - 0.922_107_485_757_139_5).abs() < 1e-12);
//! ```
//!
//! ## 8. Travel 1 km on a bearing
//!
//! From 80° N, 90° W, head 200° for 1 000 m, on a sphere of radius 6 371 km.
//!
//! ```
//! use geo3d::{Ellipsoid, Geodetic};
//!
//! let sphere = Ellipsoid::new(6_371e3, 0.0);
//! let b = sphere.great_circle_destination(Geodetic::new(80.0, -90.0, 0.0), 200.0, 1_000.0);
//! assert!((b.lat_deg - 79.991_548_67).abs() < 1e-8 && (b.lon_deg + 90.017_698_37).abs() < 1e-8);
//! ```
//!
//! On the real shape, use [`Ellipsoid::geodesic_direct`](crate::Ellipsoid::geodesic_direct).
//!
//! ## 9. Where do two paths cross?
//!
//! Path A runs from (10° N, 20° E) to (30° N, 40° E), and path B from (50° N, 60° E) to
//! (70° N, 80° E).
//!
//! ```
//! use geo3d::{Geodetic, WGS84};
//!
//! let g = |lat, lon| Geodetic::new(lat, lon, 0.0);
//! let c = WGS84.great_circle_intersection(g(10.0, 20.0), g(30.0, 40.0), g(50.0, 60.0), g(70.0, 80.0)).unwrap();
//! assert!((c.lat_deg - 40.318_643_07).abs() < 1e-8 && (c.lon_deg - 55.901_867_88).abs() < 1e-8);
//! ```
//!
//! ## 10. How far off the path?
//!
//! A path runs due north from (0°, 0°) to (10° N, 0°). How far to the side of it is (1° N, 0.1° E)?
//!
//! ```
//! use geo3d::{Ellipsoid, Geodetic};
//!
//! let sphere = Ellipsoid::new(6_371e3, 0.0);
//! let (start, end) = (Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(10.0, 0.0, 0.0));
//! let off = sphere.cross_track_distance_m(start, end, Geodetic::new(1.0, 0.1, 0.0));
//! assert!((off - 11_117.799_110_15).abs() < 1e-6);   // 11.1 km to the right
//! ```
//!
//! [`Ellipsoid::along_track_distance_m`](crate::Ellipsoid::along_track_distance_m) gives how far
//! along the path that point is.
