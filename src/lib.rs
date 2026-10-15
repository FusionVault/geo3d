//! # 🌍 geo3d
//!
//! **3D positions on Earth, made simple.**
//!
//! Give geo3d a GPS position, then ask everyday questions:
//!
//! - 📍 Where is that plane, seen from here?
//! - 📡 Which way should the antenna point?
//! - 📏 How far apart are they?
//! - 🗺️ How big is that field?
//! - 👀 Can they see each other?
//! - 💥 Are they going to collide?
//! - 🎯 How fast is it closing, and how soon can we reach it?
//! - 🌊 How high above sea level is it?
//!
//! It's pure Rust, with no dependencies.
//!
//! ```
//! use geo3d::{Geodetic, LocalFrame};
//!
//! let me = Geodetic::new(51.5, -0.12, 30.0);          // latitude, longitude, height in metres
//! let plane = Geodetic::new(51.6, -0.12, 3_000.0);
//!
//! let look = LocalFrame::new(me).look_angles(plane);   // which way to point
//! // due north, 14.9° up, 11.5 km away
//! assert!(look.elevation_deg > 14.0 && look.range_m > 11_000.0);
//! ```
//!
//! # ✨ Why geo3d?
//!
//! - 🎯 **Accurate.** Conversions are exact, and distances along the ground are good to 15
//!   nanometres. Everything is tested against GeographicLib, the reference library.
//! - 🛡️ **Reliable.** It handles the tricky places, like the poles and opposite sides of the Earth,
//!   without NaN or panics.
//! - 🧩 **Easy to use.** Degrees and metres everywhere, with named fields. Latitude always comes first.
//! - 🪶 **Light and fast.** No dependencies, no `unsafe` code, and nothing allocates. Works on Rust
//!   1.85 and newer, and on microcontrollers too, without the standard library.
//!
//! # 🗺️ The basics
//!
//! Latitude and longitude are angles, so you can't subtract two of them to get a distance. geo3d
//! turns positions into frames where normal 3D maths works.
//!
//! | Type | What it holds | Units |
//! |---|---|---|
//! | [`Geodetic`] | latitude, longitude and height, like GPS gives you | degrees, metres |
//! | [`Ecef`] | x, y, z from the centre of the Earth | metres |
//! | [`Enu`] / [`Ned`] | how far east, north and up (or north, east and down) from you | metres |
//! | [`Aer`] | azimuth, elevation and range: which way to point | degrees, metres |
//! | [`Eci`] | x, y, z that don't spin with the Earth, for satellites | metres |
//! | [`Vec3`] | a plain 3D vector, like a velocity | metres, m/s |
//! | [`Cov3`] | how uncertain a position is (a covariance) | metres², degrees² |
//! | [`Attitude`] | which way a vehicle faces: yaw, pitch and roll | degrees |
//! | [`Frd`] | an offset in a vehicle's own frame: forward, right, down | metres |
//! | [`Quaternion`] | the same turn as yaw, pitch and roll, as four numbers | none |
//! | [`UtmUps`] | a UTM or UPS grid position: zone, easting and northing | metres |
//!
//! A few rules hold everywhere:
//!
//! - Latitude comes first: `Geodetic::new(lat, lon, height)`.
//! - Bearings go clockwise from north, from 0° to 360°.
//! - Heights are above the WGS84 ellipsoid, like GPS, not above sea level. A [`Geoid`] converts
//!   between the two.
//! - [`WGS84`] is the default, but any [`Ellipsoid`] works: `GRS80.to_ecef(p)`.
//! - A few low-level helpers use radians, and their names say so (`lat_rad`, `gmst_rad`).
//!   [`Vec3::angle_to`] also returns radians.
//!
//! # 🖨️ Printing
//!
//! Every coordinate type implements `Display`. Plain `{}` keeps every digit (it reads back as
//! exactly the same number), `{:.N}` gives N decimals, and `{:#}` writes hemispheres as letters:
//!
//! ```
//! use geo3d::Geodetic;
//! let sydney = Geodetic::new(-33.87, 151.21, 20.0);
//! assert_eq!(format!("{sydney}"), "-33.87°, 151.21°, 20 m");
//! assert_eq!(format!("{sydney:.5}"), "-33.87000°, 151.21000°, 20.0 m");
//! assert_eq!(format!("{sydney:#}"), "33.87° S, 151.21° E, 20 m");
//! ```
//!
//! Where a type mixes degrees and metres, the metres get four fewer decimals: 0.0001° is about
//! 10 m on the ground.
//!
//! # 🛡️ Tricky inputs
//!
//! - The poles, the antimeridian, the Earth's centre and points far out in space all give finite,
//!   correct answers.
//! - Azimuths and bearings are always in `[0, 360)`, and longitudes come out between −180° and 180°.
//! - NaN or infinite input gives NaN out, or `None` where a function returns an `Option`.
//! - Nothing panics.
//!
//! `tests/edge_cases.rs` checks all of this for every public conversion.
//!
//! # 🧭 Find what you need
//!
//! | I want to… | Use |
//! |---|---|
//! | convert between latitude/longitude and x, y, z | [`geodetic_to_ecef`], [`ecef_to_geodetic`] |
//! | turn a GPS height into a height above sea level, and back | [`Geoid`], [`GeoidGrid`] |
//! | use UTM/UPS grid coordinates, or another transverse Mercator grid | [`Ellipsoid::to_utm_ups`], [`Ellipsoid::from_utm_ups`], [`TransverseMercator`], [`PolarStereographic`] |
//! | find where something is, relative to me | [`LocalFrame::enu_of`], [`LocalFrame::ned_of`] |
//! | work out which way to point, and back | [`LocalFrame::look_angles`], [`LocalFrame::destination`] |
//! | find where a line of sight hits the ground | [`ground_strike`], [`Ellipsoid::ray_intersect`] |
//! | measure a straight line through space | [`distance3d_m`] |
//! | measure distance along the ground | [`Ellipsoid::geodesic_inverse`] (exact), [`Ellipsoid::great_circle_distance_m`] (quick) |
//! | travel a distance on a bearing | [`Ellipsoid::geodesic_direct`], [`Ellipsoid::great_circle_destination`] |
//! | measure the area of a field or region | [`Ellipsoid::polygon_area`] |
//! | follow a route: how far off it, how far along, where routes cross | [`Ellipsoid::cross_track_distance_m`], [`Ellipsoid::along_track_distance_m`], [`Ellipsoid::great_circle_intersection`] |
//! | find the centre of a group, or a point part-way along a route | [`Ellipsoid::mean_position`], [`Ellipsoid::great_circle_interpolate`] |
//! | check if something is inside a zone | [`Sphere`], [`Cone`], [`Cylinder`], [`AltitudeBand`] (all are [`Volume`]s) |
//! | check if two objects will collide | [`Track::conflict`], [`closest_approach`] |
//! | measure how fast something is closing, for radar Doppler | [`Track::range_rate_mps`], [`LocalFrame::range_rate_mps`] |
//! | plan an intercept: when, where, and the slowest speed that works | [`Track::intercept_from`], [`Track::slowest_intercept_from`] |
//! | turn a velocity into heading, speed and climb | [`Course`], [`LocalFrame::course_of_velocity`] |
//! | work with a vehicle's yaw, pitch and roll | [`Attitude`], [`LocalFrame::geodetic_of_body`], [`LocalFrame::body_to_ecef_rotation`], [`Quaternion`] |
//! | carry measurement uncertainty between frames | [`Cov3`], [`Aer::cov_to_enu`], [`LocalFrame::enu_cov_to_ecef`] |
//! | check line of sight, or how far I can see | [`Ellipsoid::line_of_sight`], [`Ellipsoid::horizon_distance_m`], and their `_with_refraction` versions for radio |
//! | track satellites, down to their Doppler shift | the [`eci`] module, [`Track::from_eci`], and `examples/satellite.rs` with the `sgp4` crate |
//! | shift between datums | [`Helmert7`] |
//!
//! The [`navlab`] page solves ten classic navigation problems step by step. The
//! [`examples`](https://github.com/FusionVault/geo3d/tree/main/examples) folder has eight complete
//! programs for real jobs, from an airport radar to a search and rescue. More worked examples, a
//! comparison with other crates and an accuracy table are in the
//! [README](https://github.com/FusionVault/geo3d#readme).
//!
//! # 🧩 Optional features
//!
//! None of them are needed: a plain `cargo add geo3d` pulls in nothing else.
//!
//! | Feature | What it adds |
//! |---|---|
//! | `std` | on by default: the standard library's maths. Nothing to download; it comes with Rust |
//! | `libm` | runs without the standard library (`no_std`), for microcontrollers |
//! | `serde` | saving and loading every type |
//! | `mint` | conversions to `mint` types, and through them to `cgmath` and others |
//! | `glam` | conversions to `glam::DVec3`, [`Cov3`] to `glam::DMat3`, and [`Quaternion`] to `glam::DQuat` |
//! | `nalgebra` | conversions to `nalgebra` vectors and points, [`Cov3`] to `Matrix3`, and [`Quaternion`] |
//! | `geo-types` | [`Geodetic`] ↔ `geo_types::Point` (x is longitude, y is latitude) |
//!
//! For a microcontroller, `cargo add geo3d --no-default-features --features libm`. With the
//! standard library, keep `std`: its maths rounds the last digit more carefully.

#![cfg_attr(not(any(feature = "std", test)), no_std)]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

#[cfg(not(any(feature = "std", feature = "libm")))]
compile_error!(
    "geo3d needs a maths library: keep the default `std` feature, or turn on `libm` to run without the standard library."
);

// Without the standard library, the floating-point maths comes from `libm`. (Unit tests always
// have the standard library.)
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
mod math;

/// Every Rust code block in the README is compiled and run as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;

mod angle;
mod attitude;
mod auxlat;
mod covariance;
mod datum;
mod display;
pub mod eci;
mod ellipsoid;
mod frame;
mod geodesic;
mod geoid;
mod great_circle;
mod interop;
mod local;
mod motion;
pub mod navlab;
mod position;
mod projection;
mod quaternion;
mod vector;
mod visibility;
mod volume;

pub use attitude::{Attitude, Frd};
pub use covariance::Cov3;
pub use datum::Helmert7;
pub use ellipsoid::{Ellipsoid, AIRY_1830, BESSEL_1841, CLARKE_1866, GRS80, INTERNATIONAL_1924, PZ90, WGS72, WGS84};
pub use frame::LocalFrame;
pub use geodesic::{Geodesic, PolygonArea};
pub use geoid::{Geoid, GeoidGrid, GeoidGridError, GeoidInterpolation};
pub use local::{Aer, Enu, Ned};
pub use motion::{closest_approach, Course, Cpa, Intercept, Track};
pub use position::{Ecef, Eci, Geodetic};
pub use projection::{MapPoint, PolarStereographic, TransverseMercator, UtmUps};
pub use quaternion::Quaternion;
pub use vector::Vec3;
pub use visibility::STANDARD_REFRACTION_K;
pub use volume::{AltitudeBand, Cone, Cylinder, Sphere, Volume};

/// Latitude/longitude/height → x, y, z on [`WGS84`]. Shorthand for [`WGS84.to_ecef`](Ellipsoid::to_ecef).
#[doc(alias = "geodetic2ecef")]
#[doc(alias = "lla2ecef")]
#[inline]
pub fn geodetic_to_ecef(g: Geodetic) -> Ecef {
    WGS84.to_ecef(g)
}

/// x, y, z → latitude/longitude/height on [`WGS84`]. Shorthand for
/// [`WGS84.to_geodetic`](Ellipsoid::to_geodetic).
#[doc(alias = "ecef2geodetic")]
#[doc(alias = "ecef2lla")]
#[inline]
pub fn ecef_to_geodetic(e: Ecef) -> Geodetic {
    WGS84.to_geodetic(e)
}

/// The straight-line distance between two positions through space, metres — the **slant range**
/// between things at height, not the distance along the surface. On [`WGS84`].
#[doc(alias = "slant_range")]
#[inline]
pub fn distance3d_m(a: Geodetic, b: Geodetic) -> f64 {
    WGS84.to_ecef(a).distance_to(WGS84.to_ecef(b))
}

/// Where the line of sight from `observer` along `azimuth_deg` (0 = north, clockwise) and
/// `elevation_deg` (0 = horizon, negative = down) first meets the ground, on [`WGS84`]. `None` if it
/// points into space or the observer is below the surface. Shorthand for
/// [`LocalFrame::ground_strike`].
#[inline]
pub fn ground_strike(observer: Geodetic, azimuth_deg: f64, elevation_deg: f64) -> Option<Geodetic> {
    LocalFrame::new(observer).ground_strike(azimuth_deg, elevation_deg)
}

#[cfg(test)]
mod test_util {
    /// `|a − b| ≤ tol`.
    pub fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// Angles compare on the circle: 359.99999999° and 0° are 1e-8° apart.
    pub fn angle_close(a: f64, b: f64, tol: f64) -> bool {
        let d = (a - b).rem_euclid(360.0);
        d.min(360.0 - d) <= tol
    }
}
