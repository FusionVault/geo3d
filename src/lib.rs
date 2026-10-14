//! **3D positions on Earth, made simple.**
//!
//! Turn GPS coordinates (latitude, longitude, height) into 3D space and answer everyday questions:
//! where is that plane relative to me, where do I point the antenna, how far apart are they, can
//! they see each other, will they collide? Pure Rust, no dependencies.
//!
//! ```
//! use geo3d::{Geodetic, LocalFrame};
//!
//! let me = Geodetic::new(51.5, -0.12, 30.0);           // latitude°, longitude°, height in metres
//! let plane = Geodetic::new(51.6, -0.12, 3_000.0);
//!
//! let look = LocalFrame::new(me).look_angles(plane);    // where to point
//! assert!(look.azimuth_deg < 0.1);                      // due north…
//! assert!(look.elevation_deg > 14.0);                   // …15° up…
//! assert!(look.range_m > 11_000.0);                     // …11.5 km away
//! ```
//!
//! # The idea
//!
//! Latitude and longitude are angles, so you cannot subtract two positions and get a distance or a
//! direction. `geo3d` converts them into frames where ordinary 3D maths works — [`Ecef`] (x, y, z
//! from the Earth's centre) and [`Enu`] (east, north, up from you) — and does the everyday
//! calculations there.
//!
//! # Coordinate types
//!
//! | Type | What it is | Units |
//! |---|---|---|
//! | [`Geodetic`] | latitude, longitude, height — what GPS gives you | degrees, metres |
//! | [`Ecef`] | x, y, z from the Earth's centre, turning with the Earth | metres |
//! | [`Enu`] / [`Ned`] | east-north-up / north-east-down from an observer | metres |
//! | [`Aer`] | azimuth, elevation, range — where to point | degrees, metres |
//! | [`Eci`] | x, y, z from the Earth's centre, not turning — for satellites (TEME) | metres |
//! | [`Vec3`] | a plain 3D vector — a direction or a velocity | metres, m/s |
//!
//! # Find what you need
//!
//! | I want to… | Use |
//! |---|---|
//! | convert latitude/longitude/height ↔ x, y, z | [`geodetic_to_ecef`], [`ecef_to_geodetic`] |
//! | know where something is relative to me | [`LocalFrame::enu_of`], [`LocalFrame::ned_of`] |
//! | know where to point, and back | [`LocalFrame::look_angles`], [`LocalFrame::destination`] |
//! | find where a line of sight hits the ground | [`ground_strike`], [`Ellipsoid::ray_intersect`] |
//! | measure a straight-line distance | [`distance3d_m`] |
//! | measure a distance along the surface | [`Ellipsoid::geodesic_inverse`] (exact), [`Ellipsoid::great_circle_distance_m`] (quick) |
//! | travel a distance on a bearing | [`Ellipsoid::geodesic_direct`], [`Ellipsoid::great_circle_destination`] |
//! | test whether a point is inside a shape | [`Sphere`], [`Cone`], [`Cylinder`], [`AltitudeBand`] (all [`Volume`]s) |
//! | check two moving objects for a near miss | [`Track::conflict`], [`closest_approach`] |
//! | turn a velocity into course, speed and climb | [`Course`], [`LocalFrame::course_of_velocity`] |
//! | check line of sight or horizon distance | [`Ellipsoid::line_of_sight`], [`Ellipsoid::horizon_distance_m`] |
//! | work with satellites in an inertial frame | the [`eci`] module |
//! | shift positions between datums | [`Helmert7`] |
//!
//! # Conventions
//!
//! - **Degrees and metres** everywhere. The few low-level helpers that use radians say so
//!   (`lat_rad`, `gmst_rad`, [`Vec3::angle_to`]).
//! - **Latitude first**: `Geodetic::new(lat, lon, height)`.
//! - **Azimuths and bearings** run clockwise from north, 0° to 360°. Longitudes come out between
//!   −180° and 180°.
//! - **Height is above the ellipsoid** (what GPS reports), not above mean sea level.
//! - **[`WGS84`] is the default.** Every calculation is also a method on [`Ellipsoid`], so any other
//!   ellipsoid works the same way: `GRS80.to_ecef(p)`.
//! - Every type is `Copy`, nothing allocates, and there is no `unsafe` code.
//!
//! # Optional features
//!
//! All off by default; the default build has no dependencies.
//!
//! | Feature | Adds |
//! |---|---|
//! | `serde` | `Serialize` / `Deserialize` on every type |
//! | `mint` | `From` / `Into` for `mint` vectors and points (and through it `cgmath` and others) |
//! | `glam` | `From` / `Into` for `glam::DVec3` |
//! | `nalgebra` | `From` / `Into` for `nalgebra::Vector3` and `Point3` |
//! | `geo-types` | [`Geodetic`] ↔ `geo_types::Point` / `Coord` (x = longitude, y = latitude, height dropped) |

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// Every Rust code block in the README is compiled and run as a doctest.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;

mod datum;
pub mod eci;
mod ellipsoid;
mod frame;
mod geodesic;
mod great_circle;
mod interop;
mod local;
mod motion;
mod position;
mod vector;
mod visibility;
mod volume;

pub use datum::Helmert7;
pub use ellipsoid::{Ellipsoid, AIRY_1830, BESSEL_1841, CLARKE_1866, GRS80, INTERNATIONAL_1924, PZ90, WGS72, WGS84};
pub use frame::LocalFrame;
pub use geodesic::Geodesic;
pub use local::{Aer, Enu, Ned};
pub use motion::{closest_approach, Course, Cpa, Track};
pub use position::{Ecef, Eci, Geodetic};
pub use vector::Vec3;
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
