//! 3D geodesy in pure `f64`, zero dependencies.
//!
//! Real 3D geometry cannot be done on geodetic `(lat, lon, height)` directly — it is not a Euclidean
//! space — so this crate provides the metric frames and the closed-form primitives on top of them:
//!
//! - [`Ellipsoid`] — a reference ellipsoid ([`WGS84`] by default) with geodetic ↔ ECEF conversion and
//!   ray intersection.
//! - [`Geodetic`], [`Ecef`], [`Enu`], [`Ned`], [`Aer`] — typed coordinates, so latitude and longitude
//!   or east and north cannot be swapped by accident.
//! - [`LocalFrame`] — a topocentric East-North-Up frame anchored at a reference point. It precomputes
//!   its rotation once, so converting many points against one reference is a few multiplies each.
//!   Look angles ([`LocalFrame::look_angles`]) and their inverse ([`LocalFrame::destination`]) live here.
//! - [`ground_strike`] / [`Ellipsoid::ray_intersect`] — where a pointing ray meets the surface.
//! - [`Volume`]s — [`AltitudeBand`], [`Sphere`], [`Cone`], [`Cylinder`]; compose with any 2D
//!   footprint test.
//! - [`Track`] / [`closest_approach`] — dead reckoning and first-order conflict detection;
//!   [`Course`] converts between ECEF velocities and course-over-ground/speed/climb.
//! - Great circles, radii of curvature, horizon distance and geometric line of sight on
//!   [`Ellipsoid`]; [`eci`] for inertial ↔ Earth-fixed via sidereal time.
//! - [`Helmert7`] — 7-parameter datum shifts.
//!
//! Everything is `Copy`, allocation-free and `#![forbid(unsafe_code)]`. Angles in the public API are
//! **degrees**, lengths **metres**; ECEF is the frame CesiumJS renders in, so a `(lat, lon, height)`
//! track maps straight to what a globe draws.
//!
//! ```
//! use geo3d::{Geodetic, LocalFrame, Sphere, Volume, distance3d_m};
//!
//! let site = Geodetic::new(34.05, -118.24, 100.0);
//! let aircraft = Geodetic::new(34.30, -118.10, 3_000.0);
//!
//! // Straight-line (slant) distance through space, not along the surface.
//! let slant = distance3d_m(site, aircraft);
//! assert!(slant > 30_000.0 && slant < 32_000.0);
//!
//! // Where to point: azimuth from north, elevation above the horizon, range.
//! let frame = LocalFrame::new(site);
//! let aer = frame.look_angles(aircraft);
//! assert!(aer.azimuth_deg > 20.0 && aer.azimuth_deg < 30.0 && aer.elevation_deg > 4.0);
//!
//! // ...and back: the destination of that pointing is the aircraft.
//! let back = frame.destination(aer);
//! assert!((back.lat_deg - aircraft.lat_deg).abs() < 1e-9);
//!
//! // A 40 km proximity bubble around the site contains it.
//! assert!(Sphere::new(site, 40_000.0).contains(aircraft));
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod datum;
pub mod eci;
mod ellipsoid;
mod frame;
mod kinematics;
mod types;
mod volume;

pub use datum::Helmert7;
pub use ellipsoid::{
    Ellipsoid, Geodesic, AIRY_1830, BESSEL_1841, CLARKE_1866, GRS80, INTERNATIONAL_1924, PZ90, WGS72, WGS84,
};
pub use frame::{Course, LocalFrame};
pub use kinematics::{closest_approach, Cpa, Track};
pub use types::{Aer, Ecef, Eci, Enu, Geodetic, Ned, Vec3};
pub use volume::{AltitudeBand, Cone, Cylinder, Sphere, Volume};

/// WGS84 geodetic → ECEF. Shorthand for [`WGS84.to_ecef`](Ellipsoid::to_ecef).
#[inline]
pub fn geodetic_to_ecef(g: Geodetic) -> Ecef {
    WGS84.to_ecef(g)
}

/// WGS84 ECEF → geodetic. Shorthand for [`WGS84.to_geodetic`](Ellipsoid::to_geodetic).
#[inline]
pub fn ecef_to_geodetic(e: Ecef) -> Geodetic {
    WGS84.to_geodetic(e)
}

/// Straight-line **slant range** between two WGS84 geodetic points, in metres — the through-space
/// distance between two things at height (the ECEF chord), unlike a 2D great-circle distance.
#[inline]
pub fn distance3d_m(a: Geodetic, b: Geodetic) -> f64 {
    WGS84.to_ecef(a).distance_to(WGS84.to_ecef(b))
}

/// Where a ray from `observer`, pointing along compass `azimuth_deg` (0 = north, clockwise) and
/// `elevation_deg` (0 = horizon, positive up), strikes the WGS84 surface — or `None` if it points
/// above the horizon into space. Shorthand for [`LocalFrame::ground_strike`] on a fresh frame.
#[inline]
pub fn ground_strike(observer: Geodetic, azimuth_deg: f64, elevation_deg: f64) -> Option<Geodetic> {
    LocalFrame::new(observer).ground_strike(azimuth_deg, elevation_deg)
}
