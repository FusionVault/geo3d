//! Map projections: [`TransverseMercator`] and [`PolarStereographic`], and the UTM/UPS grid built
//! on them ([`UtmUps`]).

mod polar_stereographic;
mod transverse_mercator;
mod utm_ups;

pub use polar_stereographic::PolarStereographic;
pub use transverse_mercator::TransverseMercator;
pub use utm_ups::UtmUps;

/// A point on a map projection, with the map's distortion there.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MapPoint {
    /// Easting from the projection's origin, metres (no false easting).
    pub x_m: f64,
    /// Northing from the projection's origin, metres (no false northing).
    pub y_m: f64,
    /// Grid convergence: the bearing of grid north (the y axis) clockwise from true north, degrees.
    pub convergence_deg: f64,
    /// Scale factor: a short distance on the map divided by the same distance on the ellipsoid.
    pub scale: f64,
}
