//! Closed-form 3D volumes. Domain-neutral: an altitude band is airspace or a warehouse mezzanine; a
//! cone is a camera field of view, an antenna beam or a spotlight; a sphere is a delivery radius or a
//! proximity alert.

use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::frame::LocalFrame;
use crate::types::{Aer, Ecef, Geodetic, Vec3};

/// Something that can decide whether a geodetic point is inside it.
pub trait Volume {
    /// Whether `point` lies inside (boundary inclusive).
    fn contains(&self, point: Geodetic) -> bool;
}

/// A vertical band of ellipsoidal heights, metres. A `None` bound is open on that side, so a
/// floor-only or ceiling-only band works. Compose with any 2D footprint test for an extruded prism:
///
/// ```
/// use geo3d::AltitudeBand;
/// let band = AltitudeBand::new(Some(5_000.0), Some(12_000.0));
/// let footprint_contains = true; // from your 2D geometry library
/// assert!(band.within_prism(footprint_contains, 8_000.0));
/// assert!(!band.within_prism(footprint_contains, 3_000.0));
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AltitudeBand {
    /// Lowest height inside, metres; `None` = unbounded below.
    pub floor_m: Option<f64>,
    /// Highest height inside, metres; `None` = unbounded above.
    pub ceiling_m: Option<f64>,
}

impl AltitudeBand {
    /// A band between two optional bounds.
    #[inline]
    pub const fn new(floor_m: Option<f64>, ceiling_m: Option<f64>) -> Self {
        AltitudeBand { floor_m, ceiling_m }
    }

    /// Everything at or above `floor_m`.
    #[inline]
    pub const fn above(floor_m: f64) -> Self {
        AltitudeBand::new(Some(floor_m), None)
    }

    /// Everything at or below `ceiling_m`.
    #[inline]
    pub const fn below(ceiling_m: f64) -> Self {
        AltitudeBand::new(None, Some(ceiling_m))
    }

    /// Whether a height (metres) is inside the band.
    #[inline]
    pub fn contains_height(&self, height_m: f64) -> bool {
        height_m >= self.floor_m.unwrap_or(f64::NEG_INFINITY) && height_m <= self.ceiling_m.unwrap_or(f64::INFINITY)
    }

    /// The extruded-prism test: inside the 2D footprint **and** inside the band.
    #[inline]
    pub fn within_prism(&self, footprint_contains: bool, height_m: f64) -> bool {
        footprint_contains && self.contains_height(height_m)
    }
}

impl Volume for AltitudeBand {
    #[inline]
    fn contains(&self, point: Geodetic) -> bool {
        self.contains_height(point.height_m)
    }
}

/// A sphere about a geodetic centre — a range/coverage bubble. The centre's ECEF position is
/// computed once at construction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sphere {
    ellipsoid: Ellipsoid,
    center: Geodetic,
    center_ecef: Ecef,
    radius_m: f64,
}

impl Sphere {
    /// A sphere on [`WGS84`].
    #[inline]
    pub fn new(center: Geodetic, radius_m: f64) -> Self {
        Self::on(WGS84, center, radius_m)
    }

    /// A sphere on the given ellipsoid.
    pub fn on(ellipsoid: Ellipsoid, center: Geodetic, radius_m: f64) -> Self {
        Sphere {
            ellipsoid,
            center,
            center_ecef: ellipsoid.to_ecef(center),
            radius_m,
        }
    }

    /// The centre.
    #[inline]
    pub const fn center(&self) -> Geodetic {
        self.center
    }

    /// The radius, metres.
    #[inline]
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }

    /// Whether an ECEF point is inside.
    #[inline]
    pub fn contains_ecef(&self, p: Ecef) -> bool {
        self.center_ecef.distance_to(p) <= self.radius_m
    }
}

impl Volume for Sphere {
    #[inline]
    fn contains(&self, point: Geodetic) -> bool {
        self.contains_ecef(self.ellipsoid.to_ecef(point))
    }
}

/// A right-circular cone: apex at a reference point, axis along an azimuth (degrees from north,
/// clockwise) and elevation (degrees above the horizon), a half-angle spread and a slant range. A
/// directional coverage volume — sector, field of view, beam, scan arc. The apex frame, axis and
/// the half-angle cosine are computed once at construction.
///
/// ```
/// use geo3d::{Cone, Geodetic, Volume};
/// // A camera on a 100 m mast looking due east, 20° half-angle, 5 km reach.
/// let cam = Cone::new(Geodetic::new(0.0, 0.0, 100.0), 90.0, 0.0, 20.0, 5_000.0);
/// assert!(cam.contains(Geodetic::new(0.0, 0.02, 100.0)));   // ~2.2 km east
/// assert!(!cam.contains(Geodetic::new(0.02, 0.0, 100.0)));  // ~2.2 km north: outside the spread
/// assert!(!cam.contains(Geodetic::new(0.0, 0.1, 100.0)));   // ~11 km east: beyond the range
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cone {
    frame: LocalFrame,
    axis_enu: Vec3,
    cos_half_angle: f64,
    azimuth_deg: f64,
    elevation_deg: f64,
    half_angle_deg: f64,
    range_m: f64,
}

impl Cone {
    /// A cone on [`WGS84`].
    #[inline]
    pub fn new(apex: Geodetic, azimuth_deg: f64, elevation_deg: f64, half_angle_deg: f64, range_m: f64) -> Self {
        Self::on(WGS84, apex, azimuth_deg, elevation_deg, half_angle_deg, range_m)
    }

    /// A cone on the given ellipsoid.
    pub fn on(
        ellipsoid: Ellipsoid,
        apex: Geodetic,
        azimuth_deg: f64,
        elevation_deg: f64,
        half_angle_deg: f64,
        range_m: f64,
    ) -> Self {
        Cone {
            frame: LocalFrame::on(ellipsoid, apex),
            axis_enu: Aer::new(azimuth_deg, elevation_deg, 1.0).direction_enu().vec(),
            cos_half_angle: half_angle_deg.to_radians().cos(),
            azimuth_deg,
            elevation_deg,
            half_angle_deg,
            range_m,
        }
    }

    /// The apex.
    #[inline]
    pub const fn apex(&self) -> Geodetic {
        self.frame.origin()
    }

    /// The axis as a pointing (range = the cone's range).
    #[inline]
    pub const fn axis(&self) -> Aer {
        Aer::new(self.azimuth_deg, self.elevation_deg, self.range_m)
    }

    /// The half-angle, degrees.
    #[inline]
    pub const fn half_angle_deg(&self) -> f64 {
        self.half_angle_deg
    }

    /// Whether an ECEF point is inside: within range of the apex and within the half-angle of the axis.
    pub fn contains_ecef(&self, p: Ecef) -> bool {
        let v = self.frame.to_enu(p).vec();
        let r = v.norm();
        if r > self.range_m {
            return false;
        }
        if r == 0.0 {
            return true; // the apex itself
        }
        v.dot(self.axis_enu) / r >= self.cos_half_angle
    }
}

impl Volume for Cone {
    #[inline]
    fn contains(&self, point: Geodetic) -> bool {
        self.contains_ecef(self.frame.ellipsoid().to_ecef(point))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn altitude_band() {
        assert!(AltitudeBand::new(Some(3000.0), Some(10_000.0)).contains_height(5000.0));
        assert!(!AltitudeBand::new(Some(3000.0), Some(10_000.0)).contains_height(2000.0));
        assert!(AltitudeBand::default().contains_height(999_999.0));
        assert!(AltitudeBand::below(2000.0).contains_height(1000.0));
        assert!(!AltitudeBand::above(2000.0).contains_height(1000.0));
        assert!(AltitudeBand::new(Some(5000.0), Some(12_000.0)).contains(Geodetic::new(1.0, 2.0, 8000.0)));
    }

    #[test]
    fn sphere() {
        let site = Geodetic::new(34.0, -118.0, 0.0);
        let s = Sphere::new(site, 10_000.0);
        assert!(s.contains(Geodetic::new(34.0, -118.0, 1000.0)));
        assert!(!s.contains(Geodetic::new(34.2, -118.0, 0.0)));
        assert_eq!(s.center(), site);
        assert_eq!(s.radius_m(), 10_000.0);
    }

    #[test]
    fn cone_directional_volume() {
        let up = Cone::new(Geodetic::new(0.0, 0.0, 0.0), 0.0, 90.0, 30.0, 20_000.0);
        assert!(up.contains(Geodetic::new(0.0, 0.0, 5000.0)));
        assert!(!up.contains(Geodetic::new(0.0, 0.0, 25_000.0)));
        assert!(!up.contains(Geodetic::new(0.2, 0.0, 1000.0)));
        assert!(up.contains(up.apex()));
        let east = Cone::new(Geodetic::new(0.0, 0.0, 1000.0), 90.0, 0.0, 20.0, 100_000.0);
        assert!(east.contains(Geodetic::new(0.0, 0.3, 1000.0)));
        assert!(!east.contains(Geodetic::new(0.3, 0.0, 1000.0)));
        assert_eq!(east.axis().azimuth_deg, 90.0);
        assert_eq!(east.half_angle_deg(), 20.0);
    }
}
