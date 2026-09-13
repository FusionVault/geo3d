//! Topocentric (local tangent-plane) frames.

use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::types::{Aer, Ecef, Enu, Geodetic, Vec3};

/// An **East-North-Up** frame anchored at a reference point on an ellipsoid.
///
/// Building a frame computes the origin's ECEF position and the rotation once; every conversion
/// through it is then a subtraction and three dot products, which is what makes "many targets
/// against one observer" loops cheap. Up is the ellipsoid normal at the origin.
///
/// ```
/// use geo3d::{Geodetic, LocalFrame};
///
/// let radar = LocalFrame::new(Geodetic::new(51.5, -0.12, 30.0));
/// let target = Geodetic::new(51.6, -0.12, 3_000.0); // ~11 km due north, 3 km up
/// let enu = radar.enu_of(target);
/// assert!(enu.east.abs() < 1.0 && enu.north > 11_000.0 && enu.up > 2_900.0);
///
/// let aer = radar.look_angles(target);
/// assert!(aer.azimuth_deg < 0.1 || aer.azimuth_deg > 359.9);
/// assert!(aer.elevation_deg > 14.0 && aer.elevation_deg < 16.0);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalFrame {
    ellipsoid: Ellipsoid,
    origin: Geodetic,
    origin_ecef: Ecef,
    e: Vec3,
    n: Vec3,
    u: Vec3,
}

impl LocalFrame {
    /// A frame at `origin` on [`WGS84`].
    #[inline]
    pub fn new(origin: Geodetic) -> Self {
        Self::on(WGS84, origin)
    }

    /// A frame at `origin` on the given ellipsoid.
    pub fn on(ellipsoid: Ellipsoid, origin: Geodetic) -> Self {
        let (slat, clat) = origin.lat_rad().sin_cos();
        let (slon, clon) = origin.lon_rad().sin_cos();
        LocalFrame {
            ellipsoid,
            origin,
            origin_ecef: ellipsoid.to_ecef(origin),
            e: Vec3::new(-slon, clon, 0.0),
            n: Vec3::new(-slat * clon, -slat * slon, clat),
            u: Vec3::new(clat * clon, clat * slon, slat),
        }
    }

    /// The reference point.
    #[inline]
    pub const fn origin(&self) -> Geodetic {
        self.origin
    }

    /// The reference point in ECEF.
    #[inline]
    pub const fn origin_ecef(&self) -> Ecef {
        self.origin_ecef
    }

    /// The ellipsoid the frame is built on.
    #[inline]
    pub const fn ellipsoid(&self) -> Ellipsoid {
        self.ellipsoid
    }

    /// The frame's unit axes `(east, north, up)` in ECEF.
    #[inline]
    pub const fn axes(&self) -> (Vec3, Vec3, Vec3) {
        (self.e, self.n, self.u)
    }

    /// Rotate an ECEF **direction** into the frame (no translation).
    #[inline]
    pub fn dir_to_enu(&self, d: Vec3) -> Enu {
        Enu::new(self.e.dot(d), self.n.dot(d), self.u.dot(d))
    }

    /// Rotate a frame **direction** into ECEF (no translation).
    #[inline]
    pub fn dir_to_ecef(&self, d: Enu) -> Vec3 {
        self.e * d.east + self.n * d.north + self.u * d.up
    }

    /// An ECEF position as an ENU offset from the origin.
    #[inline]
    pub fn to_enu(&self, p: Ecef) -> Enu {
        self.dir_to_enu(p.vec() - self.origin_ecef.vec())
    }

    /// An ENU offset from the origin as an ECEF position.
    #[inline]
    pub fn to_ecef(&self, v: Enu) -> Ecef {
        self.origin_ecef.offset(self.dir_to_ecef(v))
    }

    /// A geodetic position as an ENU offset from the origin.
    #[inline]
    pub fn enu_of(&self, g: Geodetic) -> Enu {
        self.to_enu(self.ellipsoid.to_ecef(g))
    }

    /// An ENU offset from the origin as a geodetic position.
    #[inline]
    pub fn geodetic_of(&self, v: Enu) -> Geodetic {
        self.ellipsoid.to_geodetic(self.to_ecef(v))
    }

    /// Azimuth, elevation and range from the origin to `target` — the look angles to point at it.
    #[inline]
    pub fn look_angles(&self, target: Geodetic) -> Aer {
        self.enu_of(target).to_aer()
    }

    /// The position reached by pointing `aer` from the origin — the inverse of [`look_angles`](Self::look_angles).
    #[inline]
    pub fn destination(&self, aer: Aer) -> Geodetic {
        self.geodetic_of(aer.to_enu())
    }

    /// Where a ray from the origin along `azimuth_deg` (0 = north, clockwise) and `elevation_deg`
    /// (0 = horizon, positive up) strikes the ellipsoid surface, or `None` if it points into space.
    /// The line-of-sight ground point of a camera, antenna or bearing.
    pub fn ground_strike(&self, azimuth_deg: f64, elevation_deg: f64) -> Option<Geodetic> {
        let dir = self.dir_to_ecef(Aer::new(azimuth_deg, elevation_deg, 1.0).direction_enu());
        self.ellipsoid
            .ray_intersect(self.origin_ecef, dir)
            .map(|p| self.ellipsoid.to_geodetic(p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn axes_are_orthonormal_and_right_handed() {
        let f = LocalFrame::new(Geodetic::new(-33.87, 151.21, 20.0));
        let (e, n, u) = f.axes();
        for v in [e, n, u] {
            assert!(close(v.norm(), 1.0, 1e-12));
        }
        assert!(close(e.dot(n), 0.0, 1e-12) && close(n.dot(u), 0.0, 1e-12) && close(u.dot(e), 0.0, 1e-12));
        let c = e.cross(n);
        assert!(close(c.x, u.x, 1e-12) && close(c.y, u.y, 1e-12) && close(c.z, u.z, 1e-12));
    }

    #[test]
    fn enu_orients_east_north_up() {
        let f = LocalFrame::new(Geodetic::new(10.0, 20.0, 0.0));
        let up = f.enu_of(Geodetic::new(10.0, 20.0, 500.0));
        assert!(close(up.east, 0.0, 1e-3) && close(up.north, 0.0, 1e-3) && close(up.up, 500.0, 1e-2));
        let north = f.enu_of(Geodetic::new(10.01, 20.0, 0.0));
        assert!(north.north > 0.0 && close(north.east, 0.0, 1.0));
        let east = f.enu_of(Geodetic::new(10.0, 20.01, 0.0));
        assert!(east.east > 0.0 && close(east.north, 0.0, 1.0));
    }

    #[test]
    fn enu_and_ecef_round_trip() {
        let f = LocalFrame::new(Geodetic::new(48.85, 2.35, 35.0));
        let g = Geodetic::new(50.0, 3.0, 12_000.0);
        let back = f.geodetic_of(f.enu_of(g));
        assert!(
            close(back.lat_deg, g.lat_deg, 1e-9)
                && close(back.lon_deg, g.lon_deg, 1e-9)
                && close(back.height_m, g.height_m, 1e-4)
        );
        let d = Vec3::new(0.3, -0.4, 0.5);
        let r = f.dir_to_ecef(f.dir_to_enu(d));
        assert!(close(r.x, d.x, 1e-12) && close(r.y, d.y, 1e-12) && close(r.z, d.z, 1e-12));
    }

    #[test]
    fn look_angles_and_destination_invert() {
        let f = LocalFrame::new(Geodetic::new(0.0, 0.0, 0.0));
        let (az, ..) = {
            let a = f.look_angles(Geodetic::new(0.1, 0.0, 0.0));
            (a.azimuth_deg, a.elevation_deg, a.range_m)
        };
        assert!(az < 1.0 || az > 359.0);
        let a = f.look_angles(Geodetic::new(0.0, 0.1, 0.0));
        assert!(close(a.azimuth_deg, 90.0, 1.0));
        let a = f.look_angles(Geodetic::new(0.0, 0.0, 5000.0));
        assert!(close(a.elevation_deg, 90.0, 1e-9) && close(a.range_m, 5000.0, 1e-6));
        for g in [Geodetic::new(0.5, -0.3, 8000.0), Geodetic::new(-0.2, 0.9, 100.0)] {
            let back = f.destination(f.look_angles(g));
            assert!(
                close(back.lat_deg, g.lat_deg, 1e-9)
                    && close(back.lon_deg, g.lon_deg, 1e-9)
                    && close(back.height_m, g.height_m, 1e-4)
            );
        }
    }

    #[test]
    fn ground_strike_hits_below_or_misses() {
        let f = LocalFrame::new(Geodetic::new(0.0, 0.0, 10_000.0));
        let g = f.ground_strike(0.0, -90.0).unwrap();
        assert!(close(g.lat_deg, 0.0, 1e-6) && close(g.lon_deg, 0.0, 1e-6) && close(g.height_m, 0.0, 1e-3));
        assert!(f.ground_strike(0.0, 90.0).is_none());
        let east = f.ground_strike(90.0, -5.0).unwrap();
        assert!(east.lon_deg > 0.0 && east.height_m.abs() < 1.0);
    }
}
