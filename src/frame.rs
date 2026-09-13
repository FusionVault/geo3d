//! Topocentric (local tangent-plane) frames.

use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::types::{Aer, Ecef, Enu, Geodetic, Ned, Vec3};

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

// ── 0.2 additions: NED side, ECEF targets, courses ───────────────────────────────────────────────────

impl LocalFrame {
    /// A geodetic position as a NED offset from the origin.
    #[inline]
    pub fn ned_of(&self, g: Geodetic) -> Ned {
        self.enu_of(g).to_ned()
    }

    /// An ECEF position as a NED offset from the origin.
    #[inline]
    pub fn to_ned(&self, p: Ecef) -> Ned {
        self.to_enu(p).to_ned()
    }

    /// A NED offset from the origin as an ECEF position.
    #[inline]
    pub fn ecef_of_ned(&self, v: Ned) -> Ecef {
        self.to_ecef(v.to_enu())
    }

    /// A NED offset from the origin as a geodetic position.
    #[inline]
    pub fn geodetic_of_ned(&self, v: Ned) -> Geodetic {
        self.geodetic_of(v.to_enu())
    }

    /// Look angles from the origin to an ECEF target.
    #[inline]
    pub fn look_angles_ecef(&self, target: Ecef) -> Aer {
        self.to_enu(target).to_aer()
    }

    /// The ECEF position reached by pointing `aer` from the origin.
    #[inline]
    pub fn ecef_of_aer(&self, aer: Aer) -> Ecef {
        self.to_ecef(aer.to_enu())
    }

    /// An ECEF velocity (m/s) at the origin as course over ground, ground speed and climb rate.
    #[inline]
    pub fn course_of_velocity(&self, velocity: Vec3) -> Course {
        Course::from_enu_velocity(self.dir_to_enu(velocity))
    }

    /// A course (over ground) at the origin as an ECEF velocity vector (m/s).
    #[inline]
    pub fn velocity_of_course(&self, course: Course) -> Vec3 {
        self.dir_to_ecef(course.to_enu_velocity())
    }
}

/// Motion over the ground at a position: course (true track) in degrees from north clockwise,
/// ground speed and climb rate in m/s — the shape ADS-B, AIS and GNSS receivers report.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Course {
    /// Track over ground, degrees clockwise from north in `[0, 360)`.
    pub course_deg: f64,
    /// Horizontal speed, m/s.
    pub ground_speed_mps: f64,
    /// Vertical rate, m/s, positive up.
    pub climb_mps: f64,
}

impl Course {
    /// Build from a course in degrees and speeds in m/s.
    #[inline]
    pub const fn new(course_deg: f64, ground_speed_mps: f64, climb_mps: f64) -> Self {
        Course {
            course_deg,
            ground_speed_mps,
            climb_mps,
        }
    }

    /// From a velocity expressed in a local ENU frame (m/s). Zero horizontal speed gives course 0.
    pub fn from_enu_velocity(v: Enu) -> Course {
        let ground = (v.east * v.east + v.north * v.north).sqrt();
        let course_deg = if ground == 0.0 {
            0.0
        } else {
            v.east.atan2(v.north).to_degrees().rem_euclid(360.0)
        };
        Course {
            course_deg,
            ground_speed_mps: ground,
            climb_mps: v.up,
        }
    }

    /// As a velocity in a local ENU frame (m/s).
    #[inline]
    pub fn to_enu_velocity(self) -> Enu {
        let (s, c) = self.course_deg.to_radians().sin_cos();
        Enu::new(s * self.ground_speed_mps, c * self.ground_speed_mps, self.climb_mps)
    }

    /// Total speed through the air/water/space, m/s.
    #[inline]
    pub fn speed_mps(self) -> f64 {
        (self.ground_speed_mps * self.ground_speed_mps + self.climb_mps * self.climb_mps).sqrt()
    }
}

#[cfg(test)]
mod tests_0_2 {
    use super::*;
    use crate::WGS84;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn ned_side_mirrors_enu() {
        let f = LocalFrame::new(Geodetic::new(10.0, 20.0, 0.0));
        let g = Geodetic::new(10.01, 20.0, 500.0);
        let ned = f.ned_of(g);
        assert!(ned.north > 0.0 && close(ned.east, 0.0, 1.0) && close(ned.down, -500.0, 1.0));
        assert_eq!(f.to_ned(WGS84.to_ecef(g)), ned);
        let back = f.geodetic_of_ned(ned);
        assert!(close(back.lat_deg, g.lat_deg, 1e-9) && close(back.height_m, g.height_m, 1e-4));
        assert!(close(f.ecef_of_ned(ned).distance_to(WGS84.to_ecef(g)), 0.0, 1e-6));
    }

    #[test]
    fn ecef_targets_and_pointings() {
        let f = LocalFrame::new(Geodetic::new(0.0, 0.0, 0.0));
        let target = Geodetic::new(0.0, 0.0, 5000.0);
        let aer = f.look_angles_ecef(WGS84.to_ecef(target));
        assert!(close(aer.elevation_deg, 90.0, 1e-9) && close(aer.range_m, 5000.0, 1e-6));
        let p = f.ecef_of_aer(aer);
        assert!(close(p.distance_to(WGS84.to_ecef(target)), 0.0, 1e-6));
    }

    #[test]
    fn course_and_velocity_round_trip() {
        let f = LocalFrame::new(Geodetic::new(-33.9, 151.2, 0.0));
        let c = Course::new(45.0, 100.0, 5.0);
        let v = f.velocity_of_course(c);
        assert!(close(v.norm(), c.speed_mps(), 1e-9));
        let back = f.course_of_velocity(v);
        assert!(
            close(back.course_deg, 45.0, 1e-9)
                && close(back.ground_speed_mps, 100.0, 1e-9)
                && close(back.climb_mps, 5.0, 1e-9)
        );
        // Due north at 10 m/s: velocity is along the frame's north axis.
        let n = f.velocity_of_course(Course::new(0.0, 10.0, 0.0));
        let (_, north, _) = f.axes();
        assert!(close(n.angle_to(north), 0.0, 1e-12));
        // Hovering: course 0 by convention.
        assert_eq!(
            Course::from_enu_velocity(Enu::new(0.0, 0.0, 2.0)),
            Course::new(0.0, 0.0, 2.0)
        );
        assert!(close(Course::new(270.0, 3.0, 4.0).speed_mps(), 5.0, 1e-12));
        let w = Course::new(270.0, 10.0, 0.0).to_enu_velocity();
        assert!(close(w.east, -10.0, 1e-9) && close(w.north, 0.0, 1e-9));
    }
}
