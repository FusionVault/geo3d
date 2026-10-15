//! [`LocalFrame`]: east, north and up from an observer.

use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::local::{Aer, Enu, Ned};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::motion::{Course, Track};
use crate::position::{Ecef, Geodetic};
use crate::vector::Vec3;

/// A local **East-North-Up** frame anchored at an observer — a flat "map" tangent to the Earth at
/// that point, with up along the ellipsoid normal.
///
/// Use it to find where things are relative to the observer ([`enu_of`](Self::enu_of),
/// [`ned_of`](Self::ned_of)), where to point ([`look_angles`](Self::look_angles)), and back again.
/// Building a frame does the trigonometry once; after that an ECEF target costs a few
/// multiplications (a geodetic one adds its conversion to ECEF), so keep one frame per observer
/// when converting many targets.
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

    /// The observer's position.
    #[inline]
    pub const fn origin(&self) -> Geodetic {
        self.origin
    }

    /// The observer's position in ECEF.
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

    // ── Positions in, local offsets out (and back) ──────────────────────────────────────────────

    /// Where a position is relative to the observer, as an ENU offset.
    #[doc(alias = "geodetic2enu")]
    #[inline]
    pub fn enu_of(&self, g: Geodetic) -> Enu {
        self.to_enu(self.ellipsoid.to_ecef(g))
    }

    /// Where a position is relative to the observer, as a NED offset.
    #[doc(alias = "geodetic2ned")]
    #[inline]
    pub fn ned_of(&self, g: Geodetic) -> Ned {
        self.enu_of(g).to_ned()
    }

    /// The position at an ENU offset from the observer.
    #[doc(alias = "enu2geodetic")]
    #[inline]
    pub fn geodetic_of(&self, v: Enu) -> Geodetic {
        self.ellipsoid.to_geodetic(self.to_ecef(v))
    }

    /// The position at a NED offset from the observer.
    #[doc(alias = "ned2geodetic")]
    #[inline]
    pub fn geodetic_of_ned(&self, v: Ned) -> Geodetic {
        self.geodetic_of(v.to_enu())
    }

    /// An ECEF position as an ENU offset from the observer.
    #[doc(alias = "ecef2enu")]
    #[inline]
    pub fn to_enu(&self, p: Ecef) -> Enu {
        self.dir_to_enu(p.vec() - self.origin_ecef.vec())
    }

    /// An ECEF position as a NED offset from the observer.
    #[doc(alias = "ecef2ned")]
    #[inline]
    pub fn to_ned(&self, p: Ecef) -> Ned {
        self.to_enu(p).to_ned()
    }

    /// The ECEF position at an ENU offset from the observer.
    #[doc(alias = "enu2ecef")]
    #[inline]
    pub fn to_ecef(&self, v: Enu) -> Ecef {
        self.origin_ecef.offset(self.dir_to_ecef(v))
    }

    /// The ECEF position at a NED offset from the observer.
    #[doc(alias = "ned2ecef")]
    #[inline]
    pub fn ecef_of_ned(&self, v: Ned) -> Ecef {
        self.to_ecef(v.to_enu())
    }

    /// Rotate an ECEF **direction** (or velocity) into the frame — rotation only, no translation.
    #[doc(alias = "ecef2enuv")]
    #[inline]
    pub fn dir_to_enu(&self, d: Vec3) -> Enu {
        Enu::new(self.e.dot(d), self.n.dot(d), self.u.dot(d))
    }

    /// Rotate a frame **direction** (or velocity) into ECEF — rotation only, no translation.
    #[doc(alias = "enu2ecefv")]
    #[inline]
    pub fn dir_to_ecef(&self, d: Enu) -> Vec3 {
        self.e * d.east + self.n * d.north + self.u * d.up
    }

    // ── Pointing ────────────────────────────────────────────────────────────────────────────────

    /// Where to point to see `target`: azimuth, elevation and range from the observer.
    #[doc(alias = "geodetic2aer")]
    #[doc(alias = "azimuth")]
    #[doc(alias = "elevation")]
    #[inline]
    pub fn look_angles(&self, target: Geodetic) -> Aer {
        self.enu_of(target).to_aer()
    }

    /// Where to point to see an ECEF `target`.
    #[doc(alias = "ecef2aer")]
    #[inline]
    pub fn look_angles_ecef(&self, target: Ecef) -> Aer {
        self.to_enu(target).to_aer()
    }

    /// The position a pointing reaches — the inverse of [`look_angles`](Self::look_angles).
    #[doc(alias = "aer2geodetic")]
    #[inline]
    pub fn destination(&self, aer: Aer) -> Geodetic {
        self.geodetic_of(aer.to_enu())
    }

    /// The ECEF position a pointing reaches.
    #[doc(alias = "aer2ecef")]
    #[inline]
    pub fn ecef_of_aer(&self, aer: Aer) -> Ecef {
        self.to_ecef(aer.to_enu())
    }

    /// Where the line of sight from the observer along `azimuth_deg` (0 = north, clockwise) and
    /// `elevation_deg` (0 = horizon, negative = down) first meets the ground (the ellipsoid
    /// surface) — the ground point of a camera, antenna or bearing.
    ///
    /// Returns `None` if the line points into space, the observer is below the surface, or any
    /// input is NaN or infinite. An observer exactly on the surface strikes its own position.
    #[doc(alias = "lookAtSpheroid")]
    pub fn ground_strike(&self, azimuth_deg: f64, elevation_deg: f64) -> Option<Geodetic> {
        let Geodetic {
            lat_deg,
            lon_deg,
            height_m: h,
        } = self.origin;
        let inputs = [lat_deg, lon_deg, h, azimuth_deg, elevation_deg];
        if inputs.iter().any(|x| !x.is_finite()) || h < 0.0 {
            return None;
        }
        // The height decides, not the sign of the near root: at the surface that root is 0 give or
        // take rounding, and taking the other one would land on the far side of the Earth.
        if h == 0.0 {
            return Some(self.origin);
        }
        let dir = self.dir_to_ecef(Aer::new(azimuth_deg, elevation_deg, 1.0).direction_enu());
        let (t1, t2) = self.ellipsoid.ray_roots(self.origin_ecef, dir)?;
        if t2 < 0.0 {
            return None; // the surface is behind the line of sight
        }
        // Above the surface both roots share a sign; a near root just below zero is rounding at a
        // height so small that the origin is at the surface, whose nearest point is its own.
        Some(self.ellipsoid.to_geodetic(self.origin_ecef.offset(dir * t1.max(0.0))))
    }

    // ── Velocities ──────────────────────────────────────────────────────────────────────────────

    /// An ECEF velocity (m/s) at the observer as course over ground, ground speed and climb rate.
    #[inline]
    pub fn course_of_velocity(&self, velocity: Vec3) -> Course {
        Course::from_enu_velocity(self.dir_to_enu(velocity))
    }

    /// A course over ground at the observer as an ECEF velocity (m/s).
    #[inline]
    pub fn velocity_of_course(&self, course: Course) -> Vec3 {
        self.dir_to_ecef(course.to_enu_velocity())
    }

    /// How fast `target` is moving away from the observer, who stands still on the ground, m/s:
    /// negative while it comes closer. See [`Track::range_rate_mps`] for the Doppler shift.
    ///
    /// ```
    /// use geo3d::{geodetic_to_ecef, Course, Geodetic, LocalFrame, Track};
    /// let radar = LocalFrame::new(Geodetic::new(-33.95, 151.18, 10.0));
    /// let plane = Geodetic::new(-33.80, 151.18, 3_000.0); // 17 km north, 3 km up
    /// // Flying level, due south towards the radar, at 120 m/s.
    /// let velocity = radar.velocity_of_course(Course::new(180.0, 120.0, 0.0));
    /// let rate = radar.range_rate_mps(Track::new(geodetic_to_ecef(plane), velocity));
    /// // Closing at 118 m/s: the line of sight is 10° up, so not all of the speed is towards it.
    /// assert!((rate + 118.1).abs() < 0.1);
    /// ```
    #[doc(alias = "doppler")]
    pub fn range_rate_mps(&self, target: Track) -> f64 {
        Track::new(self.origin_ecef, Vec3::ZERO).range_rate_mps(target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::close;

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
    fn look_angles_and_destination_invert() {
        let f = LocalFrame::new(Geodetic::new(0.0, 0.0, 0.0));
        let az = f.look_angles(Geodetic::new(0.1, 0.0, 0.0)).azimuth_deg;
        assert!(!(1.0..=359.0).contains(&az));
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
    fn ecef_targets_and_pointings() {
        let f = LocalFrame::new(Geodetic::new(0.0, 0.0, 0.0));
        let target = Geodetic::new(0.0, 0.0, 5000.0);
        let aer = f.look_angles_ecef(WGS84.to_ecef(target));
        assert!(close(aer.elevation_deg, 90.0, 1e-9) && close(aer.range_m, 5000.0, 1e-6));
        let p = f.ecef_of_aer(aer);
        assert!(close(p.distance_to(WGS84.to_ecef(target)), 0.0, 1e-6));
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

    #[test]
    fn ground_strike_takes_the_nearest_root_ahead() {
        // On the surface the observer is itself the nearest point of the surface ahead (t = 0),
        // whichever way it points, never the far side of the Earth.
        let on = Geodetic::new(-33.9, 151.2, 0.0);
        for (az, el) in [(0.0, -10.0), (200.0, -30.0), (90.0, -90.0), (0.0, 45.0)] {
            assert_eq!(LocalFrame::new(on).ground_strike(az, el), Some(on), "az {az}, el {el}");
        }
        // Below the surface no line of sight reaches it.
        let below = LocalFrame::new(Geodetic::new(0.0, 0.0, -100.0));
        for el in [10.0, -10.0, -90.0] {
            assert!(below.ground_strike(0.0, el).is_none(), "el {el}");
        }
        assert!(crate::ground_strike(Geodetic::new(0.0, 0.0, -0.5), 0.0, -90.0).is_none());
        // Above it, the near root: 45° down from 1 km lands about 1 km north.
        let g = LocalFrame::new(Geodetic::new(0.0, 0.0, 1_000.0))
            .ground_strike(0.0, -45.0)
            .unwrap();
        assert!(
            g.lat_deg > 0.0089 && g.lat_deg < 0.0091 && close(g.lon_deg, 0.0, 1e-9),
            "{g:?}"
        );
        // A non-finite height or angle has no line of sight.
        let f = LocalFrame::new(Geodetic::new(0.0, 0.0, 100.0));
        assert!(f.ground_strike(f64::NAN, -90.0).is_none() && f.ground_strike(0.0, f64::INFINITY).is_none());
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
    }
}
