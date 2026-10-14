//! Moving objects: [`Course`] (heading, speed, climb), [`Track`] (position and velocity) and
//! [`closest_approach`] for near-miss checks. Motion is assumed constant — straight lines at
//! steady speed.

use crate::local::{wrap_360, Enu};
use crate::position::Ecef;
use crate::vector::Vec3;

/// Motion over the ground at a position: course in degrees clockwise from north, ground speed and
/// climb rate in m/s — what ADS-B, AIS and GNSS receivers report. Turn it into an ECEF velocity at
/// a position with [`LocalFrame::velocity_of_course`](crate::LocalFrame::velocity_of_course).
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

    /// From a velocity in a local ENU frame (m/s). Zero horizontal speed gives course 0.
    pub fn from_enu_velocity(v: Enu) -> Course {
        let ground = (v.east * v.east + v.north * v.north).sqrt();
        let course_deg = if ground == 0.0 {
            0.0
        } else {
            wrap_360(v.east.atan2(v.north).to_degrees())
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

    /// Total speed (horizontal and vertical together), m/s.
    #[inline]
    pub fn speed_mps(self) -> f64 {
        (self.ground_speed_mps * self.ground_speed_mps + self.climb_mps * self.climb_mps).sqrt()
    }
}

/// A moving object: ECEF position (metres) and velocity (m/s), assumed constant.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Track {
    /// Position now, metres.
    pub position: Ecef,
    /// Velocity, m/s.
    pub velocity: Vec3,
}

impl Track {
    /// A track from a position and velocity.
    #[inline]
    pub const fn new(position: Ecef, velocity: Vec3) -> Self {
        Track { position, velocity }
    }

    /// A track from two position fixes `dt_s` seconds apart (the second fix is "now"). A
    /// non-positive `dt_s` gives zero velocity.
    pub fn from_fixes(earlier: Ecef, now: Ecef, dt_s: f64) -> Self {
        let velocity = if dt_s > 0.0 {
            (now.vec() - earlier.vec()) / dt_s
        } else {
            Vec3::ZERO
        };
        Track {
            position: now,
            velocity,
        }
    }

    /// Where the object will be `dt_s` seconds from now (negative for the past) — dead reckoning.
    #[inline]
    pub fn at(self, dt_s: f64) -> Ecef {
        self.position.offset(self.velocity * dt_s)
    }

    /// When and how close this track gets to another.
    #[inline]
    pub fn closest_approach(self, other: Track) -> Cpa {
        closest_approach(self.position, self.velocity, other.position, other.velocity)
    }

    /// A near-miss check: `Some(cpa)` if the two tracks come within `threshold_m` of each other now
    /// or in the future, `None` otherwise.
    pub fn conflict(self, other: Track, threshold_m: f64) -> Option<Cpa> {
        let cpa = self.closest_approach(other);
        (cpa.time_s >= 0.0 && cpa.distance_m <= threshold_m).then_some(cpa)
    }
}

/// The closest point of approach of two moving objects: when, and how close.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cpa {
    /// Seconds from now until the closest approach; negative when it is already past, zero when
    /// the two move in parallel.
    pub time_s: f64,
    /// The separation at that moment, metres.
    pub distance_m: f64,
}

/// When and how close two objects get, from their ECEF positions (metres) and velocities (m/s),
/// assuming both keep going in a straight line at steady speed. Compare `distance_m` with a
/// separation threshold for aircraft, vessels or rendezvous — or use [`Track::conflict`].
///
/// ```
/// use geo3d::{closest_approach, Ecef, Vec3};
/// // Two objects converging along x, offset 5 m in y: closest in 5 s, 5 m apart.
/// let cpa = closest_approach(
///     Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0),
///     Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0),
/// );
/// assert!((cpa.time_s - 5.0).abs() < 1e-9 && (cpa.distance_m - 5.0).abs() < 1e-9);
/// ```
#[doc(alias = "cpa")]
pub fn closest_approach(p1: Ecef, v1: Vec3, p2: Ecef, v2: Vec3) -> Cpa {
    let prel = p1.vec() - p2.vec();
    let vrel = v1 - v2;
    let vv = vrel.dot(vrel);
    let time_s = if vv == 0.0 { 0.0 } else { -prel.dot(vrel) / vv };
    Cpa {
        time_s,
        distance_m: (prel + vrel * time_s).norm(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::close;

    #[test]
    fn parallel_non_closing_is_now() {
        let c = closest_approach(
            Ecef::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Ecef::new(0.0, 3.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        );
        assert_eq!(c.time_s, 0.0);
        assert!(close(c.distance_m, 3.0, 1e-9));
    }

    #[test]
    fn diverging_reports_a_past_time() {
        let c = closest_approach(
            Ecef::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Ecef::new(-10.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
        );
        assert!(c.time_s < 0.0 && c.distance_m < 1e-9);
    }

    #[test]
    fn tracks_dead_reckon_and_detect_conflicts() {
        let a = Track::new(Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(a.at(10.0), Ecef::new(10.0, 0.0, 0.0));
        assert_eq!(a.at(-2.0), Ecef::new(-2.0, 0.0, 0.0));
        let b = Track::new(Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let cpa = a.closest_approach(b);
        assert!(close(cpa.time_s, 5.0, 1e-9) && close(cpa.distance_m, 5.0, 1e-9));
        assert!(a.conflict(b, 6.0).is_some());
        assert!(a.conflict(b, 4.0).is_none());
        // A closest approach in the past is not a conflict.
        let past = Track::new(Ecef::new(-10.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        assert!(a.conflict(past, 100.0).is_none());
    }

    #[test]
    fn from_fixes_infers_velocity() {
        let t = Track::from_fixes(Ecef::new(0.0, 0.0, 0.0), Ecef::new(20.0, 0.0, 0.0), 10.0);
        assert_eq!(t.velocity, Vec3::new(2.0, 0.0, 0.0));
        assert_eq!(t.position, Ecef::new(20.0, 0.0, 0.0));
        assert_eq!(
            Track::from_fixes(Ecef::new(0.0, 0.0, 0.0), Ecef::new(20.0, 0.0, 0.0), 0.0).velocity,
            Vec3::ZERO
        );
    }

    #[test]
    fn course_conventions() {
        // Hovering: course 0 by convention.
        assert_eq!(
            Course::from_enu_velocity(Enu::new(0.0, 0.0, 2.0)),
            Course::new(0.0, 0.0, 2.0)
        );
        assert!(close(Course::new(270.0, 3.0, 4.0).speed_mps(), 5.0, 1e-12));
        let w = Course::new(270.0, 10.0, 0.0).to_enu_velocity();
        assert!(close(w.east, -10.0, 1e-9) && close(w.north, 0.0, 1e-9));
        let c = Course::from_enu_velocity(Course::new(123.0, 50.0, -3.0).to_enu_velocity());
        assert!(close(c.course_deg, 123.0, 1e-9) && close(c.ground_speed_mps, 50.0, 1e-9) && c.climb_mps == -3.0);
    }
}
