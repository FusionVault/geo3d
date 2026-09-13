//! First-order kinematics.

use crate::types::{Ecef, Vec3};

/// The closest point of approach of two objects in linear motion.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cpa {
    /// Seconds from now until the closest approach; negative when it is already past. Zero when the
    /// relative velocity is zero.
    pub time_s: f64,
    /// The separation at that moment, metres.
    pub distance_m: f64,
}

/// Closest point of approach of two objects — ECEF positions (metres) and velocities (metres per
/// second) — assuming constant velocity. A first-order conflict/collision estimate: compare
/// `distance_m` with a separation threshold (aircraft separation, vessel CPA, rendezvous).
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

    #[test]
    fn parallel_non_closing_is_now() {
        let c = closest_approach(
            Ecef::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Ecef::new(0.0, 3.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        );
        assert_eq!(c.time_s, 0.0);
        assert!((c.distance_m - 3.0).abs() < 1e-9);
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
}

// ── 0.2 additions: tracks ────────────────────────────────────────────────────────────────────────────

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

    /// A track inferred from two fixes `dt_s` seconds apart (the second fix is "now"). A
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

    /// The position `dt_s` seconds from now (negative for the past) — dead reckoning.
    #[inline]
    pub fn at(self, dt_s: f64) -> Ecef {
        self.position.offset(self.velocity * dt_s)
    }

    /// Closest point of approach with another track.
    #[inline]
    pub fn closest_approach(self, other: Track) -> Cpa {
        closest_approach(self.position, self.velocity, other.position, other.velocity)
    }

    /// Whether the two tracks come within `threshold_m` of each other at or after now, and when.
    /// Returns the CPA only if it is in the future (or now) and inside the threshold.
    pub fn conflict(self, other: Track, threshold_m: f64) -> Option<Cpa> {
        let cpa = self.closest_approach(other);
        (cpa.time_s >= 0.0 && cpa.distance_m <= threshold_m).then_some(cpa)
    }
}

#[cfg(test)]
mod tests_0_2 {
    use super::*;

    #[test]
    fn tracks_dead_reckon_and_detect_conflicts() {
        let a = Track::new(Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(a.at(10.0), Ecef::new(10.0, 0.0, 0.0));
        assert_eq!(a.at(-2.0), Ecef::new(-2.0, 0.0, 0.0));
        let b = Track::new(Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let cpa = a.closest_approach(b);
        assert!((cpa.time_s - 5.0).abs() < 1e-9 && (cpa.distance_m - 5.0).abs() < 1e-9);
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
}
