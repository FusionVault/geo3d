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
