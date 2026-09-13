//! Earth-Centred Inertial ↔ Earth-Fixed, via Greenwich mean sidereal time.
//!
//! ECI here is the frame that shares Earth's axis with ECEF but does not rotate with it: ECEF rotated
//! about Z by the sidereal angle. This is the simple model (IAU 1982 GMST, no precession, nutation
//! or polar motion) — the one satellite-tracking and ground-station code commonly uses for
//! ephemeris-to-ground work at the tens-of-metres level.
//!
//! ```
//! use geo3d::{eci, Ecef};
//! let jd = eci::julian_date(946_728_000.0); // 2000-01-01T12:00:00Z, the J2000 epoch
//! assert!((jd - 2_451_545.0).abs() < 1e-9);
//! let gmst = eci::gmst_rad(jd);
//! assert!((gmst.to_degrees() - 280.46061837).abs() < 1e-6);
//! let p = Ecef::new(7_000_000.0, 0.0, 0.0);
//! let back = eci::eci_to_ecef(eci::ecef_to_eci(p, gmst), gmst);
//! assert!((back.x - p.x).abs() < 1e-6);
//! ```

use crate::types::{Ecef, Eci};

/// Julian date (UT1 ≈ UTC) from Unix seconds since 1970-01-01T00:00:00Z.
#[inline]
pub fn julian_date(unix_seconds: f64) -> f64 {
    2_440_587.5 + unix_seconds / 86_400.0
}

/// Greenwich mean sidereal time for a Julian date (UT1), radians in `[0, 2π)`. IAU 1982.
pub fn gmst_rad(julian_date_ut1: f64) -> f64 {
    let t = (julian_date_ut1 - 2_451_545.0) / 36_525.0; // Julian centuries from J2000
    let seconds = 67_310.548_41 + (876_600.0 * 3600.0 + 8_640_184.812_866) * t + 0.093_104 * t * t - 6.2e-6 * t * t * t;
    (seconds.rem_euclid(86_400.0) / 86_400.0) * core::f64::consts::TAU
}

/// Greenwich mean sidereal time for a Unix time, radians.
#[inline]
pub fn gmst_from_unix(unix_seconds: f64) -> f64 {
    gmst_rad(julian_date(unix_seconds))
}

/// Earth-fixed → inertial: rotate about Z by `gmst_rad`.
#[inline]
pub fn ecef_to_eci(p: Ecef, gmst_rad: f64) -> Eci {
    let (s, c) = gmst_rad.sin_cos();
    Eci::new(p.x * c - p.y * s, p.x * s + p.y * c, p.z)
}

/// Inertial → Earth-fixed: rotate about Z by `-gmst_rad`.
#[inline]
pub fn eci_to_ecef(p: Eci, gmst_rad: f64) -> Ecef {
    let (s, c) = gmst_rad.sin_cos();
    Ecef::new(p.x * c + p.y * s, -p.x * s + p.y * c, p.z)
}

impl Eci {
    /// Earth-fixed position at the given sidereal angle.
    #[inline]
    pub fn to_ecef(self, gmst_rad: f64) -> Ecef {
        eci_to_ecef(self, gmst_rad)
    }
}

impl Ecef {
    /// Inertial position at the given sidereal angle.
    #[inline]
    pub fn to_eci(self, gmst_rad: f64) -> Eci {
        ecef_to_eci(self, gmst_rad)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn julian_dates() {
        assert_eq!(julian_date(0.0), 2_440_587.5);
        assert!(close(julian_date(946_728_000.0), 2_451_545.0, 1e-9));
    }

    #[test]
    fn gmst_matches_reference_values() {
        // J2000: 18h 41m 50.548s = 280.46061837°.
        assert!(close(gmst_rad(2_451_545.0).to_degrees(), 280.460_618_37, 1e-6));
        // One sidereal day later (~23h56m04.09s) GMST wraps back to the same angle.
        let sidereal_day = 86_164.090_5 / 86_400.0;
        assert!(close(
            gmst_rad(2_451_545.0 + sidereal_day).to_degrees(),
            280.460_618_37,
            1e-3
        ));
        // Always in [0, 2π).
        for jd in [2_400_000.5, 2_451_545.0, 2_470_000.25] {
            let g = gmst_rad(jd);
            assert!((0.0..core::f64::consts::TAU).contains(&g));
        }
        assert!(close(gmst_from_unix(946_728_000.0), gmst_rad(2_451_545.0), 1e-15));
    }

    #[test]
    fn rotation_conventions() {
        let p = Ecef::new(1.0, 0.0, 5.0);
        // At zero sidereal angle the frames coincide.
        assert_eq!(ecef_to_eci(p, 0.0), Eci::new(1.0, 0.0, 5.0));
        // At 90° the fixed x-axis has rotated onto the inertial y-axis.
        let q = ecef_to_eci(p, core::f64::consts::FRAC_PI_2);
        assert!(close(q.x, 0.0, 1e-15) && close(q.y, 1.0, 1e-15) && q.z == 5.0);
        // Round trip through the methods.
        let g = 1.234;
        let back = p.to_eci(g).to_ecef(g);
        assert!(close(back.x, p.x, 1e-15) && close(back.y, p.y, 1e-15) && back.z == p.z);
    }
}
