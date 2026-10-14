//! Earth-fixed ([`Ecef`]) ↔ inertial ([`Eci`]): undo the Earth's spin for satellite work.
//!
//! [`Eci`] shares the Earth's axis with [`Ecef`] but does not turn with it: it is ECEF rotated about
//! `z` by Greenwich mean sidereal time (IAU 1982), the angle the Earth has turned.
//!
//! **Which inertial frame this is.** That rotation is exactly the step between the Earth-fixed frame
//! and **TEME**, the frame SGP4 / two-line-element (TLE) propagators output, so SGP4 positions can
//! be fed straight in. It is **not** J2000 / GCRF: those also apply precession and nutation and
//! differ from TEME by about 0.4° today (hundreds of kilometres at geostationary distance).
//!
//! **Time.** [`gmst_from_unix`] treats UTC as UT1. The two differ by up to 0.9 s, which moves an
//! equatorial point by up to about 400 m. For better, pass a UT1 Julian date to [`gmst_rad`].
//! Polar motion (about 10 m at the surface) is ignored.
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

use crate::position::{Ecef, Eci};

/// Julian date from Unix seconds since 1970-01-01T00:00:00Z (UTC in, so UT1 only to within 0.9 s).
#[inline]
pub fn julian_date(unix_seconds: f64) -> f64 {
    2_440_587.5 + unix_seconds / 86_400.0
}

/// Greenwich mean sidereal time for a UT1 Julian date, radians in `[0, 2π)`. IAU 1982.
pub fn gmst_rad(julian_date_ut1: f64) -> f64 {
    let t = (julian_date_ut1 - 2_451_545.0) / 36_525.0; // Julian centuries from J2000
    let seconds = 67_310.548_41 + (876_600.0 * 3600.0 + 8_640_184.812_866) * t + 0.093_104 * t * t - 6.2e-6 * t * t * t;
    (seconds.rem_euclid(86_400.0) / 86_400.0) * core::f64::consts::TAU
}

/// Greenwich mean sidereal time for a Unix time (taking UTC as UT1), radians in `[0, 2π)`.
#[inline]
pub fn gmst_from_unix(unix_seconds: f64) -> f64 {
    gmst_rad(julian_date(unix_seconds))
}

/// Earth-fixed → inertial (TEME) at a sidereal angle in radians: rotate about `z` by `gmst_rad`.
#[inline]
pub fn ecef_to_eci(p: Ecef, gmst_rad: f64) -> Eci {
    let (s, c) = gmst_rad.sin_cos();
    Eci::new(p.x * c - p.y * s, p.x * s + p.y * c, p.z)
}

/// Inertial (TEME) → Earth-fixed at a sidereal angle in radians: rotate about `z` by `-gmst_rad`.
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
    use crate::test_util::close;

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
