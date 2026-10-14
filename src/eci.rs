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
//! **Time.** The Earth's rotation follows UT1, not the UTC that clocks and Unix time use. The two
//! differ by DUT1 = UT1 − UTC, which is never more than 0.9 s. That's up to about 400 m at the
//! equator. [`gmst_from_unix`] treats UTC as UT1. For full accuracy, use
//! [`gmst_from_unix_with_dut1`], with DUT1 from the IERS Bulletin A. Polar motion (about 10 m at
//! the surface) is ignored.
//!
//! [`era_rad`] gives the Earth Rotation Angle, the modern (IAU 2000) measure of the Earth's turn,
//! for use with the CIO-based frames.
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

use core::f64::consts::TAU;

#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::{Ecef, Eci};

/// The J2000 epoch, 2000-01-01T12:00:00, as a Julian date.
const J2000_JD: f64 = 2_451_545.0;
/// The J2000 epoch as Unix seconds.
const J2000_UNIX: f64 = 946_728_000.0;

/// Julian date from Unix seconds since 1970-01-01T00:00:00Z (UTC in, so UT1 only to within 0.9 s).
///
/// A single `f64` Julian date resolves time to about 40 µs. The `*_from_unix` functions avoid it
/// and keep sub-microsecond resolution.
#[inline]
pub fn julian_date(unix_seconds: f64) -> f64 {
    2_440_587.5 + unix_seconds / 86_400.0
}

/// Greenwich mean sidereal time for a UT1 Julian date, radians in `[0, 2π)`. IAU 1982.
#[doc(alias = "teme")]
pub fn gmst_rad(julian_date_ut1: f64) -> f64 {
    gmst_from_days(julian_date_ut1 - J2000_JD)
}

/// Greenwich mean sidereal time for a Unix time, taking UTC as UT1 (up to 0.9 s out), radians in
/// `[0, 2π)`.
#[inline]
pub fn gmst_from_unix(unix_seconds: f64) -> f64 {
    gmst_from_unix_with_dut1(unix_seconds, 0.0)
}

/// Greenwich mean sidereal time for a Unix (UTC) time, corrected to UT1 by `dut1_s` = UT1 − UTC in
/// seconds (from the IERS Bulletin A; always within ±0.9 s), radians in `[0, 2π)`.
///
/// ```
/// use geo3d::eci;
/// let t = 1_800_000_000.0;
/// let utc_as_ut1 = eci::gmst_from_unix(t);
/// let corrected = eci::gmst_from_unix_with_dut1(t, 0.25);
/// // A quarter of a second of the Earth's turn: about 116 m at the equator.
/// assert!(((corrected - utc_as_ut1) * 6_378_137.0 - 116.3).abs() < 0.5);
/// ```
pub fn gmst_from_unix_with_dut1(unix_seconds: f64, dut1_s: f64) -> f64 {
    gmst_from_days(((unix_seconds - J2000_UNIX) + dut1_s) / 86_400.0)
}

/// The Earth Rotation Angle for a UT1 Julian date, radians in `[0, 2π)` (IAU 2000). It's the
/// CIO-based successor to sidereal time.
pub fn era_rad(julian_date_ut1: f64) -> f64 {
    era_from_days(julian_date_ut1 - J2000_JD)
}

/// The Earth Rotation Angle for a Unix (UTC) time, corrected to UT1 by `dut1_s` = UT1 − UTC in
/// seconds, radians in `[0, 2π)`.
pub fn era_from_unix_with_dut1(unix_seconds: f64, dut1_s: f64) -> f64 {
    era_from_days(((unix_seconds - J2000_UNIX) + dut1_s) / 86_400.0)
}

/// GMST (IAU 1982) for `d` UT1 days since J2000. The formula turns 360.985 647 366 29° a day, so
/// the 360° whole turns are taken from the fraction of the day before anything is rounded, and
/// the rest stays small.
fn gmst_from_days(d: f64) -> f64 {
    let t = d / 36_525.0; // Julian centuries
                          // 8 640 184.812 866 s per century beyond whole turns, and the T² and T³ terms, in degrees.
    let extra_per_day = 8_640_184.812_866 / 36_525.0 / 240.0;
    // 67 310.548 41 s is 280.460 618 375°. The often-quoted 280.460 618 37° is 8.7e-11 rad short.
    let deg =
        67_310.548_41 / 240.0 + 360.0 * (d - d.floor()) + extra_per_day * d + (0.093_104 - 6.2e-6 * t) * t * t / 240.0;
    wrap_turn(deg / 360.0)
}

/// ERA = 2π (0.779 057 273 264 + 1.002 737 811 911 354 48 Du), with the whole turns of `Du` taken
/// first, as SOFA does.
fn era_from_days(du: f64) -> f64 {
    wrap_turn((du - du.floor()) + 0.779_057_273_264 + 0.002_737_811_911_354_48 * du)
}

/// A number of turns as an angle in radians in `[0, 2π)`.
fn wrap_turn(turns: f64) -> f64 {
    let r = (turns - turns.floor()) * TAU;
    if r >= TAU {
        0.0
    } else {
        r
    }
}

/// Earth-fixed → inertial (TEME) at a sidereal angle in radians: rotate about `z` by `gmst_rad`.
#[doc(alias = "teme")]
#[inline]
pub fn ecef_to_eci(p: Ecef, gmst_rad: f64) -> Eci {
    let (s, c) = gmst_rad.sin_cos();
    Eci::new(p.x * c - p.y * s, p.x * s + p.y * c, p.z)
}

/// Inertial (TEME) → Earth-fixed at a sidereal angle in radians: rotate about `z` by `-gmst_rad`.
#[doc(alias = "teme")]
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

    /// The IAU 1982 formula as it was written before (seconds of a day from Julian centuries).
    fn gmst_old(jd: f64) -> f64 {
        let t = (jd - 2_451_545.0) / 36_525.0;
        let seconds =
            67_310.548_41 + (876_600.0 * 3600.0 + 8_640_184.812_866) * t + 0.093_104 * t * t - 6.2e-6 * t * t * t;
        (seconds.rem_euclid(86_400.0) / 86_400.0) * TAU
    }

    #[test]
    fn matches_sofa_reference_values() {
        // SOFA's own test values (t_sofa_c.c): iauGmst82(2400000.5, 53736.0) and
        // iauEra00(2400000.5, 54388.0).
        assert!(close(gmst_rad(2_400_000.5 + 53_736.0), 1.754_174_981_860_675, 1e-12));
        assert!(close(era_rad(2_400_000.5 + 54_388.0), 0.402_283_724_002_815_8, 1e-12));
        // ERA and GMST agree at J2000 to within the 0.0003° their constants differ by.
        assert!(close(era_rad(2_451_545.0), gmst_rad(2_451_545.0), 1e-5));
    }

    #[test]
    fn rearranged_gmst_matches_the_formula_and_keeps_microseconds() {
        for jd in [2_400_000.5, 2_451_545.0, 2_460_000.25, 2_470_000.75] {
            assert!(close(gmst_rad(jd), gmst_old(jd), 1e-9), "{jd}");
        }
        // A 10 µs step in 2027 changes GMST by 10 µs of the Earth's turn (7.29e-10 rad). A single
        // f64 Julian date (40 µs resolution) can't resolve that; the Unix-time path can.
        let t = 1_800_000_000.0;
        let step = gmst_from_unix(t + 1e-5) - gmst_from_unix(t);
        assert!(close(step, 7.292_115_9e-10, 0.05 * 7.29e-10), "{step}");
        // DUT1 shifts the angle by DUT1 seconds of turning.
        let shift = gmst_from_unix_with_dut1(t, 0.5) - gmst_from_unix(t);
        assert!(close(shift, 0.5 * 7.292_115_9e-5, 1e-9), "{shift}");
        assert!(close(
            era_from_unix_with_dut1(946_728_000.0, 0.0),
            era_rad(2_451_545.0),
            1e-15
        ));
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
