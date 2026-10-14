//! Angles in degrees, computed so that exact inputs stay exact: sine and cosine of whole multiples
//! of 90°, arctangents of the cardinal directions, and differences of longitudes. Shared by the
//! geodesics, the attitude and route tools, and the map projections.

#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;

/// Error-free sum: `s + t == u + v` exactly.
#[inline]
pub(crate) fn two_sum(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;
    (s, -((up - u) + (vpp - v)))
}

/// Tiny angles round to zero; the rest keep a resolution of 2⁻⁵⁷ (0.7 pm on the Earth), so
/// near-singular inputs like 1e-200° need no special handling.
pub(crate) fn ang_round(x: f64) -> f64 {
    const Z: f64 = 1.0 / 16.0;
    let mut y = x.abs();
    if y < Z {
        y = Z - (Z - y);
    }
    if x == 0.0 {
        0.0
    } else {
        y.copysign(x)
    }
}

/// An angle in degrees reduced to `(-180, 180]`.
pub(crate) fn ang_normalize(x: f64) -> f64 {
    if x > -180.0 && x <= 180.0 {
        return x; // the common case, without the cost of a remainder
    }
    let z = if x.is_finite() { x % 360.0 } else { f64::NAN };
    let z = if x == 0.0 { x } else { z };
    let y = if z < -180.0 {
        z + 360.0
    } else if z < 180.0 {
        z
    } else {
        z - 360.0
    };
    if y == -180.0 {
        180.0
    } else {
        y
    }
}

/// NaN for a latitude beyond ±90°.
pub(crate) fn lat_fix(x: f64) -> f64 {
    if x.abs() > 90.0 {
        f64::NAN
    } else {
        x
    }
}

/// `y − x` reduced to `[-180, 180]`, as an exact sum `(d, e)`.
pub(crate) fn ang_diff(x: f64, y: f64) -> (f64, f64) {
    let (d, t) = two_sum(ang_normalize(-x), ang_normalize(y));
    let d = ang_normalize(d);
    if d == 180.0 && t > 0.0 {
        two_sum(-180.0, t)
    } else {
        two_sum(d, t)
    }
}

/// Sine and cosine of an angle in degrees, exact at multiples of 90°: the angle is reduced to
/// `|r| ≤ 45°` before converting to radians.
pub(crate) fn sincosd(x: f64) -> (f64, f64) {
    // `%` is a slow library call, so it is only used for large angles; below 720° every step here
    // is exact (Sterbenz's lemma), so the reduction loses nothing. NaN and ±∞ give NaN.
    let r = if x.abs() < 720.0 { x } else { x % 360.0 };
    let mut q = (r / 90.0) as i64;
    let mut r = r - 90.0 * q as f64;
    if r > 45.0 {
        q += 1;
        r -= 90.0;
    } else if r < -45.0 {
        q -= 1;
        r += 90.0;
    }
    // ±45° and ±30° come out correctly rounded (as in GeographicLib): sin 45° = cos 45° = √½
    // exactly, and sin 30° = ½.
    let (s, c) = if 2.0 * r.abs() == 90.0 {
        (
            core::f64::consts::FRAC_1_SQRT_2.copysign(r),
            core::f64::consts::FRAC_1_SQRT_2,
        )
    } else if 3.0 * r.abs() == 90.0 {
        (0.5f64.copysign(r), 3f64.sqrt() / 2.0)
    } else {
        r.to_radians().sin_cos()
    };
    let (s, c) = match q & 3 {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        _ => (-c, s),
    };
    // No negative zeros except sin(−0).
    (if s == 0.0 { s.copysign(x) } else { s }, c + 0.0)
}

/// `atan2(y, x)` in degrees, in `(-180, 180]`, exact for the cardinal directions.
pub(crate) fn atan2d(mut y: f64, mut x: f64) -> f64 {
    let mut q = if y.abs() > x.abs() {
        core::mem::swap(&mut x, &mut y);
        2
    } else {
        0
    };
    if x < 0.0 {
        q += 1;
        x = -x;
    }
    let ang = y.atan2(x).to_degrees();
    match q {
        1 => (if y >= 0.0 { 180.0 } else { -180.0 }) - ang,
        2 => 90.0 - ang,
        3 => -90.0 + ang,
        _ => ang,
    }
}

/// The tangent of an angle in degrees. At ±90° it is ±1/ε² (a large finite number), not infinite,
/// so callers can keep working with it.
pub(crate) fn tand(x: f64) -> f64 {
    const OVERFLOW: f64 = 1.0 / (f64::EPSILON * f64::EPSILON);
    let (s, c) = sincosd(x);
    (s / c).clamp(-OVERFLOW, OVERFLOW)
}

/// The arctangent in degrees, exact at 0, ±1 and ±∞.
pub(crate) fn atand(x: f64) -> f64 {
    atan2d(x, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_where_they_should_be() {
        for (deg, s, c) in [
            (0.0, 0.0, 1.0),
            (90.0, 1.0, 0.0),
            (180.0, 0.0, -1.0),
            (-90.0, -1.0, 0.0),
            (270.0, -1.0, 0.0),
            (720.0, 0.0, 1.0),
        ] {
            assert_eq!(sincosd(deg), (s, c), "{deg}");
        }
        // Correctly rounded at 45° and 30°.
        let h = core::f64::consts::FRAC_1_SQRT_2;
        assert_eq!(sincosd(45.0), (h, h));
        assert_eq!(sincosd(-135.0), (-h, -h));
        assert_eq!(sincosd(30.0), (0.5, 3f64.sqrt() / 2.0));
        assert_eq!(sincosd(150.0).0, 0.5);
        assert_eq!(atand(1.0), 45.0);
        assert_eq!(atand(f64::INFINITY), 90.0);
        assert_eq!(tand(45.0), 1.0);
        assert_eq!(tand(0.0), 0.0);
        assert!(tand(90.0) > 1e30 && tand(-90.0) < -1e30);
        assert!(tand(f64::NAN).is_nan());
        assert!(lat_fix(90.5).is_nan());
        assert_eq!(ang_normalize(-180.0), 180.0);
        let (d, e) = ang_diff(-179.0, 179.0);
        assert_eq!(d + e, -2.0);
    }
}
