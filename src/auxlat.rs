//! Auxiliary latitudes. For now the conformal latitude χ, which the conformal map projections
//! (transverse Mercator, polar stereographic, and later Mercator and rhumb lines) are built on. Both
//! directions use exact formulas, not series, after Karney (2011) and GeographicLib's
//! `Math::taupf` and `Math::tauf`. They work with tangents, so they stay precise right up to the
//! poles.

#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;

/// The signed eccentricity: `e` for an oblate ellipsoid, `−|e|` for a prolate one (flattening < 0).
pub(crate) fn signed_eccentricity(e2: f64, f: f64) -> f64 {
    e2.abs().sqrt().copysign(f)
}

/// `e·atanh(e·x)`, continued to a prolate ellipsoid (`es < 0`) as `|e|·atan(|e|·x)` with the sign
/// flipped, so one formula serves both.
pub(crate) fn eatanhe(x: f64, es: f64) -> f64 {
    if es > 0.0 {
        es * (es * x).atanh()
    } else {
        -es * (es * x).atan()
    }
}

/// tan χ, the tangent of the conformal latitude, from tan φ.
pub(crate) fn taupf(tau: f64, es: f64) -> f64 {
    if tau.is_finite() {
        let tau1 = 1f64.hypot(tau);
        let sig = eatanhe(tau / tau1, es).sinh();
        1f64.hypot(sig) * tau - sig * tau1
    } else {
        tau // ±∞ stays ±∞; the formula would give NaN
    }
}

/// tan φ from tan χ: the reverse of [`taupf`], by Newton's method (one or two steps; at most five).
pub(crate) fn tauf(taup: f64, es: f64) -> f64 {
    const TOL: f64 = 1.490_116_119_384_765_6e-8 / 10.0; // √ε / 10
    const TAUMAX: f64 = 2.0 / 1.490_116_119_384_765_6e-8; // 2 / √ε
                                                          // 1 − e², with e² negative on a prolate ellipsoid. (GeographicLib uses 1 − es², which is the
                                                          // same when oblate but gives Newton the wrong slope when prolate, so it stops short.)
    let e2m = 1.0 - es * es.abs();
    // tan φ ≈ tan χ / (1 − e²) near the equator; for |tan χ| > 70 (|φ| > 89°) tan χ grows like
    // tan φ · exp(−e·atanh e), and that start also avoids overflow.
    let mut tau = if taup.abs() > 70.0 {
        taup * eatanhe(1.0, es).exp()
    } else {
        taup / e2m
    };
    if tau.is_nan() || tau.abs() >= TAUMAX {
        return tau; // ±∞, NaN, and |tan φ| so large that 1 + tan² φ = tan² φ
    }
    let stol = TOL * taup.abs().max(1.0);
    for _ in 0..5 {
        let taupa = taupf(tau, es);
        let dtau = (taup - taupa) * (1.0 + e2m * tau * tau) / (e2m * 1f64.hypot(tau) * 1f64.hypot(taupa));
        tau += dtau;
        if dtau.is_nan() || dtau.abs() < stol {
            break;
        }
    }
    tau
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::WGS84;

    /// The textbook conformal latitude, written independently:
    /// χ = 2·atan(tan(π/4 + φ/2)·((1 − e sin φ)/(1 + e sin φ))^(e/2)) − π/2.
    fn textbook_chi(phi: f64, e: f64) -> f64 {
        let s = phi.sin();
        2.0 * ((core::f64::consts::FRAC_PI_4 + phi / 2.0).tan() * ((1.0 - e * s) / (1.0 + e * s)).powf(e / 2.0)).atan()
            - core::f64::consts::FRAC_PI_2
    }

    #[test]
    fn conformal_latitude_matches_the_textbook_and_reverses() {
        let es = signed_eccentricity(WGS84.e2(), WGS84.f());
        for i in -89..=89 {
            let phi = (i as f64 + 0.37).to_radians();
            let chi = taupf(phi.tan(), es).atan();
            assert!(
                (chi - textbook_chi(phi, es)).abs() < 2e-15,
                "{i}: {chi} vs {}",
                textbook_chi(phi, es)
            );
        }
        // Back again, to a few units in the last place, from the equator to 1e-12° off the pole.
        for k in 0..2_000 {
            let lat = 90.0 * (1.0 - 10f64.powf(-12.0 * k as f64 / 2_000.0));
            for t in [
                lat.to_radians().tan(),
                -lat.to_radians().tan(),
                (lat / 90.0).to_radians().tan(),
            ] {
                let back = tauf(taupf(t, es), es);
                assert!(
                    (back - t).abs() <= 4.0 * f64::EPSILON * t.abs().max(1e-300),
                    "{t}: {back}"
                );
            }
        }
    }

    #[test]
    fn special_cases() {
        // A sphere: the conformal latitude is the latitude.
        assert_eq!(taupf(0.7, 0.0), 0.7);
        assert_eq!(tauf(0.7, 0.0), 0.7);
        // The poles and the equator.
        let es = signed_eccentricity(WGS84.e2(), WGS84.f());
        assert_eq!(taupf(f64::INFINITY, es), f64::INFINITY);
        assert_eq!(tauf(f64::NEG_INFINITY, es), f64::NEG_INFINITY);
        assert_eq!(taupf(0.0, es), 0.0);
        assert!(taupf(f64::NAN, es).is_nan() && tauf(f64::NAN, es).is_nan());
        // A prolate ellipsoid works too, and reverses.
        let es = signed_eccentricity(-0.01, -0.005);
        assert!(es < 0.0 && (tauf(taupf(1.3, es), es) - 1.3).abs() < 4e-16);
    }
}
