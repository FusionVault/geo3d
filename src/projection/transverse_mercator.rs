//! The transverse Mercator projection, by Krüger's series to eighth order in the third flattening
//! (Karney 2011, *Transverse Mercator with an accuracy of a few nanometers*,
//! <https://doi.org/10.1007/s00190-011-0445-3>), ported from GeographicLib.

use core::f64::consts::{FRAC_PI_2, PI};

use super::MapPoint;
use crate::angle::{ang_diff, ang_normalize, atan2d, atand, lat_fix, sincosd};
use crate::auxlat::{eatanhe, signed_eccentricity, tauf, taupf};
use crate::ellipsoid::{Ellipsoid, WGS84};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;

/// Series order: 8, two beyond GeographicLib's default, so the truncation error (n⁹ ≈ 10⁻²⁵) is far
/// below rounding.
const ORDER: usize = 8;

/// `b1·(1 + n)` as a polynomial in n², over the last entry (GeographicLib, order 8).
const B1_COEFF: [f64; 6] = [25.0, 64.0, 256.0, 4096.0, 16384.0, 16384.0];

/// Krüger's α (forward) and β (reverse) coefficients: for each order j, `α_j / nʲ` as a polynomial
/// in n, over its last entry (GeographicLib, order 8).
#[rustfmt::skip]
const ALP_COEFF: [f64; 44] = [
    -75900428.0, 37884525.0, 42422016.0, -89611200.0, 46287360.0, 63504000.0, -135475200.0,
    101606400.0, 203212800.0,
    148003883.0, 83274912.0, -178508970.0, 77690880.0, 67374720.0, -104509440.0, 47174400.0,
    174182400.0,
    318729724.0, -738126169.0, 294981280.0, 178924680.0, -234938880.0, 81164160.0, 319334400.0,
    -40176129013.0, 14967552000.0, 6971354016.0, -8165836800.0, 2355138720.0, 7664025600.0,
    10421654396.0, 3997835751.0, -4266773472.0, 1072709352.0, 2490808320.0,
    175214326799.0, -171950693600.0, 38652967262.0, 58118860800.0,
    -67039739596.0, 13700311101.0, 12454041600.0,
    1424729850961.0, 743921418240.0,
];
#[rustfmt::skip]
const BET_COEFF: [f64; 44] = [
    31777436.0, -37845269.0, 43097152.0, -42865200.0, -752640.0, 104428800.0, -180633600.0,
    135475200.0, 270950400.0,
    24749483.0, 14930208.0, -100683990.0, 152616960.0, -105719040.0, 23224320.0, 7257600.0,
    348364800.0,
    -232468668.0, 101880889.0, 39205760.0, -29795040.0, -28131840.0, 22619520.0, 638668800.0,
    324154477.0, 1433121792.0, -876745056.0, -167270400.0, 208945440.0, 7664025600.0,
    457888660.0, -312227409.0, -67920528.0, 70779852.0, 2490808320.0,
    -19841813847.0, -3665348512.0, 3758062126.0, 116237721600.0,
    -1989295244.0, 1979471673.0, 49816166400.0,
    191773887257.0, 3719607091200.0,
];

/// The polynomial `p[o]·xᵐ + … + p[o+m]`, by Horner's rule; `const`, so the coefficients for WGS84
/// are worked out at compile time.
const fn polyval(m: usize, p: &[f64], o: usize, x: f64) -> f64 {
    let (mut y, mut i) = (p[o], 1);
    while i <= m {
        y = y * x + p[o + i];
        i += 1;
    }
    y
}

/// A complex number as `(re, im)`.
type C = (f64, f64);

#[inline]
fn cmul(a: C, b: C) -> C {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

/// `(sinh x, cosh x)` from a single exponential, both to full relative precision (`expm1` keeps
/// sinh precise near 0).
#[inline]
fn sinh_cosh(x: f64) -> (f64, f64) {
    let em = x.exp_m1();
    let ep = em + 1.0;
    (em * (1.0 + 1.0 / ep) / 2.0, (ep + 1.0 / ep) / 2.0)
}

/// The transverse Mercator projection on an ellipsoid, with scale `k0` on the central meridian.
/// It's conformal (shapes and angles are kept locally) and accurate to about 5 nm within 3 900 km of
/// the central meridian. UTM is this projection with `k0 = 0.9996`, as
/// [`TransverseMercator::UTM`]; many national grids (British National Grid, Gauss–Krüger zones)
/// are too.
///
/// It gives plain x and y from the central meridian and the equator, without the false easting
/// and northing a grid adds.
///
/// ```
/// use geo3d::{Geodetic, TransverseMercator};
/// let p = TransverseMercator::UTM.forward(45.0, Geodetic::new(33.3, 44.4, 0.0));
/// assert!((p.x_m + 55_859.455_081_574).abs() < 1e-8);   // 55.9 km west of 45° E
/// let back = TransverseMercator::UTM.reverse(45.0, p.x_m, p.y_m);
/// assert!((back.lat_deg - 33.3).abs() < 1e-14);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransverseMercator {
    k0: f64,
    f: f64,
    e2: f64,
    e2m: f64,
    /// The rectifying radius: the meridian's length is 2π times this.
    a1: f64,
    b1: f64,
    alp: [f64; ORDER + 1],
    bet: [f64; ORDER + 1],
}

impl TransverseMercator {
    /// UTM's projection: WGS84 with `k0 = 0.9996`. Built at compile time.
    pub const UTM: TransverseMercator = TransverseMercator::new(WGS84, 0.9996);

    /// The projection on `ellipsoid` with scale `k0` on the central meridian (positive).
    pub const fn new(ellipsoid: Ellipsoid, k0: f64) -> Self {
        let (a, f) = (ellipsoid.a(), ellipsoid.f());
        let n = f / (2.0 - f);
        let m = ORDER / 2;
        let b1 = polyval(m, &B1_COEFF, 0, n * n) / (B1_COEFF[m + 1] * (1.0 + n));
        let (mut alp, mut bet) = ([0.0; ORDER + 1], [0.0; ORDER + 1]);
        let (mut o, mut d, mut l) = (0, n, 1);
        while l <= ORDER {
            let m = ORDER - l;
            alp[l] = d * polyval(m, &ALP_COEFF, o, n) / ALP_COEFF[o + m + 1];
            bet[l] = d * polyval(m, &BET_COEFF, o, n) / BET_COEFF[o + m + 1];
            o += m + 2;
            d *= n;
            l += 1;
        }
        let e2 = f * (2.0 - f);
        TransverseMercator {
            k0,
            f,
            e2,
            e2m: 1.0 - e2,
            a1: b1 * a,
            b1,
            alp,
            bet,
        }
    }

    /// The scale on the central meridian.
    #[inline]
    pub const fn k0(&self) -> f64 {
        self.k0
    }

    fn es(&self) -> f64 {
        signed_eccentricity(self.e2, self.f)
    }

    /// The scale factor at a pole, `√(1 − e²)·exp(e·atanh e)`.
    fn c(&self) -> f64 {
        self.e2m.sqrt() * eatanhe(1.0, self.es()).exp()
    }

    /// Project `g` (height ignored) with the central meridian at `central_meridian_deg`.
    pub fn forward(&self, central_meridian_deg: f64, g: Geodetic) -> MapPoint {
        self.project::<true>(central_meridian_deg, g)
    }

    /// [`forward`](Self::forward), with the convergence and scale left as NaN unless `SCALE`: they
    /// cost a second series and several arctangents that a grid conversion doesn't need.
    #[inline]
    pub(crate) fn project<const SCALE: bool>(&self, central_meridian_deg: f64, g: Geodetic) -> MapPoint {
        let lat = lat_fix(g.lat_deg);
        let lon = ang_diff(central_meridian_deg, g.lon_deg).0;
        // Work in the first quadrant (lat ≥ 0, 0 ≤ lon ≤ 90) and restore the signs at the end.
        let mut latsign = if lat.is_sign_negative() { -1.0 } else { 1.0 };
        let lonsign = if lon.is_sign_negative() { -1.0 } else { 1.0 };
        let (lat, mut lon) = (lat * latsign, lon * lonsign);
        let backside = lon > 90.0;
        if backside {
            if lat == 0.0 {
                latsign = -1.0;
            }
            lon = 180.0 - lon;
        }
        let (sphi, cphi) = sincosd(lat);
        let (slam, clam) = sincosd(lon);
        // ξ' and η': the Gauss–Schreiber coordinates, from the conformal latitude.
        let (mut xip, mut etap, mut gamma, mut k) = (FRAC_PI_2, 0.0, f64::NAN, f64::NAN);
        if lat != 90.0 {
            let tau = sphi / cphi;
            let taup = taupf(tau, self.es());
            xip = taup.atan2(clam);
            etap = (slam / taup.hypot(clam)).asinh();
            if SCALE {
                gamma = atan2d(slam * taup, clam * 1f64.hypot(taup));
                // A form whose rounding errors cancel (Karney 2011).
                k = (self.e2m + self.e2 * cphi * cphi).sqrt() * 1f64.hypot(tau) / taup.hypot(clam);
            }
        } else if SCALE {
            (gamma, k) = (lon, self.c());
        }
        // ζ = ζ' + Σ αⱼ sin(2jζ'), and its derivative, by Clenshaw summation in complex numbers.
        let (s0, c0) = (2.0 * xip).sin_cos();
        let (sh0, ch0) = sinh_cosh(2.0 * etap);
        let a: C = (2.0 * c0 * ch0, -2.0 * s0 * sh0); // 2 cos 2ζ'
        let (mut y0, mut y1, mut z0, mut z1): (C, C, C, C) = ((0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0));
        let mut j = ORDER; // even
        while j > 0 {
            let ay = cmul(a, y0);
            y1 = (ay.0 - y1.0 + self.alp[j], ay.1 - y1.1);
            if SCALE {
                let az = cmul(a, z0);
                z1 = (az.0 - z1.0 + 2.0 * j as f64 * self.alp[j], az.1 - z1.1);
            }
            j -= 1;
            let ay = cmul(a, y1);
            y0 = (ay.0 - y0.0 + self.alp[j], ay.1 - y0.1);
            if SCALE {
                let az = cmul(a, z1);
                z0 = (az.0 - z0.0 + 2.0 * j as f64 * self.alp[j], az.1 - z0.1);
            }
            j -= 1;
        }
        let ay = cmul((s0 * ch0, c0 * sh0), y0); // sin 2ζ' · y0
        let (xi, eta) = (xip + ay.0, etap + ay.1);
        let y = self.a1 * self.k0 * (if backside { PI - xi } else { xi }) * latsign;
        let x = self.a1 * self.k0 * eta * lonsign;
        if SCALE {
            let az = cmul((a.0 / 2.0, a.1 / 2.0), z0); // cos 2ζ' · z0
            let z1 = (1.0 - z1.0 + az.0, -z1.1 + az.1);
            gamma -= atan2d(z1.1, z1.0);
            k *= self.b1 * z1.0.hypot(z1.1) * self.k0;
            if backside {
                gamma = 180.0 - gamma;
            }
            gamma = ang_normalize(gamma * latsign * lonsign);
        }
        MapPoint {
            x_m: x,
            y_m: y,
            convergence_deg: gamma,
            scale: k,
        }
    }

    /// The position (height 0) at `x_m`, `y_m` from the central meridian at
    /// `central_meridian_deg` and the equator: the reverse of [`forward`](Self::forward).
    pub fn reverse(&self, central_meridian_deg: f64, x_m: f64, y_m: f64) -> Geodetic {
        self.unproject::<false>(central_meridian_deg, x_m, y_m).0
    }

    /// [`reverse`](Self::reverse), with the grid convergence (degrees) and scale there.
    #[cfg(test)]
    pub(crate) fn reverse_full(&self, lon0: f64, x_m: f64, y_m: f64) -> (Geodetic, f64, f64) {
        self.unproject::<true>(lon0, x_m, y_m)
    }

    /// The reverse projection; the convergence and scale are NaN unless `SCALE`.
    #[inline]
    fn unproject<const SCALE: bool>(&self, lon0: f64, x_m: f64, y_m: f64) -> (Geodetic, f64, f64) {
        let (xi, eta) = (y_m / (self.a1 * self.k0), x_m / (self.a1 * self.k0));
        let xisign = if xi.is_sign_negative() { -1.0 } else { 1.0 };
        let etasign = if eta.is_sign_negative() { -1.0 } else { 1.0 };
        let (mut xi, eta) = (xi * xisign, eta * etasign);
        let backside = xi > FRAC_PI_2;
        if backside {
            xi = PI - xi;
        }
        // ζ' = ζ − Σ βⱼ sin(2jζ), and its derivative.
        let (s0, c0) = (2.0 * xi).sin_cos();
        let (sh0, ch0) = sinh_cosh(2.0 * eta);
        let a: C = (2.0 * c0 * ch0, -2.0 * s0 * sh0);
        let (mut y0, mut y1, mut z0, mut z1): (C, C, C, C) = ((0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0));
        let mut j = ORDER;
        while j > 0 {
            let ay = cmul(a, y0);
            y1 = (ay.0 - y1.0 - self.bet[j], ay.1 - y1.1);
            if SCALE {
                let az = cmul(a, z0);
                z1 = (az.0 - z1.0 - 2.0 * j as f64 * self.bet[j], az.1 - z1.1);
            }
            j -= 1;
            let ay = cmul(a, y1);
            y0 = (ay.0 - y0.0 - self.bet[j], ay.1 - y0.1);
            if SCALE {
                let az = cmul(a, z1);
                z0 = (az.0 - z0.0 - 2.0 * j as f64 * self.bet[j], az.1 - z0.1);
            }
            j -= 1;
        }
        let ay = cmul((s0 * ch0, c0 * sh0), y0);
        let (xip, etap) = (xi + ay.0, eta + ay.1);
        let (mut gamma, mut k) = (f64::NAN, f64::NAN);
        if SCALE {
            let az = cmul((a.0 / 2.0, a.1 / 2.0), z0);
            let z1 = (1.0 - z1.0 + az.0, -z1.1 + az.1);
            gamma = atan2d(z1.1, z1.0);
            k = self.b1 / z1.0.hypot(z1.1);
        }
        let (s, ch) = sinh_cosh(etap);
        let (sxip, cxip) = xip.sin_cos();
        let c = cxip.max(0.0); // cos(π/2) can come out a hair negative
        let r = s.hypot(c);
        let (mut lat, mut lon);
        if r != 0.0 {
            lon = atan2d(s, c);
            let tau = tauf(sxip / r, self.es());
            lat = atand(tau);
            if SCALE {
                gamma += atan2d(sxip * (s / ch), c);
                k *= (self.e2m + self.e2 / (1.0 + tau * tau)).sqrt() * 1f64.hypot(tau) * r;
            }
        } else {
            lat = 90.0;
            lon = 0.0;
            if SCALE {
                k *= self.c();
            }
        }
        lat *= xisign;
        if backside {
            lon = 180.0 - lon;
        }
        lon *= etasign;
        lon = ang_normalize(lon + lon0);
        if SCALE {
            if backside {
                gamma = 180.0 - gamma;
            }
            gamma = ang_normalize(gamma * xisign * etasign);
            k *= self.k0;
        }
        (Geodetic::new(lat, lon, 0.0), gamma, k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apart(a: Geodetic, b: Geodetic) -> f64 {
        WGS84.to_ecef(a).distance_to(WGS84.to_ecef(b))
    }

    #[test]
    fn symmetries_poles_and_the_central_meridian() {
        let tm = TransverseMercator::UTM;
        let (e, w) = (
            tm.forward(9.0, Geodetic::new(40.0, 11.0, 0.0)),
            tm.forward(9.0, Geodetic::new(40.0, 7.0, 0.0)),
        );
        assert_eq!((e.x_m, e.y_m, e.scale), (-w.x_m, w.y_m, w.scale));
        assert_eq!(e.convergence_deg, -w.convergence_deg);
        let s = tm.forward(9.0, Geodetic::new(-40.0, 11.0, 0.0));
        assert_eq!((s.x_m, s.y_m), (e.x_m, -e.y_m));
        // On the central meridian: x = 0, convergence 0, scale exactly k0 at the equator.
        let c = tm.forward(9.0, Geodetic::new(0.0, 9.0, 0.0));
        assert_eq!((c.x_m, c.y_m, c.convergence_deg), (0.0, 0.0, 0.0));
        assert!((c.scale - 0.9996).abs() < 2e-16, "{c:?}");
        // The pole: a quarter of the meridian (10 001 965.729 m) times k0, from any longitude.
        for lon in [0.0, 37.0, 179.0] {
            let p = tm.forward(9.0, Geodetic::new(90.0, lon, 0.0));
            assert!(
                p.x_m == 0.0 && (p.y_m - 0.9996 * 10_001_965.729_312_7).abs() < 1e-6,
                "{p:?}"
            );
        }
        let pole = tm.reverse(9.0, 0.0, 0.9996 * 10_001_965.729_312_7);
        assert!((pole.lat_deg - 90.0).abs() < 1e-12, "{pole:?}");
    }

    /// Accurate to a few nanometres within 3 900 km of the central meridian (Karney 2011), on both
    /// sides of the equator, out to 90° and back round the far side.
    #[test]
    fn round_trips_far_from_the_central_meridian() {
        let tm = TransverseMercator::UTM;
        let mut worst = 0.0f64;
        for i in 0..=36 {
            for j in 0..=60 {
                let g = Geodetic::new(-89.5 + 4.97 * i as f64, -30.0 + j as f64, 0.0);
                let p = tm.forward(0.0, g);
                if p.x_m.abs() > 3.9e6 {
                    continue;
                }
                let back = tm.reverse(0.0, p.x_m, p.y_m);
                worst = worst.max(apart(g, back));
            }
        }
        assert!(worst < 5e-9, "{worst}");
        // The far side (more than 90° from the central meridian) is folded back consistently.
        let g = Geodetic::new(10.0, 120.0, 0.0);
        let p = tm.forward(0.0, g);
        let back = tm.reverse(0.0, p.x_m, p.y_m);
        assert!(back.lat_deg.is_finite() && back.lon_deg.is_finite(), "{back:?}");
    }

    #[test]
    fn reverse_reports_the_same_scale_and_convergence() {
        let tm = TransverseMercator::new(crate::GRS80, 1.0);
        for &(lat, lon) in &[(51.2, 3.7), (-12.0, -8.0), (0.0, 25.0), (75.0, 1.0)] {
            let p = tm.forward(0.0, Geodetic::new(lat, lon, 0.0));
            let (_, gamma, k) = tm.reverse_full(0.0, p.x_m, p.y_m);
            assert!(
                (gamma - p.convergence_deg).abs() < 1e-12 && (k - p.scale).abs() < 1e-14,
                "{lat} {lon}"
            );
        }
    }
}
