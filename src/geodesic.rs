//! The shortest path along the ellipsoid surface — Karney's algorithm, as in GeographicLib.
//!
//! Accurate to about 15 nanometres on WGS84, and it converges for every pair of points, including
//! ones on nearly opposite sides of the Earth (where Vincenty's older method fails).
//!
//! Ported from GeographicLib by Charles F. F. Karney, via its Rust port geographiclib-rs; both are
//! MIT-licensed (see the NOTICE file). Reference: C. F. F. Karney, *Algorithms for geodesics*,
//! J. Geodesy 87, 43–55 (2013), <https://doi.org/10.1007/s00190-012-0578-z>.

use core::f64::consts::PI;

use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::local::wrap_360;
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;
use crate::vector::fma;

/// The shortest path between two positions along the ellipsoid surface.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Geodesic {
    /// Length along the surface, metres.
    pub distance_m: f64,
    /// Bearing at the start, degrees clockwise from north in `[0, 360)`.
    pub initial_bearing_deg: f64,
    /// Bearing on arrival, degrees clockwise from north in `[0, 360)`.
    pub final_bearing_deg: f64,
}

/// The size of a polygon on the ellipsoid surface, from [`Ellipsoid::polygon_area`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PolygonArea {
    /// Area, square metres. Positive when the vertices go counter-clockwise (seen from above),
    /// negative when they go clockwise.
    pub area_m2: f64,
    /// Length of the boundary, metres.
    pub perimeter_m: f64,
}

impl Ellipsoid {
    /// The surface distance and bearings between two positions, accurate to about 15 nm. Heights
    /// are ignored (the path is on the surface). Works for every pair of points, including ones on
    /// opposite sides of the Earth.
    ///
    /// Returns `None` only if a latitude or longitude is NaN or infinite, or a latitude is beyond
    /// ±90°. Heights are ignored, so an unknown (NaN) height still gets a distance. Coincident
    /// points give a zero-length geodesic with bearings 0.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// // A quarter of the equator is exactly a·π/2 on the ellipsoid.
    /// let g = WGS84.geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 90.0, 0.0)).unwrap();
    /// assert!((g.distance_m - WGS84.a() * std::f64::consts::FRAC_PI_2).abs() < 1e-8);
    ///
    /// // Exactly opposite points on the equator: the shortest path runs over a pole.
    /// let half = WGS84.geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 180.0, 0.0)).unwrap();
    /// assert!((half.distance_m - 20_003_931.458_6).abs() < 1e-3);
    /// ```
    #[doc(alias = "vincenty")]
    #[doc(alias = "karney")]
    #[doc(alias = "vdist")]
    #[doc(alias = "distance")]
    pub fn geodesic_inverse(self, from: Geodetic, to: Geodetic) -> Option<Geodesic> {
        if from.lat_deg == to.lat_deg && from.lon_deg == to.lon_deg {
            return Some(Geodesic::default());
        }
        let (s12, salp1, calp1, salp2, calp2, _) = Karney::with(self, |k| {
            k.inverse(from.lat_deg, from.lon_deg, to.lat_deg, to.lon_deg, false)
        });
        s12.is_finite().then(|| Geodesic {
            distance_m: s12,
            initial_bearing_deg: wrap_360(atan2d(salp1, calp1)),
            final_bearing_deg: wrap_360(atan2d(salp2, calp2)),
        })
    }

    /// Travel `distance_m` along the surface from `from` on `bearing_deg`: returns where you arrive
    /// and the bearing on arrival, accurate to about 15 nm. Height is carried over. A negative
    /// distance travels backwards. Non-finite input gives NaN.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// let start = Geodetic::new(-37.95, 144.42, 0.0);
    /// let (end, _) = WGS84.geodesic_direct(start, 306.87, 54_972.0);
    /// let back = WGS84.geodesic_inverse(start, end).unwrap();
    /// assert!((back.distance_m - 54_972.0).abs() < 1e-8 && (back.initial_bearing_deg - 306.87).abs() < 1e-10);
    /// ```
    #[doc(alias = "vincenty")]
    #[doc(alias = "karney")]
    #[doc(alias = "vreckon")]
    #[doc(alias = "reckon")]
    #[doc(alias = "destination")]
    pub fn geodesic_direct(self, from: Geodetic, bearing_deg: f64, distance_m: f64) -> (Geodetic, f64) {
        let (lat2, lon2, azi2) = Karney::with(self, |k| k.direct(from.lat_deg, from.lon_deg, bearing_deg, distance_m));
        (Geodetic::new(lat2, lon2, from.height_m), wrap_360(azi2))
    }

    /// The area and perimeter of a polygon on the ellipsoid surface, matching GeographicLib to
    /// within 0.1 m².
    /// Each edge is the shortest path (a geodesic) between neighbouring vertices, and the last
    /// vertex joins back to the first. Heights are ignored.
    ///
    /// - The area is positive when the vertices go counter-clockwise (seen from above), and
    ///   negative when they go clockwise. Use `area_m2.abs()` if the order doesn't matter to you.
    /// - Polygons around a pole or across the antimeridian work as they are.
    /// - Repeating the first vertex at the end, as GeoJSON does, changes nothing.
    /// - The area is always between −½ and +½ of the Earth's surface. Going counter-clockwise
    ///   round a region bigger than half the Earth is the same as going clockwise round the rest,
    ///   so you get the rest, negative.
    ///
    /// Returns `None` if a latitude or longitude is NaN or infinite, or a latitude is beyond ±90°.
    /// Fewer than three vertices give zero area.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// // A 1° × 1° square on the equator, counter-clockwise.
    /// let square = [
    ///     Geodetic::new(0.0, 0.0, 0.0),
    ///     Geodetic::new(0.0, 1.0, 0.0),
    ///     Geodetic::new(1.0, 1.0, 0.0),
    ///     Geodetic::new(1.0, 0.0, 0.0),
    /// ];
    /// let p = WGS84.polygon_area(&square).unwrap();
    /// assert!((p.area_m2 - 12_308_778_361.469).abs() < 0.01);   // about 12 309 km²
    /// assert!((p.perimeter_m - 443_770.917).abs() < 0.001);
    /// ```
    #[doc(alias = "area")]
    #[doc(alias = "planimeter")]
    #[doc(alias = "areaint")]
    pub fn polygon_area(self, vertices: &[Geodetic]) -> Option<PolygonArea> {
        if vertices
            .iter()
            .any(|v| v.lat_deg.is_nan() || v.lat_deg.abs() > 90.0 || !v.lon_deg.is_finite())
        {
            return None;
        }
        Karney::with(self, |k| {
            let (mut area, mut perimeter, mut crossings) = (Accumulator::default(), Accumulator::default(), false);
            for (i, v1) in vertices.iter().enumerate() {
                let v2 = vertices[(i + 1) % vertices.len()];
                let (s12, _, _, _, _, s12_area) = k.inverse(v1.lat_deg, v1.lon_deg, v2.lat_deg, v2.lon_deg, true);
                perimeter.add(s12);
                area.add(s12_area);
                crossings ^= crosses_prime_meridian(v1.lon_deg, v2.lon_deg);
            }
            // Each edge adds the area between it and the equator, so the sum is the polygon's area
            // up to whole Earths, plus half an Earth for each pole it goes round (Karney 2013, §6).
            let earth = 4.0 * PI * k.c2();
            area.remainder(earth);
            if crossings {
                area.add(if area.s < 0.0 { earth / 2.0 } else { -earth / 2.0 });
            }
            // The sum is clockwise-positive; flip it, then keep it within ±half an Earth.
            let mut area = Accumulator { s: -area.s, t: -area.t };
            if area.s > earth / 2.0 {
                area.add(-earth);
            } else if area.s <= -earth / 2.0 {
                area.add(earth);
            }
            Some(PolygonArea {
                area_m2: 0.0 + area.s,
                perimeter_m: perimeter.s,
            })
        })
    }
}

/// A running sum kept to about twice `f64` precision (GeographicLib's `Accumulator`), so a
/// polygon's many edge terms, which can be far larger than its area, add up without rounding loss.
#[derive(Clone, Copy, Default)]
struct Accumulator {
    s: f64,
    t: f64,
}

impl Accumulator {
    fn add(&mut self, y: f64) {
        let (y, u) = two_sum(y, self.t);
        (self.s, self.t) = two_sum(y, self.s);
        if self.s == 0.0 {
            self.s = u;
        } else {
            self.t += u;
        }
    }

    /// Reduce the sum to about `[-y/2, y/2]` by whole multiples of `y`, exactly.
    fn remainder(&mut self, y: f64) {
        let r = self.s % y; // exact
        self.s = if r.abs() > y / 2.0 { r - y.copysign(r) } else { r }; // exact, by Sterbenz's lemma
        self.add(0.0);
    }
}

/// Whether the edge from `lon1` to `lon2` crosses the prime meridian, in either direction
/// (GeographicLib's `transit`; only the parity of the count matters). Longitude ±0 counts as east.
fn crosses_prime_meridian(lon1: f64, lon2: f64) -> bool {
    let (lon12, _) = ang_diff(lon1, lon2);
    let (lon1, lon2) = (ang_normalize(lon1), ang_normalize(lon2));
    (lon12 > 0.0 && ((lon1 < 0.0 && lon2 >= 0.0) || (lon1 > 0.0 && lon2 == 0.0)))
        || (lon12 < 0.0 && lon1 >= 0.0 && lon2 < 0.0)
}

// ── Karney's algorithm ───────────────────────────────────────────────────────────────────────────
//
// Names follow GeographicLib so the port can be checked against it line by line: bet = reduced
// latitude, alp = azimuth, sig = arc length on the auxiliary sphere, omg = longitude on it, lam =
// longitude, eps = the expansion parameter; s/c prefixes are sines and cosines.

/// Order of the series expansions (GeographicLib's default; full double precision).
const ORDER: usize = 6;
/// Number of coefficients in the C3 series.
const N_C3: usize = 15;
/// Number of coefficients in the C4 series (for areas).
const N_C4: usize = 21;

/// √(f64::MIN_POSITIVE) = 2⁻⁵¹¹: below this, quantities are treated as zero.
const TINY: f64 = 1.491_668_146_240_041_3e-154;
const TOL0: f64 = f64::EPSILON;
const TOL1: f64 = 200.0 * TOL0;
/// √ε = 2⁻²⁶.
const TOL2: f64 = 1.490_116_119_384_765_6e-8;
const TOLB: f64 = TOL0 * TOL2;
const XTHRESH: f64 = 1000.0 * TOL2;
/// Newton iterations before falling back to bisection, and the overall cap.
const MAXIT1: usize = 20;
const MAXIT2: usize = MAXIT1 + f64::MANTISSA_DIGITS as usize + 10;

const COEFF_A3: [f64; 18] = [
    -3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0, 2.0, 1.0, 1.0,
];

const COEFF_C3: [f64; 45] = [
    3.0, 128.0, 2.0, 5.0, 128.0, -1.0, 3.0, 3.0, 64.0, -1.0, 0.0, 1.0, 8.0, -1.0, 1.0, 4.0, 5.0, 256.0, 1.0, 3.0,
    128.0, -3.0, -2.0, 3.0, 64.0, 1.0, -3.0, 2.0, 32.0, 7.0, 512.0, -10.0, 9.0, 384.0, 5.0, -9.0, 5.0, 192.0, 7.0,
    512.0, -14.0, 7.0, 512.0, 21.0, 2560.0,
];

const COEFF_C4: [f64; 77] = [
    97.0, 15015.0, 1088.0, 156.0, 45045.0, -224.0, -4784.0, 1573.0, 45045.0, -10656.0, 14144.0, -4576.0, -858.0,
    45045.0, 64.0, 624.0, -4576.0, 6864.0, -3003.0, 15015.0, 100.0, 208.0, 572.0, 3432.0, -12012.0, 30030.0, 45045.0,
    1.0, 9009.0, -2944.0, 468.0, 135135.0, 5792.0, 1040.0, -1287.0, 135135.0, 5952.0, -11648.0, 9152.0, -2574.0,
    135135.0, -64.0, -624.0, 4576.0, -6864.0, 3003.0, 135135.0, 8.0, 10725.0, 1856.0, -936.0, 225225.0, -8448.0,
    4992.0, -1144.0, 225225.0, -1440.0, 4160.0, -4576.0, 1716.0, 225225.0, -136.0, 63063.0, 1024.0, -208.0, 105105.0,
    3584.0, -3328.0, 1144.0, 315315.0, -128.0, 135135.0, -2560.0, 832.0, 405405.0, 128.0, 99099.0,
];

/// The ellipsoid's derived constants and series coefficients.
#[derive(Debug, PartialEq)]
struct Karney {
    a: f64,
    f: f64,
    f1: f64,
    ep2: f64,
    n: f64,
    b: f64,
    e2: f64,
    a3x: [f64; ORDER],
    c3x: [f64; N_C3],
    c4x: [f64; N_C4],
}

/// The result of one evaluation of λ₁₂(α₁) in the inverse iteration.
struct Lambda {
    lam12: f64,
    salp2: f64,
    calp2: f64,
    sig12: f64,
    ssig1: f64,
    csig1: f64,
    ssig2: f64,
    csig2: f64,
    eps: f64,
    /// λ₁₂ − η: the ellipsoid's correction to the longitude difference.
    domg12: f64,
    dlam12: f64,
}

/// The A3, C3 and C4 series coefficients for third flattening `n`. A `const fn` (plain
/// multiply-add, so compile-time and run-time results are identical), so WGS84's are built into
/// the binary.
const fn series_coefficients(n: f64) -> ([f64; ORDER], [f64; N_C3], [f64; N_C4]) {
    const fn polyval_at(m: usize, p: &[f64], o: usize, x: f64) -> f64 {
        let (mut y, mut i) = (p[o], 1);
        while i <= m {
            y = y * x + p[o + i];
            i += 1;
        }
        y
    }
    const fn min(a: usize, b: usize) -> usize {
        if a < b {
            a
        } else {
            b
        }
    }
    let mut a3x = [0.0; ORDER];
    let (mut o, mut k) = (0, 0);
    while k < ORDER {
        let j = ORDER - 1 - k;
        let m = min(j, ORDER - j - 1);
        a3x[k] = polyval_at(m, &COEFF_A3, o, n) / COEFF_A3[o + m + 1];
        o += m + 2;
        k += 1;
    }
    let mut c3x = [0.0; N_C3];
    let (mut o, mut k, mut l) = (0, 0, 1);
    while l < ORDER {
        let mut j = ORDER - 1;
        while j >= l {
            let m = min(j, ORDER - j - 1);
            c3x[k] = polyval_at(m, &COEFF_C3, o, n) / COEFF_C3[o + m + 1];
            k += 1;
            o += m + 2;
            j -= 1;
        }
        l += 1;
    }
    let mut c4x = [0.0; N_C4];
    let (mut o, mut k, mut l) = (0, 0, 0);
    while l < ORDER {
        let mut j = ORDER - 1;
        while j >= l {
            let m = ORDER - j - 1;
            c4x[k] = polyval_at(m, &COEFF_C4, o, n) / COEFF_C4[o + m + 1];
            k += 1;
            o += m + 2;
            if j == 0 {
                break;
            }
            j -= 1;
        }
        l += 1;
    }
    (a3x, c3x, c4x)
}

static WGS84_KARNEY: Karney = Karney::new(WGS84.a(), WGS84.f());

impl Karney {
    /// The solver for `ell`, handed to `run`. WGS84's is built at compile time.
    fn with<R>(ell: Ellipsoid, run: impl FnOnce(&Karney) -> R) -> R {
        if ell == WGS84 {
            run(&WGS84_KARNEY)
        } else {
            run(&Karney::new(ell.a(), ell.f()))
        }
    }

    const fn new(a: f64, f: f64) -> Self {
        let f1 = 1.0 - f;
        let e2 = f * (2.0 - f);
        let n = f / (2.0 - f);
        let (a3x, c3x, c4x) = series_coefficients(n);
        Karney {
            a,
            f,
            f1,
            ep2: e2 / (f1 * f1),
            n,
            b: a * f1,
            e2,
            a3x,
            c3x,
            c4x,
        }
    }

    /// The threshold below which a short line is solved directly on the auxiliary sphere.
    fn etol2(&self) -> f64 {
        0.1 * TOL2 / (self.f.abs().max(0.001) * (1.0 - self.f / 2.0).min(1.0) / 2.0).sqrt()
    }

    /// The authalic radius squared: the ellipsoid's surface area is `4π·c2`.
    fn c2(&self) -> f64 {
        if self.e2 == 0.0 {
            return sq(self.a);
        }
        // e·atanh(e) / e², with e imaginary (so atanh becomes atan) on a prolate ellipsoid.
        let e = self.e2.abs().sqrt();
        let ratio = if self.e2 > 0.0 { e.atanh() / e } else { e.atan() / e };
        (sq(self.a) + sq(self.b) * ratio) / 2.0
    }

    fn a3f(&self, eps: f64) -> f64 {
        polyval(ORDER - 1, &self.a3x, eps)
    }

    fn c3f(&self, eps: f64, c: &mut [f64; ORDER]) {
        let (mut mult, mut o) = (1.0, 0);
        for (l, cl) in c.iter_mut().enumerate().skip(1) {
            let m = ORDER - l - 1;
            mult *= eps;
            *cl = mult * polyval(m, &self.c3x[o..], eps);
            o += m + 1;
        }
    }

    fn c4f(&self, eps: f64, c: &mut [f64; ORDER]) {
        let (mut mult, mut o) = (1.0, 0);
        for (l, cl) in c.iter_mut().enumerate() {
            let m = ORDER - l - 1;
            *cl = mult * polyval(m, &self.c4x[o..], eps);
            o += m + 1;
            mult *= eps;
        }
    }

    /// The scaled distance `s12/b` (when `distance`), and the scaled reduced length `m12/b` and
    /// `m0` (when `reduced`); NaN for anything not asked for.
    #[allow(clippy::too_many_arguments)]
    fn lengths(
        &self,
        eps: f64,
        sig12: f64,
        (ssig1, csig1, dn1): (f64, f64, f64),
        (ssig2, csig2, dn2): (f64, f64, f64),
        distance: bool,
        reduced: bool,
        c1a: &mut [f64; ORDER + 1],
        c2a: &mut [f64; ORDER + 1],
    ) -> (f64, f64, f64) {
        let (mut s12b, mut m12b, mut m0) = (f64::NAN, f64::NAN, f64::NAN);
        let (mut a1, mut a2, mut m0x, mut j12) = (0.0, 0.0, 0.0, 0.0);
        if distance || reduced {
            a1 = a1m1f(eps);
            c1f(eps, c1a);
            if reduced {
                a2 = a2m1f(eps);
                c2f(eps, c2a);
                m0x = a1 - a2;
                a2 += 1.0;
            }
            a1 += 1.0;
        }
        if distance {
            let b1 = sin_cos_series(true, ssig2, csig2, c1a) - sin_cos_series(true, ssig1, csig1, c1a);
            s12b = a1 * (sig12 + b1);
            if reduced {
                let b2 = sin_cos_series(true, ssig2, csig2, c2a) - sin_cos_series(true, ssig1, csig1, c2a);
                j12 = m0x * sig12 + (a1 * b1 - a2 * b2);
            }
        } else if reduced {
            for l in 1..=ORDER {
                c2a[l] = a1 * c1a[l] - a2 * c2a[l];
            }
            j12 = m0x * sig12 + (sin_cos_series(true, ssig2, csig2, c2a) - sin_cos_series(true, ssig1, csig1, c2a));
        }
        if reduced {
            m0 = m0x;
            m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * j12;
        }
        (s12b, m12b, m0)
    }

    /// The inverse problem: `(s12, salp1, calp1, salp2, calp2, S12)` for two positions in degrees.
    /// `S12` is the area between the geodesic and the equator in m², computed only when `area`
    /// (NaN otherwise).
    fn inverse(&self, lat1: f64, lon1: f64, lat2: f64, lon2: f64, area: bool) -> (f64, f64, f64, f64, f64, f64) {
        let Karney {
            a, f, f1, ep2, n, b, ..
        } = *self;
        let (mut c1a, mut c2a, mut c3a) = ([0.0; ORDER + 1], [0.0; ORDER + 1], [0.0; ORDER]);

        // Reduce to lon12 in [0, 180] and lat1 <= -|lat2| (Karney 2013, §4), remembering the signs.
        let (lon12, lon12s) = ang_diff(lon1, lon2);
        let mut lonsign = if lon12 >= 0.0 { 1.0 } else { -1.0 };
        let lon12 = lonsign * ang_round(lon12);
        let lon12s = ang_round((180.0 - lon12) - lonsign * lon12s);
        let lam12 = lon12.to_radians();
        let (slam12, clam12) = if lon12 > 90.0 {
            let (s, c) = sincosd(lon12s);
            (s, -c)
        } else {
            sincosd(lon12)
        };
        let mut lat1 = ang_round(lat_fix(lat1));
        let mut lat2 = ang_round(lat_fix(lat2));
        let swapp = if lat1.abs() < lat2.abs() { -1.0 } else { 1.0 };
        if swapp < 0.0 {
            lonsign = -lonsign;
            core::mem::swap(&mut lat1, &mut lat2);
        }
        let latsign = if lat1 < 0.0 { 1.0 } else { -1.0 };
        lat1 *= latsign;
        lat2 *= latsign;

        let (sbet1, cbet1) = {
            let (s, c) = sincosd(lat1);
            let (s, c) = norm2(f1 * s, c);
            (s, c.max(TINY))
        };
        let (mut sbet2, mut cbet2) = {
            let (s, c) = sincosd(lat2);
            let (s, c) = norm2(f1 * s, c);
            (s, c.max(TINY))
        };
        // Make |bet1| = |bet2| exact when it should be, so the symmetric cases come out exactly.
        if cbet1 < -sbet1 {
            if cbet2 == cbet1 {
                sbet2 = if sbet2 < 0.0 { sbet1 } else { -sbet1 };
            }
        } else if sbet2.abs() == -sbet1 {
            cbet2 = cbet1;
        }
        let dn1 = (1.0 + ep2 * sq(sbet1)).sqrt();
        let dn2 = (1.0 + ep2 * sq(sbet2)).sqrt();

        let (mut salp1, mut calp1, mut salp2, mut calp2) = (0.0, 0.0, 0.0, 0.0);
        let mut s12x = 0.0;
        let mut meridian = lat1 == -90.0 || slam12 == 0.0;
        // The longitude difference on the auxiliary sphere, for the area: set as an angle, or as
        // its sine and cosine (somg12 > 1 means "not set").
        let (mut omg12, mut somg12, mut comg12) = (0.0, 2.0, 0.0);

        if meridian {
            // Along a meridian the geodesic is the meridian itself, unless the line is long enough
            // to be beaten by a path over the pole (then fall through to the general case).
            calp1 = clam12;
            salp1 = slam12;
            calp2 = 1.0;
            salp2 = 0.0;
            let (ssig1, csig1) = (sbet1, calp1 * cbet1);
            let (ssig2, csig2) = (sbet2, calp2 * cbet2);
            let sig12 = (csig1 * ssig2 - ssig1 * csig2)
                .max(0.0)
                .atan2(csig1 * csig2 + ssig1 * ssig2);
            let (s, m12x, _) = self.lengths(
                n,
                sig12,
                (ssig1, csig1, dn1),
                (ssig2, csig2, dn2),
                true,
                true,
                &mut c1a,
                &mut c2a,
            );
            if sig12 < 1.0 || m12x >= 0.0 {
                s12x = if sig12 < 3.0 * TINY { 0.0 } else { s * b };
            } else {
                meridian = false;
            }
        }

        if !meridian && sbet1 == 0.0 && (f <= 0.0 || lon12s >= f * 180.0) {
            // Along the equator.
            calp1 = 0.0;
            calp2 = 0.0;
            salp1 = 1.0;
            salp2 = 1.0;
            s12x = a * lam12;
            omg12 = lam12 / f1;
        } else if !meridian {
            let (sig12, sa1, ca1, sa2, ca2, dnm) = self.inverse_start(
                sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, &mut c1a, &mut c2a,
            );
            (salp1, calp1, salp2, calp2) = (sa1, ca1, sa2, ca2);
            if sig12 >= 0.0 {
                // A short line, solved directly.
                s12x = sig12 * b * dnm;
                omg12 = lam12 / (f1 * dnm);
            } else {
                // Newton's method on alp1, falling back to bisection on a bracket.
                let (mut tripn, mut tripb) = (false, false);
                let (mut salp1a, mut calp1a, mut salp1b, mut calp1b) = (TINY, 1.0, TINY, -1.0);
                let mut last = None;
                for numit in 0..MAXIT2 {
                    let l = self.lambda12(
                        (sbet1, cbet1, dn1),
                        (sbet2, cbet2, dn2),
                        salp1,
                        calp1,
                        slam12,
                        clam12,
                        numit < MAXIT1,
                        &mut c1a,
                        &mut c2a,
                        &mut c3a,
                    );
                    let (v, dv) = (l.lam12, l.dlam12);
                    salp2 = l.salp2;
                    calp2 = l.calp2;
                    last = Some(l);
                    let tol = if tripn { 8.0 } else { 1.0 } * TOL0;
                    if tripb || v.is_nan() || v.abs() < tol {
                        break;
                    }
                    // Tighten the bracket [alp1a, alp1b] around the root.
                    if v > 0.0 && (numit > MAXIT1 || calp1 / salp1 > calp1b / salp1b) {
                        salp1b = salp1;
                        calp1b = calp1;
                    } else if v < 0.0 && (numit > MAXIT1 || calp1 / salp1 < calp1a / salp1a) {
                        salp1a = salp1;
                        calp1a = calp1;
                    }
                    if numit < MAXIT1 && dv > 0.0 {
                        let dalp1 = -v / dv;
                        let (sdalp1, cdalp1) = dalp1.sin_cos();
                        let nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
                        if nsalp1 > 0.0 && dalp1.abs() < PI {
                            (salp1, calp1) = norm2(nsalp1, calp1 * cdalp1 - salp1 * sdalp1);
                            tripn = v.abs() <= 16.0 * TOL0;
                            continue;
                        }
                    }
                    // Newton overshot or stalled: bisect the bracket.
                    (salp1, calp1) = norm2((salp1a + salp1b) / 2.0, (calp1a + calp1b) / 2.0);
                    tripn = false;
                    tripb = (salp1a - salp1).abs() + (calp1a - calp1) < TOLB
                        || (salp1 - salp1b).abs() + (calp1 - calp1b) < TOLB;
                }
                let l = last.expect("the loop runs at least once");
                let (s, _, _) = self.lengths(
                    l.eps,
                    l.sig12,
                    (l.ssig1, l.csig1, dn1),
                    (l.ssig2, l.csig2, dn2),
                    true,
                    false,
                    &mut c1a,
                    &mut c2a,
                );
                s12x = s * b;
                if area {
                    let (sdomg12, cdomg12) = l.domg12.sin_cos();
                    somg12 = slam12 * cdomg12 - clam12 * sdomg12;
                    comg12 = clam12 * cdomg12 + slam12 * sdomg12;
                }
            }
        }

        // The area between the geodesic and the equator (Karney 2013, §6).
        let mut s12_area = f64::NAN;
        if area {
            let salp0 = salp1 * cbet1;
            let calp0 = calp1.hypot(salp1 * sbet1);
            s12_area = if calp0 != 0.0 && salp0 != 0.0 {
                let (ssig1, csig1) = norm2(sbet1, calp1 * cbet1);
                let (ssig2, csig2) = norm2(sbet2, calp2 * cbet2);
                let k2 = sq(calp0) * ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                let mut c4a = [0.0; ORDER];
                self.c4f(eps, &mut c4a);
                let b4 = sin_cos_series(false, ssig2, csig2, &c4a) - sin_cos_series(false, ssig1, csig1, &c4a);
                sq(a) * calp0 * salp0 * self.e2 * b4
            } else {
                0.0
            };
            if !meridian && somg12 > 1.0 {
                (somg12, comg12) = omg12.sin_cos();
            }
            // The azimuth change along the line. For short lines this form avoids the cancellation
            // in alp2 − alp1.
            let alp12 = if !meridian && comg12 > -core::f64::consts::FRAC_1_SQRT_2 && sbet2 - sbet1 < 1.75 {
                let (domg12, dbet1, dbet2) = (1.0 + comg12, 1.0 + cbet1, 1.0 + cbet2);
                2.0 * (somg12 * (sbet1 * dbet2 + sbet2 * dbet1)).atan2(domg12 * (sbet1 * sbet2 + dbet1 * dbet2))
            } else {
                let (mut salp12, mut calp12) = (salp2 * calp1 - calp2 * salp1, calp2 * calp1 + salp2 * salp1);
                if salp12 == 0.0 && calp12 < 0.0 {
                    // Exactly antipodal azimuths: pick the sign the way GeographicLib does.
                    salp12 = TINY * calp1;
                    calp12 = -1.0;
                }
                salp12.atan2(calp12)
            };
            s12_area = (s12_area + self.c2() * alp12) * (swapp * lonsign * latsign) + 0.0;
        }

        // Undo the reductions.
        if swapp < 0.0 {
            core::mem::swap(&mut salp1, &mut salp2);
            core::mem::swap(&mut calp1, &mut calp2);
        }
        salp1 *= swapp * lonsign;
        calp1 *= swapp * latsign;
        salp2 *= swapp * lonsign;
        calp2 *= swapp * latsign;
        (0.0 + s12x, salp1, calp1, salp2, calp2, s12_area)
    }

    /// A starting azimuth for the inverse iteration, or the whole answer for short lines
    /// (`sig12 >= 0`): `(sig12, salp1, calp1, salp2, calp2, dnm)`.
    #[allow(clippy::too_many_arguments)]
    fn inverse_start(
        &self,
        sbet1: f64,
        cbet1: f64,
        dn1: f64,
        sbet2: f64,
        cbet2: f64,
        dn2: f64,
        lam12: f64,
        slam12: f64,
        clam12: f64,
        c1a: &mut [f64; ORDER + 1],
        c2a: &mut [f64; ORDER + 1],
    ) -> (f64, f64, f64, f64, f64, f64) {
        let Karney { f, f1, ep2, n, .. } = *self;
        let mut sig12 = -1.0;
        let (mut salp2, mut calp2, mut dnm) = (f64::NAN, f64::NAN, f64::NAN);
        let sbet12 = sbet2 * cbet1 - cbet2 * sbet1;
        let cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
        let sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
        let shortline = cbet12 >= 0.0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
        let (mut somg12, mut comg12) = if shortline {
            let mut sbetm2 = sq(sbet1 + sbet2);
            sbetm2 /= sbetm2 + sq(cbet1 + cbet2);
            dnm = (1.0 + ep2 * sbetm2).sqrt();
            (lam12 / (f1 * dnm)).sin_cos()
        } else {
            (slam12, clam12)
        };
        let mut salp1 = cbet2 * somg12;
        let mut calp1 = if comg12 >= 0.0 {
            sbet12 + cbet2 * sbet1 * sq(somg12) / (1.0 + comg12)
        } else {
            sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12)
        };
        let ssig12 = salp1.hypot(calp1);
        let csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

        if shortline && ssig12 < self.etol2() {
            // Really short: the great-circle solution on the auxiliary sphere is the answer.
            let c = if comg12 >= 0.0 {
                sq(somg12) / (1.0 + comg12)
            } else {
                1.0 - comg12
            };
            (salp2, calp2) = norm2(cbet1 * somg12, sbet12 - cbet1 * sbet2 * c);
            sig12 = ssig12.atan2(csig12);
        } else if n.abs() > 0.1 || csig12 >= 0.0 || ssig12 >= 6.0 * n.abs() * PI * sq(cbet1) {
            // Not nearly antipodal: the spherical estimate above is a good start.
        } else {
            // Nearly antipodal: start from the solution of the astroid problem (Karney 2013, §5).
            let lam12x = (-slam12).atan2(-clam12);
            let (x, y, lamscale);
            if f >= 0.0 {
                let k2 = sq(sbet1) * ep2;
                let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
                lamscale = f * cbet1 * self.a3f(eps) * PI;
                let betscale = lamscale * cbet1;
                x = lam12x / lamscale;
                y = sbet12a / betscale;
            } else {
                let cbet12a = cbet2 * cbet1 - sbet2 * sbet1;
                let bet12a = sbet12a.atan2(cbet12a);
                let (_, m12b, m0) = self.lengths(
                    n,
                    PI + bet12a,
                    (sbet1, -cbet1, dn1),
                    (sbet2, cbet2, dn2),
                    false,
                    true,
                    c1a,
                    c2a,
                );
                x = -1.0 + m12b / (cbet1 * cbet2 * m0 * PI);
                let betscale = if x < -0.01 { sbet12a / x } else { -f * sq(cbet1) * PI };
                lamscale = betscale / cbet1;
                y = lam12x / lamscale;
            }
            if y > -TOL1 && x > -1.0 - XTHRESH {
                if f >= 0.0 {
                    salp1 = (-x).min(1.0);
                    calp1 = -(1.0 - sq(salp1)).sqrt();
                } else {
                    calp1 = x.max(if x > -TOL1 { 0.0 } else { -1.0 });
                    salp1 = (1.0 - sq(calp1)).sqrt();
                }
            } else {
                let k = astroid(x, y);
                let omg12a = lamscale
                    * if f >= 0.0 {
                        -x * k / (1.0 + k)
                    } else {
                        -y * (1.0 + k) / k
                    };
                (somg12, comg12) = omg12a.sin_cos();
                comg12 = -comg12;
                salp1 = cbet2 * somg12;
                calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1.0 - comg12);
            }
        }
        if salp1 > 0.0 || salp1.is_nan() {
            (salp1, calp1) = norm2(salp1, calp1);
        } else {
            salp1 = 1.0;
            calp1 = 0.0;
        }
        (sig12, salp1, calp1, salp2, calp2, dnm)
    }

    /// λ₁₂ as a function of α₁ (and its derivative when `diffp`), for the inverse iteration.
    #[allow(clippy::too_many_arguments)]
    fn lambda12(
        &self,
        (sbet1, cbet1, dn1): (f64, f64, f64),
        (sbet2, cbet2, dn2): (f64, f64, f64),
        salp1: f64,
        calp1: f64,
        slam120: f64,
        clam120: f64,
        diffp: bool,
        c1a: &mut [f64; ORDER + 1],
        c2a: &mut [f64; ORDER + 1],
        c3a: &mut [f64; ORDER],
    ) -> Lambda {
        let Karney { f, f1, ep2, .. } = *self;
        // Break the degeneracy of an equatorial line heading due north.
        let calp1 = if sbet1 == 0.0 && calp1 == 0.0 { -TINY } else { calp1 };
        let salp0 = salp1 * cbet1;
        let calp0 = calp1.hypot(salp1 * sbet1);
        let somg1 = salp0 * sbet1;
        let comg1 = calp1 * cbet1;
        let (ssig1, csig1) = norm2(sbet1, comg1);

        let salp2 = if cbet2 != cbet1 { salp0 / cbet2 } else { salp1 };
        let calp2 = if cbet2 != cbet1 || sbet2.abs() != -sbet1 {
            let d = if cbet1 < -sbet1 {
                (cbet2 - cbet1) * (cbet1 + cbet2)
            } else {
                (sbet1 - sbet2) * (sbet1 + sbet2)
            };
            (sq(calp1 * cbet1) + d).sqrt() / cbet2
        } else {
            calp1.abs()
        };
        let somg2 = salp0 * sbet2;
        let comg2 = calp2 * cbet2;
        let (ssig2, csig2) = norm2(sbet2, comg2);

        let sig12 = (csig1 * ssig2 - ssig1 * csig2)
            .max(0.0)
            .atan2(csig1 * csig2 + ssig1 * ssig2);
        let somg12 = (comg1 * somg2 - somg1 * comg2).max(0.0);
        let comg12 = comg1 * comg2 + somg1 * somg2;
        let eta = (somg12 * clam120 - comg12 * slam120).atan2(comg12 * clam120 + somg12 * slam120);

        let k2 = sq(calp0) * ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);
        self.c3f(eps, c3a);
        let b312 = sin_cos_series(true, ssig2, csig2, c3a) - sin_cos_series(true, ssig1, csig1, c3a);
        let domg12 = -f * self.a3f(eps) * salp0 * (sig12 + b312);
        let lam12 = eta + domg12;

        let dlam12 = if !diffp {
            f64::NAN
        } else if calp2 == 0.0 {
            -2.0 * f1 * dn1 / sbet1
        } else {
            let (_, m12b, _) = self.lengths(
                eps,
                sig12,
                (ssig1, csig1, dn1),
                (ssig2, csig2, dn2),
                false,
                true,
                c1a,
                c2a,
            );
            m12b * f1 / (calp2 * cbet2)
        };
        Lambda {
            lam12,
            salp2,
            calp2,
            sig12,
            ssig1,
            csig1,
            ssig2,
            csig2,
            eps,
            domg12,
            dlam12,
        }
    }

    /// The direct problem: `(lat2, lon2, azi2)` in degrees after `s12` metres from `lat1, lon1` on
    /// azimuth `azi1`.
    fn direct(&self, lat1: f64, lon1: f64, azi1: f64, s12: f64) -> (f64, f64, f64) {
        let Karney { f, f1, ep2, b, .. } = *self;
        let (salp1, calp1) = sincosd(ang_round(ang_normalize(azi1)));
        let (sbet1, cbet1) = {
            let (s, c) = sincosd(ang_round(lat_fix(lat1)));
            let (s, c) = norm2(f1 * s, c);
            (s, c.max(TINY))
        };
        // The equatorial crossing (alp0) and the start on the auxiliary sphere (sig1, omg1).
        let salp0 = salp1 * cbet1;
        let calp0 = calp1.hypot(salp1 * sbet1);
        let somg1 = salp0 * sbet1;
        let comg1 = if sbet1 != 0.0 || calp1 != 0.0 {
            cbet1 * calp1
        } else {
            1.0
        };
        let (ssig1, csig1) = norm2(sbet1, comg1);
        let k2 = sq(calp0) * ep2;
        let eps = k2 / (2.0 * (1.0 + (1.0 + k2).sqrt()) + k2);

        let (mut c1a, mut c1pa, mut c3a) = ([0.0; ORDER + 1], [0.0; ORDER + 1], [0.0; ORDER]);
        let a1m1 = a1m1f(eps);
        c1f(eps, &mut c1a);
        c1pf(eps, &mut c1pa);
        self.c3f(eps, &mut c3a);
        let b11 = sin_cos_series(true, ssig1, csig1, &c1a);
        let (s, c) = b11.sin_cos();
        let (stau1, ctau1) = (ssig1 * c + csig1 * s, csig1 * c - ssig1 * s);
        let a3c = -f * salp0 * self.a3f(eps);
        let b31 = sin_cos_series(true, ssig1, csig1, &c3a);

        // Distance → arc length on the auxiliary sphere.
        let tau12 = s12 / (b * (1.0 + a1m1));
        let (s, c) = tau12.sin_cos();
        let b12 = -sin_cos_series(true, stau1 * c + ctau1 * s, ctau1 * c - stau1 * s, &c1pa);
        let mut sig12 = tau12 - (b12 - b11);
        let (mut ssig12, mut csig12) = sig12.sin_cos();
        if f.abs() > 0.01 {
            // One Newton step for strongly flattened ellipsoids, where the reverted series is short.
            let ssig2 = ssig1 * csig12 + csig1 * ssig12;
            let csig2 = csig1 * csig12 - ssig1 * ssig12;
            let b12 = sin_cos_series(true, ssig2, csig2, &c1a);
            let serr = (1.0 + a1m1) * (sig12 + (b12 - b11)) - s12 / b;
            sig12 -= serr / (1.0 + k2 * sq(ssig2)).sqrt();
            (ssig12, csig12) = sig12.sin_cos();
        }

        let ssig2 = ssig1 * csig12 + csig1 * ssig12;
        let mut csig2 = csig1 * csig12 - ssig1 * ssig12;
        let sbet2 = calp0 * ssig2;
        let mut cbet2 = salp0.hypot(calp0 * csig2);
        if cbet2 == 0.0 {
            // Arriving exactly at a pole.
            cbet2 = TINY;
            csig2 = TINY;
        }
        let (salp2, calp2) = (salp0, calp0 * csig2);
        let (somg2, comg2) = (salp0 * ssig2, csig2);
        let omg12 = (somg2 * comg1 - comg2 * somg1).atan2(comg2 * comg1 + somg2 * somg1);
        let lam12 = omg12 + a3c * (sig12 + (sin_cos_series(true, ssig2, csig2, &c3a) - b31));
        let lon2 = ang_normalize(ang_normalize(lon1) + ang_normalize(lam12.to_degrees()));
        (atan2d(sbet2, f1 * cbet2), lon2, atan2d(salp2, calp2))
    }
}

// ── Series and angle helpers ─────────────────────────────────────────────────────────────────────

#[inline]
fn sq(x: f64) -> f64 {
    x * x
}

/// `(x, y)` scaled to unit length.
#[inline]
fn norm2(x: f64, y: f64) -> (f64, f64) {
    let r = x.hypot(y);
    (x / r, y / r)
}

/// Error-free sum: `s + t == u + v` exactly.
#[inline]
fn two_sum(u: f64, v: f64) -> (f64, f64) {
    let s = u + v;
    let up = s - v;
    let vpp = s - up;
    (s, -((up - u) + (vpp - v)))
}

/// The polynomial `p[0]·xⁿ + … + p[n]`, by Horner's rule.
#[inline]
fn polyval(n: usize, p: &[f64], x: f64) -> f64 {
    p[1..=n].iter().fold(p[0], |y, &c| fma(y, x, c))
}

/// Tiny angles round to zero; the rest keep a resolution of 2⁻⁵⁷ (0.7 pm on the Earth), so
/// near-singular inputs like 1e-200° need no special handling.
fn ang_round(x: f64) -> f64 {
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
fn ang_normalize(x: f64) -> f64 {
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
fn lat_fix(x: f64) -> f64 {
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
fn sincosd(x: f64) -> (f64, f64) {
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
    let (s, c) = r.to_radians().sin_cos();
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

/// Clenshaw summation of `Σ c[l]·sin(2l·x)` (`sinp`) or `Σ c[l]·cos((2l+1)·x)`.
fn sin_cos_series<const N: usize>(sinp: bool, sinx: f64, cosx: f64, c: &[f64; N]) -> f64 {
    let mut k = N;
    let mut n = k - usize::from(sinp);
    let ar = 2.0 * (cosx - sinx) * (cosx + sinx);
    let mut y1 = 0.0;
    let mut y0 = if n % 2 == 1 {
        k -= 1;
        c[k]
    } else {
        0.0
    };
    n /= 2;
    while n > 0 {
        n -= 1;
        k -= 1;
        y1 = ar * y0 - y1 + c[k];
        k -= 1;
        y0 = ar * y1 - y0 + c[k];
    }
    if sinp {
        2.0 * sinx * cosx * y0
    } else {
        cosx * (y0 - y1)
    }
}

/// The positive root `k` of `k⁴ + 2k³ − (x² + y² − 1)k² − 2y²k − y² = 0` (Karney 2013, eq. 55).
fn astroid(x: f64, y: f64) -> f64 {
    let p = sq(x);
    let q = sq(y);
    let r = (p + q - 1.0) / 6.0;
    if q == 0.0 && r <= 0.0 {
        return 0.0;
    }
    let s = p * q / 4.0;
    let r2 = sq(r);
    let r3 = r * r2;
    let disc = s * (s + 2.0 * r3);
    let mut u = r;
    if disc >= 0.0 {
        let mut t3 = s + r3;
        // Pick the sign of the square root to avoid cancellation.
        t3 += if t3 < 0.0 { -disc.sqrt() } else { disc.sqrt() };
        let t = t3.cbrt();
        u += t + if t != 0.0 { r2 / t } else { 0.0 };
    } else {
        let ang = (-disc).sqrt().atan2(-(s + r3));
        u += 2.0 * r * (ang / 3.0).cos();
    }
    let v = (sq(u) + q).sqrt();
    let uv = if u < 0.0 { q / (v - u) } else { u + v };
    let w = (uv - q) / (2.0 * v);
    uv / ((uv + sq(w)).sqrt() + w)
}

/// `A₁ − 1` (Karney 2013, eq. 17).
fn a1m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [1.0, 4.0, 64.0, 0.0, 256.0];
    let m = ORDER / 2;
    let t = polyval(m, &COEFF, sq(eps)) / COEFF[m + 1];
    (t + eps) / (1.0 - eps)
}

/// `A₂ − 1` (Karney 2013, eq. 42).
fn a2m1f(eps: f64) -> f64 {
    const COEFF: [f64; 5] = [-11.0, -28.0, -192.0, 0.0, 256.0];
    let m = ORDER / 2;
    let t = polyval(m, &COEFF, sq(eps)) / COEFF[m + 1];
    (t - eps) / (1.0 + eps)
}

/// Fill `c[1..=ORDER]` with an eps series whose coefficients are packed in `coeff`.
fn eps_series(coeff: &[f64; 18], eps: f64, c: &mut [f64; ORDER + 1]) {
    let eps2 = sq(eps);
    let (mut d, mut o) = (eps, 0);
    for (l, cl) in c.iter_mut().enumerate().skip(1) {
        let m = (ORDER - l) / 2;
        *cl = d * polyval(m, &coeff[o..], eps2) / coeff[o + m + 1];
        o += m + 2;
        d *= eps;
    }
}

/// The `C₁` coefficients (Karney 2013, eq. 18).
fn c1f(eps: f64, c: &mut [f64; ORDER + 1]) {
    const COEFF: [f64; 18] = [
        -1.0, 6.0, -16.0, 32.0, -9.0, 64.0, -128.0, 2048.0, 9.0, -16.0, 768.0, 3.0, -5.0, 512.0, -7.0, 1280.0, -7.0,
        2048.0,
    ];
    eps_series(&COEFF, eps, c);
}

/// The `C′₁` coefficients, for turning distance back into arc length (Karney 2013, eq. 21).
fn c1pf(eps: f64, c: &mut [f64; ORDER + 1]) {
    const COEFF: [f64; 18] = [
        205.0, -432.0, 768.0, 1536.0, 4005.0, -4736.0, 3840.0, 12288.0, -225.0, 116.0, 384.0, -7173.0, 2695.0, 7680.0,
        3467.0, 7680.0, 38081.0, 61440.0,
    ];
    eps_series(&COEFF, eps, c);
}

/// The `C₂` coefficients (Karney 2013, eq. 43).
fn c2f(eps: f64, c: &mut [f64; ORDER + 1]) {
    const COEFF: [f64; 18] = [
        1.0, 2.0, 16.0, 32.0, 35.0, 64.0, 384.0, 2048.0, 15.0, 80.0, 768.0, 7.0, 35.0, 512.0, 63.0, 1280.0, 77.0,
        2048.0,
    ];
    eps_series(&COEFF, eps, c);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::{angle_close, close};

    /// The classic Flinders Peak → Buninyong line, with GeographicLib (Karney) as the reference on
    /// both WGS84 and Vincenty's own Australian National Spheroid: s = 54 972.271139 m /
    /// 54 972.469016 m, azimuths 306.868159203° → 307.173630629°.
    #[test]
    fn flinders_peak_to_buninyong_matches_geographiclib() {
        let flinders = Geodetic::new(
            -(37.0 + 57.0 / 60.0 + 3.72030 / 3600.0),
            144.0 + 25.0 / 60.0 + 29.52440 / 3600.0,
            0.0,
        );
        let buninyong = Geodetic::new(
            -(37.0 + 39.0 / 60.0 + 10.15610 / 3600.0),
            143.0 + 55.0 / 60.0 + 35.38390 / 3600.0,
            0.0,
        );
        for (ell, s12, azi1, azi2) in [
            (WGS84, 54_972.271_139, 306.868_159_203, 307.173_630_629),
            (
                Ellipsoid::new(6_378_160.0, 1.0 / 298.25),
                54_972.469_016,
                306.868_156_394,
                307.173_627_820,
            ),
        ] {
            let g = ell.geodesic_inverse(flinders, buninyong).unwrap();
            // The references are rounded to 1e-6 m and 1e-9°.
            assert!(close(g.distance_m, s12, 1e-6), "{}", g.distance_m);
            assert!(close(g.initial_bearing_deg, azi1, 1e-9), "{}", g.initial_bearing_deg);
            assert!(close(g.final_bearing_deg, azi2, 1e-9), "{}", g.final_bearing_deg);
            let (end, a2) = ell.geodesic_direct(flinders, g.initial_bearing_deg, g.distance_m);
            assert!(
                close(end.lat_deg, buninyong.lat_deg, 1e-12) && close(end.lon_deg, buninyong.lon_deg, 1e-12),
                "{end:?}"
            );
            assert!(close(a2, g.final_bearing_deg, 1e-11));
        }
    }

    /// Direct problems pinned to GeographicLib: 1 m and 1 km to 9 000 km, including a line that
    /// crosses the antimeridian.
    #[test]
    fn direct_matches_geographiclib() {
        let cases = [
            (0.0, 0.0, 0.0, 1.0, 9.043_694_770_503_82e-6, 0.0, 0.0),
            (
                51.5,
                -0.12,
                37.5,
                1_000.0,
                51.507_130_431_849_15,
                -0.111_231_971_285_650_47,
                37.506_862_270_433_38,
            ),
            (
                -33.9,
                151.2,
                200.0,
                500_000.0,
                -38.119_248_371_387_85,
                149.251_459_112_268_35,
                201.146_440_193_158_72,
            ),
            (
                10.0,
                179.0,
                271.0,
                9_000_000.0,
                2.513_907_928_159_937,
                97.852_725_808_884_27,
                260.298_633_547_235_97,
            ),
        ];
        for (lat, lon, brg, dist, lat2, lon2, azi2) in cases {
            let (to, a2) = WGS84.geodesic_direct(Geodetic::new(lat, lon, 0.0), brg, dist);
            assert!(
                close(to.lat_deg, lat2, 1e-12) && close(to.lon_deg, lon2, 1e-12),
                "{lat},{lon} {brg} {dist}: {to:?}"
            );
            assert!(angle_close(a2, azi2, 1e-11), "{a2} vs {azi2}");
        }
    }

    #[test]
    fn direct_and_inverse_round_trip_worldwide() {
        for (lat, lon) in [
            (0.0, 0.0),
            (51.5, -0.12),
            (-33.9, 151.2),
            (89.0, 45.0),
            (-80.0, -170.0),
            (10.0, 179.0),
        ] {
            let from = Geodetic::new(lat, lon, 12.0);
            for (bearing, dist) in [
                (0.0, 1.0),
                (37.5, 1_000.0),
                (200.0, 500_000.0),
                (271.0, 9_000_000.0),
                // Long, but short of the cut locus near the antipode: beyond it a geodesic is no
                // longer the shortest path, and the inverse rightly finds a shorter one.
                (123.0, 15_000_000.0),
            ] {
                let (to, a2) = WGS84.geodesic_direct(from, bearing, dist);
                assert_eq!(to.height_m, 12.0);
                let g = WGS84.geodesic_inverse(from, to).unwrap();
                assert!(
                    close(g.distance_m, dist, 1e-8),
                    "{from:?} {bearing} {dist}: {}",
                    g.distance_m
                );
                // Over 1 m the bearing is limited by position resolution (~1e-15° of latitude).
                let tol = if dist < 100.0 { 1e-6 } else { 1e-10 };
                assert!(
                    angle_close(g.initial_bearing_deg, bearing, tol),
                    "{from:?} {bearing} {dist}: {}",
                    g.initial_bearing_deg
                );
                assert!(angle_close(g.final_bearing_deg, a2, tol));
            }
        }
    }

    #[test]
    fn exact_lengths_on_the_equator_and_meridian() {
        let eq = WGS84
            .geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 90.0, 0.0))
            .unwrap();
        assert!(
            close(eq.distance_m, WGS84.a() * core::f64::consts::FRAC_PI_2, 1e-8),
            "{}",
            eq.distance_m
        );
        assert!(eq.initial_bearing_deg == 90.0 && eq.final_bearing_deg == 90.0);
        // Quarter meridian (equator to pole) on WGS84 = 10 001 965.729 312 7 m (GeographicLib).
        let mer = WGS84
            .geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(90.0, 0.0, 0.0))
            .unwrap();
        assert!(close(mer.distance_m, 10_001_965.729_312_7, 1e-6), "{}", mer.distance_m);
        assert_eq!(mer.initial_bearing_deg, 0.0);
    }

    /// Where Vincenty fails: points on (nearly) opposite sides of the Earth.
    #[test]
    fn antipodal_points_converge() {
        // Exactly antipodal on the equator: over a pole, twice the quarter meridian.
        let g = WGS84
            .geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 180.0, 0.0))
            .unwrap();
        assert!(
            close(g.distance_m, 2.0 * 10_001_965.729_312_7, 1e-6),
            "{}",
            g.distance_m
        );
        // Nearly antipodal pairs anywhere: the answer is found, and travelling it lands on the target.
        for (lat, lon, dlat, dlon) in [
            (0.0, 0.0, 0.5, 179.7),
            (30.0, 10.0, 0.1, 179.9),
            (-41.3, 174.8, 0.01, 179.99),
            (0.1, 0.0, 0.0, 179.5),
        ] {
            let a = Geodetic::new(lat, lon, 0.0);
            let b = Geodetic::new(-lat + dlat, lon + dlon, 0.0).normalized();
            let g = WGS84.geodesic_inverse(a, b).expect("Karney converges for every pair");
            assert!(g.distance_m > 19_000_000.0 && g.distance_m < 20_003_932.0, "{g:?}");
            let (end, _) = WGS84.geodesic_direct(a, g.initial_bearing_deg, g.distance_m);
            let miss = WGS84.to_ecef(end).distance_to(WGS84.to_ecef(b));
            assert!(miss < 1e-6, "{a:?} -> {b:?}: missed by {miss} m");
        }
    }

    #[test]
    fn degenerate_cases() {
        let p = Geodetic::new(10.0, 20.0, 0.0);
        assert_eq!(WGS84.geodesic_inverse(p, p), Some(Geodesic::default()));
        // Non-finite input or an impossible latitude has no geodesic.
        assert!(WGS84.geodesic_inverse(p, Geodetic::new(f64::NAN, 0.0, 0.0)).is_none());
        assert!(WGS84.geodesic_inverse(p, Geodetic::new(95.0, 0.0, 0.0)).is_none());
        assert!(WGS84.geodesic_direct(p, f64::NAN, 1.0).0.lat_deg.is_nan());
        // Pole to pole along any meridian.
        let g = WGS84
            .geodesic_inverse(Geodetic::new(90.0, 0.0, 0.0), Geodetic::new(-90.0, 0.0, 0.0))
            .unwrap();
        assert!(close(g.distance_m, 2.0 * 10_001_965.729_312_7, 1e-6));
        // Going backwards is going the other way.
        let (fwd, _) = WGS84.geodesic_direct(p, 30.0, -1_000.0);
        let (back, _) = WGS84.geodesic_direct(p, 210.0, 1_000.0);
        assert!(close(fwd.lat_deg, back.lat_deg, 1e-13) && close(fwd.lon_deg, back.lon_deg, 1e-13));
        // The spherical helper stays within 0.5 % of the geodesic on a long line.
        let (a, b) = (
            Geodetic::new(33.9425, -118.4081, 0.0),
            Geodetic::new(40.6413, -73.7781, 0.0),
        );
        let geo = WGS84.geodesic_inverse(a, b).unwrap().distance_m;
        let sph = WGS84.great_circle_distance_m(a, b);
        assert!((geo - sph).abs() / geo < 0.005, "geodesic {geo} vs spherical {sph}");
    }

    #[test]
    fn other_ellipsoids_including_a_sphere_and_prolate() {
        // On a sphere the geodesic is the great circle: R·central angle.
        let sphere = Ellipsoid::new(6_371_000.0, 0.0);
        let g = sphere
            .geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 60.0, 0.0))
            .unwrap();
        assert!(close(g.distance_m, 6_371_000.0 * core::f64::consts::FRAC_PI_3, 1e-8));
        // A prolate ellipsoid (negative flattening) round-trips too.
        let prolate = Ellipsoid::new(6_378_137.0, -1.0 / 150.0);
        let from = Geodetic::new(-20.0, 30.0, 0.0);
        for (brg, dist) in [(10.0, 1e5), (95.0, 1e7), (250.0, 1.95e7)] {
            let (to, _) = prolate.geodesic_direct(from, brg, dist);
            let back = prolate.geodesic_inverse(from, to).unwrap();
            assert!(close(back.distance_m, dist, 1e-7), "{brg} {dist}: {}", back.distance_m);
        }
    }

    #[test]
    fn the_compile_time_wgs84_solver_matches_one_built_at_run_time() {
        let (a, f) = std::hint::black_box((WGS84.a(), WGS84.f()));
        assert_eq!(WGS84_KARNEY, Karney::new(a, f));
    }

    #[test]
    fn helpers_are_exact_where_they_should_be() {
        assert_eq!(sincosd(90.0), (1.0, 0.0));
        assert_eq!(sincosd(-180.0), (-0.0, -1.0));
        assert_eq!(sincosd(-0.0), (-0.0, 1.0));
        assert_eq!(sincosd(810.0), (1.0, 0.0)); // through the large-angle path
        assert_eq!(sincosd(-45.0).0, -sincosd(45.0).0);
        assert!(sincosd(f64::INFINITY).0.is_nan() && sincosd(f64::NAN).1.is_nan());
        assert_eq!(atan2d(1.0, 0.0), 90.0);
        assert_eq!(atan2d(0.0, -1.0), 180.0);
        assert_eq!(atan2d(-1.0, -1.0), -135.0);
        assert_eq!(ang_normalize(540.0), 180.0);
        assert_eq!(ang_normalize(-180.0), 180.0);
        assert_eq!(ang_normalize(-190.0), 170.0);
        let (d, e) = ang_diff(179.0, -179.0);
        assert_eq!(d + e, 2.0);
    }

    fn polygon(points: &[(f64, f64)]) -> PolygonArea {
        let v: Vec<Geodetic> = points.iter().map(|&(lat, lon)| Geodetic::new(lat, lon, 0.0)).collect();
        WGS84.polygon_area(&v).unwrap()
    }

    /// GeographicLib's own planimeter tests (geographiclib-octave `geographiclib_test.m` and
    /// geographiclib-js `geodesictest.js`), with their published areas and perimeters.
    #[test]
    fn polygon_area_matches_geographiclib_planimeter_tests() {
        let check = |points: &[(f64, f64)], perimeter: f64, area: f64, tol: f64| {
            let p = polygon(points);
            assert!(
                (p.perimeter_m - perimeter).abs() < tol,
                "{points:?}: perimeter {}",
                p.perimeter_m
            );
            assert!((p.area_m2 - area).abs() < tol, "{points:?}: area {}", p.area_m2);
        };
        // Around the north pole, and the same around the south pole (clockwise seen from above).
        check(
            &[(89.0, 0.0), (89.0, 90.0), (89.0, 180.0), (89.0, 270.0)],
            631_819.874_5,
            24_952_305_678.0,
            1.0,
        );
        check(
            &[(-89.0, 0.0), (-89.0, 90.0), (-89.0, 180.0), (-89.0, 270.0)],
            631_819.874_5,
            -24_952_305_678.0,
            1.0,
        );
        check(
            &[(0.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 0.0)],
            627_598.273_1,
            24_619_419_146.0,
            1.0,
        );
        check(
            &[(90.0, 0.0), (0.0, 0.0), (0.0, 90.0)],
            30_022_685.0,
            63_758_202_715_511.0,
            1.0,
        );
        check(
            &[(89.0, 0.1), (89.0, 90.1), (89.0, -179.9)],
            539_297.0,
            12_476_152_838.5,
            1.0,
        );
        // Degenerate "polygons" that run along a parallel and back.
        check(&[(9.0, -1e-14), (9.0, 180.0), (9.0, 0.0)], 36_026_861.0, 0.0, 1.0);
        check(&[(9.0, 1e-14), (9.0, 0.0), (9.0, 180.0)], 36_026_861.0, 0.0, 1.0);
        check(&[(9.0, 1e-14), (9.0, 180.0), (9.0, 0.0)], 36_026_861.0, 0.0, 1.0);
        check(
            &[(66.562222222, 0.0), (66.562222222, 180.0), (66.562222222, 360.0)],
            10_465_729.0,
            0.0,
            1.0,
        );
        check(
            &[(66.562222222, -0.0), (66.562222222, -180.0), (66.562222222, -360.0)],
            10_465_729.0,
            0.0,
            1.0,
        );
        // Longitudes given outside [-180, 180].
        let ring = [
            (89.0, -360.0),
            (89.0, -240.0),
            (89.0, -120.0),
            (89.0, 0.0),
            (89.0, 120.0),
            (89.0, 240.0),
        ];
        check(&ring, 1_160_741.0, 32_415_230_256.0, 1.0);
        // Order decides the sign.
        assert!((polygon(&[(2.0, 1.0), (1.0, 2.0), (3.0, 3.0)]).area_m2 - 18_454_562_325.451_2).abs() < 1e-4);
        assert!((polygon(&[(1.0, 2.0), (2.0, 1.0), (3.0, 3.0)]).area_m2 + 18_454_562_325.451_2).abs() < 1e-4);
    }

    /// A triangle on the 45° parallel walked round three or four times: each circuit adds the same
    /// area (GeographicLib's planimeter21 test).
    #[test]
    fn polygon_area_counts_repeated_circuits() {
        let r = 39_433_884_866_571.43; // the area of one circuit
        for circuits in [3, 4] {
            let ring: Vec<(f64, f64)> = (0..circuits)
                .flat_map(|_| [(45.0, 60.0), (45.0, 180.0), (45.0, -60.0)])
                .collect();
            assert!((polygon(&ring).area_m2 - circuits as f64 * r).abs() < 0.5);
        }
    }

    #[test]
    fn polygon_area_edge_cases() {
        let square = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)];
        let p = polygon(&square);
        // Closing the ring explicitly (GeoJSON style) changes nothing.
        assert_eq!(
            polygon(&[(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]),
            p
        );
        // Reversing the order flips the sign.
        let mut reversed = square;
        reversed.reverse();
        assert_eq!(polygon(&reversed).area_m2, -p.area_m2);
        assert_eq!(polygon(&[]), PolygonArea::default());
        assert_eq!(polygon(&[(10.0, 20.0)]), PolygonArea::default());
        let line = polygon(&[(10.0, 20.0), (11.0, 21.0)]);
        assert_eq!(line.area_m2, 0.0);
        let one_way = WGS84
            .geodesic_inverse(Geodetic::new(10.0, 20.0, 0.0), Geodetic::new(11.0, 21.0, 0.0))
            .unwrap();
        assert!(close(line.perimeter_m, 2.0 * one_way.distance_m, 1e-9));
        // A whole hemisphere: the equator, counter-clockwise seen from the north.
        let earth = 4.0 * PI * WGS84_KARNEY.c2();
        let hemisphere = polygon(&[(0.0, 0.0), (0.0, 120.0), (0.0, -120.0)]);
        assert!((hemisphere.area_m2 - earth / 2.0).abs() < 1.0, "{hemisphere:?}");
        assert!(close(hemisphere.perimeter_m, 2.0 * PI * WGS84.a(), 1e-6));
        // Bad input.
        for bad in [f64::NAN, f64::INFINITY, 90.5] {
            assert!(WGS84
                .polygon_area(&[Geodetic::new(bad, 0.0, 0.0), Geodetic::new(1.0, 1.0, 0.0)])
                .is_none());
        }
        assert!(WGS84
            .polygon_area(&[Geodetic::new(1.0, f64::NAN, 0.0), Geodetic::new(1.0, 1.0, 0.0)])
            .is_none());
        // Heights are ignored.
        let raised: Vec<Geodetic> = square.iter().map(|&(lat, lon)| Geodetic::new(lat, lon, 1e4)).collect();
        assert_eq!(WGS84.polygon_area(&raised).unwrap(), p);
    }

    /// The Earth's area, 4π·c2: 510 065 621.7 km² for WGS84 (GeographicLib's `Geodesic::EllipsoidArea`).
    #[test]
    fn earth_area_matches_geographiclib() {
        assert!((4.0 * PI * WGS84_KARNEY.c2() - 510_065_621_724_088.5).abs() < 1.0);
        let sphere = Karney::new(6_371_000.0, 0.0);
        assert!(close(4.0 * PI * sphere.c2(), 4.0 * PI * 6_371_000.0f64.powi(2), 1e-6));
    }
}
