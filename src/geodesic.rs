//! The shortest path along the ellipsoid surface (Vincenty's formulae).

use crate::ellipsoid::Ellipsoid;
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

const VINCENTY_TOL: f64 = 1e-14;
const VINCENTY_MAX_ITER: usize = 200;

/// The σ-dependent state of Vincenty's inverse for one value of λ.
#[derive(Clone, Copy)]
struct InverseState {
    sin_sigma: f64,
    cos_sigma: f64,
    sigma: f64,
    sin_alpha: f64,
    cos2_alpha: f64,
    cos2_sigma_m: f64,
}

impl Ellipsoid {
    /// The surface distance and bearings between two positions — Vincenty's inverse formula,
    /// accurate to about 0.5 mm. Heights are ignored (the path is on the surface).
    ///
    /// Returns `None` for nearly antipodal points (on opposite sides of the Earth), where Vincenty's
    /// iteration does not converge; use
    /// [`great_circle_distance_m`](Self::great_circle_distance_m) there for an estimate.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// // A quarter of the equator is exactly a·π/2 on the ellipsoid.
    /// let g = WGS84.geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 90.0, 0.0)).unwrap();
    /// assert!((g.distance_m - WGS84.a() * std::f64::consts::FRAC_PI_2).abs() < 1e-3);
    /// ```
    #[doc(alias = "vincenty")]
    #[doc(alias = "vdist")]
    #[doc(alias = "distance")]
    pub fn geodesic_inverse(self, from: Geodetic, to: Geodetic) -> Option<Geodesic> {
        let f = self.f();
        let big_l = to.lon_rad() - from.lon_rad();
        let (su1, cu1) = self.reduced_latitude(from.lat_rad()).sin_cos();
        let (su2, cu2) = self.reduced_latitude(to.lat_rad()).sin_cos();
        let mut lambda = big_l;
        let mut state;
        let mut iter = 0;
        loop {
            state = match Self::inverse_state(lambda, su1, cu1, su2, cu2) {
                None => return Some(Geodesic::default()), // coincident points
                Some(st) => st,
            };
            let c = self.vincenty_c(state.cos2_alpha);
            let next = big_l
                + (1.0 - c)
                    * f
                    * state.sin_alpha
                    * (state.sigma
                        + c * state.sin_sigma
                            * (state.cos2_sigma_m
                                + c * state.cos_sigma * (-1.0 + 2.0 * state.cos2_sigma_m * state.cos2_sigma_m)));
            iter += 1;
            if !next.is_finite() || iter > VINCENTY_MAX_ITER {
                return None;
            }
            let done = (next - lambda).abs() < VINCENTY_TOL;
            lambda = next;
            if done {
                // Re-evaluate at the converged λ so nothing below uses a stale iterate.
                state = Self::inverse_state(lambda, su1, cu1, su2, cu2)?;
                break;
            }
        }
        let (big_a, big_b) = self.vincenty_ab(state.cos2_alpha);
        let delta_sigma = Self::vincenty_delta_sigma(big_b, state.sin_sigma, state.cos_sigma, state.cos2_sigma_m);
        let (sl, cl) = lambda.sin_cos();
        let a1 = (cu2 * sl).atan2(cu1 * su2 - su1 * cu2 * cl);
        let a2 = (cu1 * sl).atan2(-su1 * cu2 + cu1 * su2 * cl);
        Some(Geodesic {
            distance_m: self.b() * big_a * (state.sigma - delta_sigma),
            initial_bearing_deg: a1.to_degrees().rem_euclid(360.0),
            final_bearing_deg: a2.to_degrees().rem_euclid(360.0),
        })
    }

    /// Travel `distance_m` along the surface from `from` on `bearing_deg`: returns where you arrive
    /// and the bearing on arrival — Vincenty's direct formula (about 0.5 mm). Height is carried over.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// let start = Geodetic::new(-37.95, 144.42, 0.0);
    /// let (end, _) = WGS84.geodesic_direct(start, 306.87, 54_972.0);
    /// let back = WGS84.geodesic_inverse(start, end).unwrap();
    /// assert!((back.distance_m - 54_972.0).abs() < 1e-6 && (back.initial_bearing_deg - 306.87).abs() < 1e-9);
    /// ```
    #[doc(alias = "vincenty")]
    #[doc(alias = "vreckon")]
    #[doc(alias = "destination")]
    pub fn geodesic_direct(self, from: Geodetic, bearing_deg: f64, distance_m: f64) -> (Geodetic, f64) {
        let f = self.f();
        let b = self.b();
        let (sa1, ca1) = bearing_deg.to_radians().sin_cos();
        let u1 = self.reduced_latitude(from.lat_rad());
        let (su1, cu1) = u1.sin_cos();
        let sigma1 = u1.tan().atan2(ca1);
        let sin_alpha = cu1 * sa1;
        let cos2_alpha = 1.0 - sin_alpha * sin_alpha;
        let (big_a, big_b) = self.vincenty_ab(cos2_alpha);
        let base = distance_m / (b * big_a);
        let mut sigma = base;
        let mut iter = 0;
        loop {
            let (s, c) = sigma.sin_cos();
            let next = base + Self::vincenty_delta_sigma(big_b, s, c, (2.0 * sigma1 + sigma).cos());
            iter += 1;
            let done = (next - sigma).abs() < VINCENTY_TOL || iter >= VINCENTY_MAX_ITER;
            sigma = next;
            if done {
                break;
            }
        }
        // Everything below is evaluated at the converged σ.
        let (sin_sigma, cos_sigma) = sigma.sin_cos();
        let cos2_sigma_m = (2.0 * sigma1 + sigma).cos();
        let t = su1 * sin_sigma - cu1 * cos_sigma * ca1;
        let lat2 = (su1 * cos_sigma + cu1 * sin_sigma * ca1).atan2((1.0 - f) * (sin_alpha * sin_alpha + t * t).sqrt());
        let lambda = (sin_sigma * sa1).atan2(cu1 * cos_sigma - su1 * sin_sigma * ca1);
        let c = self.vincenty_c(cos2_alpha);
        let big_l = lambda
            - (1.0 - c)
                * f
                * sin_alpha
                * (sigma + c * sin_sigma * (cos2_sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m * cos2_sigma_m)));
        let lon2 = from.lon_rad() + big_l;
        let a2 = sin_alpha.atan2(-t);
        (
            Geodetic::new(lat2.to_degrees(), lon2.to_degrees(), from.height_m).normalized(),
            a2.to_degrees().rem_euclid(360.0),
        )
    }

    /// Reduced (parametric) latitude, radians.
    #[inline]
    fn reduced_latitude(self, lat_rad: f64) -> f64 {
        ((1.0 - self.f()) * lat_rad.tan()).atan()
    }

    #[inline]
    fn vincenty_ab(self, cos2_alpha: f64) -> (f64, f64) {
        let b = self.b();
        let u2 = cos2_alpha * (self.a() * self.a() - b * b) / (b * b);
        // Horner form of the series, fused where the target has FMA (one rounding per term).
        let big_a = fma(
            u2 / 16384.0,
            fma(u2, fma(u2, fma(-175.0, u2, 320.0), -768.0), 4096.0),
            1.0,
        );
        let big_b = (u2 / 1024.0) * fma(u2, fma(u2, fma(-47.0, u2, 74.0), -128.0), 256.0);
        (big_a, big_b)
    }

    #[inline]
    fn vincenty_c(self, cos2_alpha: f64) -> f64 {
        let f = self.f();
        f / 16.0 * cos2_alpha * fma(f, fma(-3.0, cos2_alpha, 4.0), 4.0)
    }

    #[inline]
    fn vincenty_delta_sigma(big_b: f64, sin_sigma: f64, cos_sigma: f64, cos2_sigma_m: f64) -> f64 {
        let c2m2 = cos2_sigma_m * cos2_sigma_m;
        let t_cos = fma(2.0, c2m2, -1.0);
        let t_sin = fma(4.0, sin_sigma * sin_sigma, -3.0);
        let t_c2m = fma(4.0, c2m2, -3.0);
        let inner = fma(cos_sigma, t_cos, -(big_b / 6.0 * cos2_sigma_m * t_sin * t_c2m));
        big_b * sin_sigma * fma(big_b / 4.0, inner, cos2_sigma_m)
    }

    /// One evaluation of the inverse state at `lambda`; `None` for coincident points.
    #[inline]
    fn inverse_state(lambda: f64, su1: f64, cu1: f64, su2: f64, cu2: f64) -> Option<InverseState> {
        let (sl, cl) = lambda.sin_cos();
        let t = cu1 * su2 - su1 * cu2 * cl;
        let sin_sigma = ((cu2 * sl) * (cu2 * sl) + t * t).sqrt();
        if sin_sigma == 0.0 {
            return None;
        }
        let cos_sigma = su1 * su2 + cu1 * cu2 * cl;
        let sigma = sin_sigma.atan2(cos_sigma);
        let sin_alpha = cu1 * cu2 * sl / sin_sigma;
        let cos2_alpha = 1.0 - sin_alpha * sin_alpha;
        let cos2_sigma_m = if cos2_alpha == 0.0 {
            0.0
        } else {
            cos_sigma - 2.0 * su1 * su2 / cos2_alpha
        };
        Some(InverseState {
            sin_sigma,
            cos_sigma,
            sigma,
            sin_alpha,
            cos2_alpha,
            cos2_sigma_m,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::{angle_close, close};
    use crate::WGS84;

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
            assert!(close(g.distance_m, s12, 5e-4), "{}", g.distance_m);
            assert!(close(g.initial_bearing_deg, azi1, 1e-8), "{}", g.initial_bearing_deg);
            assert!(close(g.final_bearing_deg, azi2, 1e-8), "{}", g.final_bearing_deg);
            let (end, a2) = ell.geodesic_direct(flinders, g.initial_bearing_deg, g.distance_m);
            assert!(
                close(end.lat_deg, buninyong.lat_deg, 1e-10) && close(end.lon_deg, buninyong.lon_deg, 1e-10),
                "{end:?}"
            );
            assert!(close(a2, g.final_bearing_deg, 1e-9));
        }
    }

    /// Direct problems pinned to GeographicLib: 1 m and 1 km to 9 000 km, including a line that
    /// crosses the antimeridian. Vincenty and Karney agree to ~0.5 mm ≈ 5e-9°.
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
                close(to.lat_deg, lat2, 1e-8) && close(to.lon_deg, lon2, 1e-8),
                "{lat},{lon} {brg} {dist}: {to:?}"
            );
            assert!(angle_close(a2, azi2, 1e-7), "{a2} vs {azi2}");
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
            for (bearing, dist) in [(0.0, 1.0), (37.5, 1_000.0), (200.0, 500_000.0), (271.0, 9_000_000.0)] {
                let (to, a2) = WGS84.geodesic_direct(from, bearing, dist);
                assert_eq!(to.height_m, 12.0);
                let g = WGS84.geodesic_inverse(from, to).unwrap();
                assert!(
                    close(g.distance_m, dist, 1e-6),
                    "{from:?} {bearing} {dist}: {}",
                    g.distance_m
                );
                // Bearing resolution is limited by position resolution over the line length:
                // ~3e-10 m lateral over 1 m is 2e-8°, so the tolerance scales with distance.
                let tol = if dist < 100.0 { 1e-6 } else { 1e-8 };
                assert!(
                    angle_close(g.initial_bearing_deg, bearing, tol),
                    "{}",
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
            close(eq.distance_m, WGS84.a() * core::f64::consts::FRAC_PI_2, 1e-4),
            "{}",
            eq.distance_m
        );
        assert!(close(eq.initial_bearing_deg, 90.0, 1e-9) && close(eq.final_bearing_deg, 90.0, 1e-9));
        // Quarter meridian (equator to pole) on WGS84 = 10 001 965.729 m.
        let mer = WGS84
            .geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(90.0, 0.0, 0.0))
            .unwrap();
        assert!(close(mer.distance_m, 10_001_965.729, 1e-3), "{}", mer.distance_m);
        assert!(close(mer.initial_bearing_deg, 0.0, 1e-9));
    }

    #[test]
    fn degenerate_cases() {
        let p = Geodetic::new(10.0, 20.0, 0.0);
        assert_eq!(WGS84.geodesic_inverse(p, p), Some(Geodesic::default()));
        // Exactly antipodal on the equator: either Vincenty gives up (None) or it returns the
        // over-the-pole path, twice the quarter meridian.
        match WGS84.geodesic_inverse(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 180.0, 0.0)) {
            None => {}
            Some(g) => assert!(close(g.distance_m, 2.0 * 10_001_965.729, 1.0), "{}", g.distance_m),
        }
        // The spherical helper stays within 0.5 % of the geodesic on a long line.
        let (a, b) = (
            Geodetic::new(33.9425, -118.4081, 0.0),
            Geodetic::new(40.6413, -73.7781, 0.0),
        );
        let geo = WGS84.geodesic_inverse(a, b).unwrap().distance_m;
        let sph = WGS84.great_circle_distance_m(a, b);
        assert!((geo - sph).abs() / geo < 0.005, "geodesic {geo} vs spherical {sph}");
    }
}
