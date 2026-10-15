//! Rays against the ellipsoid: where a line meets the surface, line of sight and horizon distance.
//! The surface is smooth, with no terrain. Refraction (the air bending the line of sight) is
//! optional.

use crate::ellipsoid::Ellipsoid;
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Ecef;
use crate::vector::{Dd, Vec3};

/// The refraction factor for radio waves (radar, radio links) in a standard atmosphere: the
/// effective Earth radius is 4/3 of the true one (ITU-R P.834). For visible light use about 7/6,
/// and 1 for no refraction.
pub const STANDARD_REFRACTION_K: f64 = 4.0 / 3.0;

impl Ellipsoid {
    /// Where the ray `origin + t·dir` (t ≥ 0) first meets the ellipsoid surface, or `None` if it
    /// misses or the surface is behind it. `dir` need not be unit length.
    ///
    /// For a pointing from an observer (azimuth and elevation) use
    /// [`LocalFrame::ground_strike`](crate::LocalFrame::ground_strike).
    pub fn ray_intersect(self, origin: Ecef, dir: Vec3) -> Option<Ecef> {
        let (t1, t2) = self.ray_roots(origin, dir)?;
        let t = if t1 >= 0.0 {
            t1
        } else if t2 >= 0.0 {
            t2
        } else {
            return None;
        };
        Some(origin.offset(dir * t))
    }

    /// The roots `t1 ≤ t2` where the line `origin + t·dir` crosses the ellipsoid surface, or `None`
    /// when it misses, `dir` is zero, or a coordinate is not finite. Scales space so the ellipsoid
    /// becomes the unit sphere and solves the ray–sphere quadratic in its cancellation-free form,
    /// so the near root stays accurate when the origin is close to the surface.
    pub(crate) fn ray_roots(self, origin: Ecef, dir: Vec3) -> Option<(f64, f64)> {
        let s = Vec3::new(1.0 / self.a(), 1.0 / self.a(), 1.0 / self.b());
        let o = Vec3::new(origin.x * s.x, origin.y * s.y, origin.z * s.z);
        let d = Vec3::new(dir.x * s.x, dir.y * s.y, dir.z * s.z);
        // The terms to 32 digits: o·o − 1 cancels near the surface, and the discriminant cancels
        // for rays that graze it.
        let (dd, od, oo1) = (Dd::dot(d, d), Dd::dot(o, d), Dd::dot(o, o) - Dd(1.0, 0.0));
        let disc = 4.0 * (od * od - dd * oo1).f64();
        let (qa, qb, qc) = (dd.f64(), 2.0 * od.f64(), oo1.f64());
        if disc.is_nan() || disc < 0.0 || qa == 0.0 {
            return None;
        }
        // q = −(b + sign(b)·√disc)/2 never subtracts nearly equal numbers; the roots are q/a and c/q.
        let q = -0.5 * (qb + disc.sqrt().copysign(qb));
        if q == 0.0 {
            return Some((0.0, 0.0)); // tangent at the origin, which is on the surface
        }
        let (r1, r2) = (q / qa, qc / q);
        Some((r1.min(r2), r1.max(r2)))
    }

    /// Whether the straight line between two points clears the ellipsoid: can `a` see `b`? Both
    /// points should be at or above the surface. A point on the surface sees another only if the
    /// line between them stays above the surface (to within a micrometre). For example, two
    /// points on the ground 1 km apart can't see each other: the Earth bulges 2 cm between them.
    pub fn line_of_sight(self, a: Ecef, b: Ecef) -> bool {
        let d = b.vec() - a.vec();
        if d.dot(d) == 0.0 {
            return true;
        }
        // The line is inside the ellipsoid between its two crossings t1 and t2. If that stretch
        // overlaps the segment [0, 1], check the height at its deepest point. The deepest point
        // is the middle of the stretch, clamped to the segment.
        let Some((t1, t2)) = self.ray_roots(a, d) else {
            return true; // the line misses the ellipsoid entirely
        };
        let (lo, hi) = (t1.max(0.0), t2.min(1.0));
        if lo >= hi {
            return true;
        }
        let deepest = ((t1 + t2) / 2.0).clamp(lo, hi);
        self.to_geodetic(a.offset(d * deepest)).height_m >= -1e-6
    }

    /// How far you can see from `height_m` above the surface before the Earth curves away, metres
    /// (on a sphere of the mean radius). Negative heights count as zero. The line of sight is
    /// straight; [`horizon_distance_with_refraction_m`](Self::horizon_distance_with_refraction_m)
    /// lets the air bend it.
    #[inline]
    pub fn horizon_distance_m(self, height_m: f64) -> f64 {
        let (r, h) = (self.mean_radius(), height_m.max(0.0));
        (2.0 * r * h + h * h).sqrt()
    }

    /// Like [`horizon_distance_m`](Self::horizon_distance_m), but with the air bending the line of
    /// sight around the Earth. `k` is the effective-Earth-radius factor:
    /// [`STANDARD_REFRACTION_K`] (4/3) for radio and radar, about 7/6 for visible light, and 1 for
    /// none.
    ///
    /// The bending happens in the lower atmosphere, as in
    /// [`line_of_sight_with_refraction`](Self::line_of_sight_with_refraction): from heights above
    /// 7.35 km the line of sight runs straight down to that layer, then bends through it.
    ///
    /// ```
    /// use geo3d::{STANDARD_REFRACTION_K, WGS84};
    /// let straight = WGS84.horizon_distance_m(100.0);                                    // 35.7 km
    /// let radio = WGS84.horizon_distance_with_refraction_m(100.0, STANDARD_REFRACTION_K); // 41.2 km
    /// assert!(radio > straight * 1.15);
    /// ```
    pub fn horizon_distance_with_refraction_m(self, height_m: f64, k: f64) -> f64 {
        let (r, h, c) = (self.mean_radius(), height_m.max(0.0), REFRACTION_CEILING_M);
        let kr = k * r;
        if h <= c {
            return (2.0 * kr * h + h * h).sqrt();
        }
        // Through the air, straight on an Earth k times larger, up to the top of the layer...
        let through_air = (2.0 * kr * c + c * c).sqrt();
        // ...where the path climbs at an angle θ, sin θ = through_air / (kr + c). From there it's
        // straight: the distance u from radius r + c to r + h, solved without cancelling.
        let rc = r + c;
        let rise = rc * through_air / (kr + c);
        let span = (h - c) * (2.0 * r + h + c);
        through_air + span / (rise + (rise * rise + span).sqrt())
    }

    /// Like [`line_of_sight`](Self::line_of_sight), but with the air bending the line of sight
    /// around the Earth, so radar and radio links see further. `k` is the effective-Earth-radius
    /// factor: [`STANDARD_REFRACTION_K`] for radio, about 7/6 for visible light, and 1 for none.
    ///
    /// It uses the standard effective-Earth model: the path curves around the Earth with curvature
    /// `(1 − 1/k)/R`, so between two points in the air its clearance gains `d₁·d₂·(1 − 1/k) / 2R` at
    /// distances d₁ and d₂ from the ends. The path is clear if its lowest point stays above the
    /// surface (to within a micrometre).
    ///
    /// The air bends the path only where there is enough of it: below 7.35 km, the refractivity
    /// scale height of ITU-R P.453's reference atmosphere. A path that rises above it, to a
    /// satellite say, bends only along its lower part. A distant satellite is then seen down to
    /// 0.8° below the horizon for radio; ITU-R P.834 gives 0.76° there.
    ///
    /// ```
    /// use geo3d::{Geodetic, STANDARD_REFRACTION_K, WGS84};
    /// // Two 100 m masts 76 km apart: beyond the straight-line horizon, within the radio one.
    /// let a = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 100.0));
    /// let b = WGS84.to_ecef(Geodetic::new(0.0, 0.683, 100.0));
    /// assert!(!WGS84.line_of_sight(a, b));
    /// assert!(WGS84.line_of_sight_with_refraction(a, b, STANDARD_REFRACTION_K));
    ///
    /// // A geostationary satellite 18° below the horizon stays out of sight.
    /// let ground = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 0.0));
    /// let satellite = WGS84.to_ecef(Geodetic::new(0.0, 100.0, 35_786_000.0));
    /// assert!(!WGS84.line_of_sight_with_refraction(ground, satellite, STANDARD_REFRACTION_K));
    /// ```
    pub fn line_of_sight_with_refraction(self, a: Ecef, b: Ecef, k: f64) -> bool {
        let d = b.vec() - a.vec();
        let height = |t: f64| self.to_geodetic(a.offset(d * t)).height_m;
        // The path is pinned at both ends and curves with κ = (1 − 1/k)/R along the stretch [p, q]
        // that is inside the bending layer, so its rise above the straight line is
        // κ·L²·∫ G(t, u) du over [p, q], with G(t, u) = u(1 − t) for u ≤ t and t(1 − u) for u ≥ t.
        // Over the whole path that is the classic κ·L²·t(1 − t)/2.
        let kl2 = (1.0 - 1.0 / k) * d.dot(d) / self.mean_radius();
        let rise = |t: f64, (p, q): (f64, f64)| {
            if p >= q {
                return 0.0;
            }
            let m = t.clamp(p, q);
            kl2 * ((1.0 - t) * (m * m - p * p) / 2.0 + t * ((q - m) - (q * q - m * m) / 2.0))
        };
        // Where the straight line is in the layer...
        let layer = Ellipsoid::from_axes(self.a() + REFRACTION_CEILING_M, self.b() + REFRACTION_CEILING_M);
        let mut span = match layer.ray_roots(a, d) {
            Some((t1, t2)) => (t1.max(0.0), t2.min(1.0)),
            None => (1.0, 0.0),
        };
        // ...then where the bent path is, which leaves the layer sooner. A path that stays in the
        // air throughout (the usual radio link) is done already.
        if span.0 < span.1 && (span.0 > 0.0 || span.1 < 1.0) {
            for _ in 0..4 {
                let above = |t: f64| height(t) + rise(t, span) >= REFRACTION_CEILING_M;
                let edge = |inside: f64, outside: f64| {
                    let (mut i, mut o) = (inside, outside);
                    for _ in 0..40 {
                        let mid = 0.5 * (i + o);
                        if above(mid) {
                            o = mid;
                        } else {
                            i = mid;
                        }
                    }
                    0.5 * (i + o)
                };
                let mid = 0.5 * (span.0 + span.1);
                let p = if span.0 > 0.0 { edge(mid, 0.0) } else { 0.0 };
                let q = if span.1 < 1.0 { edge(mid, 1.0) } else { 1.0 };
                span = (p, q);
            }
        }
        let clearance = |t: f64| height(t) + rise(t, span);
        // The clearance can have low points both near the ends and in the middle, so sample the
        // path first, then close in on the lowest sample with a golden-section search.
        const SAMPLES: usize = 32;
        let mut lowest = f64::INFINITY;
        let mut at = 0;
        for i in 0..=SAMPLES {
            let f = clearance(i as f64 / SAMPLES as f64);
            if f.is_nan() {
                return false;
            }
            if f < lowest {
                (lowest, at) = (f, i);
            }
        }
        const G: f64 = 0.618_033_988_749_894_8; // (√5 − 1) / 2
        let (mut lo, mut hi) = (
            at.saturating_sub(1) as f64 / SAMPLES as f64,
            ((at + 1).min(SAMPLES)) as f64 / SAMPLES as f64,
        );
        let (mut x1, mut x2) = (hi - G * (hi - lo), lo + G * (hi - lo));
        let (mut f1, mut f2) = (clearance(x1), clearance(x2));
        for _ in 0..48 {
            if f1 <= f2 {
                (hi, x2, f2) = (x2, x1, f1);
                x1 = hi - G * (hi - lo);
                f1 = clearance(x1);
            } else {
                (lo, x1, f1) = (x1, x2, f2);
                x2 = lo + G * (hi - lo);
                f2 = clearance(x2);
            }
        }
        lowest.min(f1).min(f2) >= -1e-6
    }
}

/// Where the air stops bending radio and light: the refractivity scale height of ITU-R P.453's
/// reference atmosphere, metres.
const REFRACTION_CEILING_M: f64 = 7_350.0;

#[cfg(test)]
mod tests {
    use super::STANDARD_REFRACTION_K;
    use crate::test_util::close;
    use crate::{Ecef, Ellipsoid, Geodetic, Vec3, WGS84};

    /// Rays from geostationary distance that skim the Earth: the near root matches one refined in
    /// 32-digit arithmetic to a few micrometres along the ray. (In plain f64 it was 0.5 mm out
    /// skimming by 1e-9, and 7.5 mm by 1e-12.)
    #[test]
    fn grazing_rays_from_far_away_keep_their_precision() {
        use crate::vector::Dd;
        let ell = Ellipsoid::new(1.0, 0.0); // the unit sphere: the roots solve |o + t·d|² = 1
        let o = Vec3::new(6.6, 0.0, 0.0); // 42 000 km out, in Earth radii
        for skim in [1e-3, 1e-6, 1e-9, 1e-12] {
            // Aim so the ray passes 1 − skim from the centre: just inside the surface.
            let miss = 1.0 - skim;
            let d = Vec3::new(-(o.x * o.x - miss * miss).sqrt(), miss, 0.0) * (1.0 / o.x);
            let (t, _) = ell.ray_roots(Ecef::new(o.x, o.y, o.z), d).unwrap();
            // Newton on |o + t·d|² − 1, evaluated in double-double.
            let mut exact = t;
            for _ in 0..3 {
                let p =
                    [o.x, o.y, o.z]
                        .map(Dd::from)
                        .iter()
                        .zip([d.x, d.y, d.z])
                        .fold(Dd(-1.0, 0.0), |acc, (&oi, di)| {
                            let x = oi + Dd::prod(exact, di);
                            acc + x * x
                        });
                let slope = 2.0 * (o + d * exact).dot(d);
                exact -= p.f64() / slope;
            }
            // 1e-13 of the distance is about 4 µm.
            assert!((t - exact).abs() < 1e-13 * exact, "skim {skim}: {t} vs {exact}");
        }
    }

    #[test]
    fn ray_roots_have_no_cancellation() {
        // Vieta: the roots' product is c/a. The textbook formula loses the small root to
        // cancellation when the origin is near the surface; the stable form keeps it, and c is
        // worked out without cancelling too.
        let ell = Ellipsoid::new(1.0, 0.0);
        let o = Ecef::new(1.0 + 1e-9, 0.0, 0.0);
        for deg in [10.0f64, 60.0, 89.0, 120.0, 179.0] {
            let d = Vec3::new(deg.to_radians().cos(), deg.to_radians().sin(), 0.0);
            let (t1, t2) = ell.ray_roots(o, d).unwrap();
            let c = (o.x - 1.0) * (o.x + 1.0); // |o|² − 1, without cancelling
            assert!(close(t1 * t2 / c, 1.0, 1e-13), "{deg}°: product {} vs {c}", t1 * t2);
            assert!(t1 <= t2);
        }
        // Tangent from a point on the surface: a double root at 0.
        assert_eq!(
            ell.ray_roots(Ecef::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)),
            Some((0.0, 0.0))
        );
    }

    #[test]
    fn ray_intersect_finds_the_near_surface_or_misses() {
        // Straight down from 10 km over the origin hits (a, 0, 0).
        let o = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 10_000.0));
        let hit = WGS84.ray_intersect(o, Vec3::new(-1.0, 0.0, 0.0)).unwrap();
        assert!(close(hit.x, WGS84.a(), 1e-6) && close(hit.y, 0.0, 1e-9));
        // Pointing away: the surface is behind the ray.
        assert!(WGS84.ray_intersect(o, Vec3::new(1.0, 0.0, 0.0)).is_none());
        // Sideways from far out in space: misses.
        assert!(WGS84
            .ray_intersect(Ecef::new(1e8, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0))
            .is_none());
        // Zero direction: no ray.
        assert!(WGS84.ray_intersect(o, Vec3::ZERO).is_none());
    }

    /// Masts in the air keep the classic effective-Earth answer: visible just inside the radio
    /// horizon, not just beyond it.
    #[test]
    fn links_through_the_air_keep_the_classic_model() {
        let r = WGS84.a(); // along the equator the curvature radius is a
        for (h1, h2) in [(10.0, 10.0), (100.0, 50.0), (30.0, 3_000.0), (500.0, 7_000.0)] {
            let reach = (2.0 * STANDARD_REFRACTION_K * r * h1).sqrt() + (2.0 * STANDARD_REFRACTION_K * r * h2).sqrt();
            for (factor, visible) in [(0.995, true), (1.005, false)] {
                let lon = (reach * factor / r).to_degrees();
                let a = WGS84.to_ecef(Geodetic::new(0.0, 0.0, h1));
                let b = WGS84.to_ecef(Geodetic::new(0.0, lon, h2));
                assert_eq!(
                    WGS84.line_of_sight_with_refraction(a, b, STANDARD_REFRACTION_K),
                    visible,
                    "{h1}, {h2}, {factor}"
                );
            }
        }
    }

    /// Satellites: radio refraction lets a ground station see a distant satellite about 0.8° below
    /// the horizon, near ITU-R P.834's 0.76°.
    #[test]
    fn satellites_are_seen_a_little_below_the_horizon() {
        let station = crate::LocalFrame::new(Geodetic::new(-27.5, 153.0, 0.0));
        for range_m in [2.0e6, 4.0e6, 4.0e7] {
            // (Nearer than this, part of the path is above the air, and the bending is less.)
            let visible = |el: f64| {
                let sat = station.ecef_of_aer(crate::Aer::new(40.0, el, range_m));
                WGS84.line_of_sight_with_refraction(station.origin_ecef(), sat, STANDARD_REFRACTION_K)
            };
            // Find where it stops being visible, to a hundredth of a degree.
            let (mut lo, mut hi) = (-3.0, 0.0);
            while hi - lo > 0.01 {
                let mid = 0.5 * (lo + hi);
                if visible(mid) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            assert!((-0.85..=-0.7).contains(&hi), "{range_m} m: seen down to {hi}°");
            assert!(!visible(-5.0) && visible(0.5));
        }
    }

    #[test]
    fn refraction_horizon_is_consistent() {
        let k = STANDARD_REFRACTION_K;
        let r = WGS84.mean_radius();
        // Classic below the bending layer, continuous at its top, straight with k = 1.
        for h in [1.0, 100.0, 7_000.0, 7_350.0] {
            let classic = (2.0 * k * r * h + h * h).sqrt();
            assert!(close(WGS84.horizon_distance_with_refraction_m(h, k), classic, 1e-6));
        }
        let (below, above) = (
            WGS84.horizon_distance_with_refraction_m(7_350.0 - 1e-6, k),
            WGS84.horizon_distance_with_refraction_m(7_350.0 + 1e-6, k),
        );
        assert!((above - below).abs() < 1e-3);
        for h in [10.0, 1e4, 4e5, 3.6e7] {
            let straight = WGS84.horizon_distance_m(h);
            assert!(close(
                WGS84.horizon_distance_with_refraction_m(h, 1.0),
                straight,
                1e-6 * straight
            ));
            // Refraction helps, but above the air only by the bending's worth: about 0.8° of arc,
            // under 100 km.
            let bent = WGS84.horizon_distance_with_refraction_m(h, k);
            assert!(
                bent > straight && (h < 7_350.0 || bent - straight < 100_000.0),
                "{h}: {bent} vs {straight}"
            );
        }
    }

    #[test]
    fn refraction_reaches_further_and_k_1_is_straight() {
        // With k = 1 the refracted test agrees with the straight one, pair by pair.
        let heights = [0.0, 10.0, 1_000.0, 10_000.0, 400_000.0];
        for (i, &ha) in heights.iter().enumerate() {
            for &hb in &heights {
                for dlon in [0.01, 0.3, 1.0, 3.0, 10.0, 40.0, 120.0] {
                    let a = WGS84.to_ecef(Geodetic::new(5.0 * i as f64, 0.0, ha));
                    let b = WGS84.to_ecef(Geodetic::new(-3.0, dlon, hb));
                    assert_eq!(
                        WGS84.line_of_sight_with_refraction(a, b, 1.0),
                        WGS84.line_of_sight(a, b),
                        "{ha} m to {hb} m, {dlon}° apart"
                    );
                }
            }
        }
        // Radio reaches about 15 % further than a straight line.
        let ratio =
            WGS84.horizon_distance_with_refraction_m(100.0, STANDARD_REFRACTION_K) / WGS84.horizon_distance_m(100.0);
        assert!(close(ratio, (4.0f64 / 3.0).sqrt(), 1e-5), "{ratio}"); // √k, apart from the tiny h² term
                                                                       // NaN in: no line of sight.
        let a = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 100.0));
        assert!(!WGS84.line_of_sight_with_refraction(a, Ecef::new(f64::NAN, 0.0, 0.0), 1.2));
    }

    #[test]
    fn horizon_and_line_of_sight() {
        // ~357 km from 10 km up; zero from the surface; negative heights clamp to zero.
        assert!(close(WGS84.horizon_distance_m(10_000.0), 357_000.0, 1_000.0));
        assert_eq!(WGS84.horizon_distance_m(0.0), 0.0);
        assert_eq!(WGS84.horizon_distance_m(-5.0), 0.0);
        let a = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 10_000.0));
        let near = WGS84.to_ecef(Geodetic::new(0.0, 1.0, 10_000.0)); // 111 km: well within the horizon
        let far = WGS84.to_ecef(Geodetic::new(0.0, 10.0, 10_000.0)); // 1113 km: over the horizon
        assert!(WGS84.line_of_sight(a, near) && WGS84.line_of_sight(near, a));
        assert!(!WGS84.line_of_sight(a, far) && !WGS84.line_of_sight(far, a));
        assert!(WGS84.line_of_sight(a, a));
        // Two satellites on opposite sides of the Earth are blocked; two close together are not.
        let s1 = Ecef::new(2.0e7, 0.0, 0.0);
        assert!(!WGS84.line_of_sight(s1, Ecef::new(-2.0e7, 0.0, 0.0)));
        assert!(WGS84.line_of_sight(s1, Ecef::new(0.0, 2.0e7, 0.0)));
        // Two points on the ground can't see each other: the Earth bulges between them.
        let g1 = WGS84.to_ecef(Geodetic::new(10.0, 20.0, 0.0));
        for dlon in [0.01, 1.0] {
            let g2 = WGS84.to_ecef(Geodetic::new(10.0, 20.0 + dlon, 0.0));
            assert!(!WGS84.line_of_sight(g1, g2), "{dlon}°");
        }
        // From the surface, a point 1 km up 5 km away is visible.
        let ground = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 0.0));
        assert!(WGS84.line_of_sight(ground, WGS84.to_ecef(Geodetic::new(0.0, 0.045, 1000.0))));
    }
}
