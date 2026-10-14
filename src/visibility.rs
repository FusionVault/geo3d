//! Rays against the ellipsoid: where a line meets the surface, line of sight and horizon distance.
//! The surface is smooth, with no terrain. Refraction (the air bending the line of sight) is
//! optional.

use crate::ellipsoid::Ellipsoid;
use crate::position::Ecef;
use crate::vector::Vec3;

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
        let (qa, qb, qc) = (d.dot(d), 2.0 * o.dot(d), o.dot(o) - 1.0);
        let disc = qb * qb - 4.0 * qa * qc;
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
        self.horizon_distance_with_refraction_m(height_m, 1.0)
    }

    /// Like [`horizon_distance_m`](Self::horizon_distance_m), but with the air bending the line of
    /// sight around the Earth. `k` is the effective-Earth-radius factor:
    /// [`STANDARD_REFRACTION_K`] (4/3) for radio and radar, about 7/6 for visible light, and 1 for
    /// none.
    ///
    /// ```
    /// use geo3d::{STANDARD_REFRACTION_K, WGS84};
    /// let straight = WGS84.horizon_distance_m(100.0);                                    // 35.7 km
    /// let radio = WGS84.horizon_distance_with_refraction_m(100.0, STANDARD_REFRACTION_K); // 41.2 km
    /// assert!(radio > straight * 1.15);
    /// ```
    #[inline]
    pub fn horizon_distance_with_refraction_m(self, height_m: f64, k: f64) -> f64 {
        let r = k * self.mean_radius();
        let h = height_m.max(0.0);
        (2.0 * r * h + h * h).sqrt()
    }

    /// Like [`line_of_sight`](Self::line_of_sight), but with the air bending the line of sight
    /// around the Earth, so radar and radio links see further. `k` is the effective-Earth-radius
    /// factor: [`STANDARD_REFRACTION_K`] for radio, about 7/6 for visible light, and 1 for none.
    ///
    /// It uses the standard effective-Earth model: at distances d₁ and d₂ from the two ends, the
    /// path's clearance above the surface gains `d₁·d₂·(1 − 1/k) / 2R`. The path is clear if its
    /// lowest point stays above the surface (to within a micrometre).
    ///
    /// ```
    /// use geo3d::{Geodetic, STANDARD_REFRACTION_K, WGS84};
    /// // Two 100 m masts 76 km apart: beyond the straight-line horizon, within the radio one.
    /// let a = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 100.0));
    /// let b = WGS84.to_ecef(Geodetic::new(0.0, 0.683, 100.0));
    /// assert!(!WGS84.line_of_sight(a, b));
    /// assert!(WGS84.line_of_sight_with_refraction(a, b, STANDARD_REFRACTION_K));
    /// ```
    pub fn line_of_sight_with_refraction(self, a: Ecef, b: Ecef, k: f64) -> bool {
        let d = b.vec() - a.vec();
        let bend = (1.0 - 1.0 / k) * d.dot(d) / (2.0 * self.mean_radius());
        let clearance = |t: f64| self.to_geodetic(a.offset(d * t)).height_m + bend * t * (1.0 - t);
        // The clearance along the path has a single low point, so a golden-section search finds
        // it. Sixty-four steps narrow it to 1e-13 of the path's length.
        const G: f64 = 0.618_033_988_749_894_8; // (√5 − 1) / 2
        let (mut lo, mut hi) = (0.0, 1.0);
        let (mut x1, mut x2) = (1.0 - G, G);
        let (mut f1, mut f2) = (clearance(x1), clearance(x2));
        for _ in 0..64 {
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
        let lowest = [f1, f2, clearance(0.0), clearance(1.0)]
            .into_iter()
            .fold(f64::INFINITY, |m, f| {
                if f.is_nan() || m.is_nan() {
                    f64::NAN
                } else {
                    m.min(f)
                }
            });
        lowest >= -1e-6
    }
}

#[cfg(test)]
mod tests {
    use super::STANDARD_REFRACTION_K;
    use crate::test_util::close;
    use crate::{Ecef, Ellipsoid, Geodetic, Vec3, WGS84};

    #[test]
    fn ray_roots_have_no_cancellation() {
        // Vieta: the roots' product is c/a. The textbook formula loses the small root to
        // cancellation when the origin is near the surface; the stable form keeps it.
        let ell = Ellipsoid::new(1.0, 0.0);
        let o = Ecef::new(1.0 + 1e-9, 0.0, 0.0);
        for deg in [10.0f64, 60.0, 89.0, 120.0, 179.0] {
            let d = Vec3::new(deg.to_radians().cos(), deg.to_radians().sin(), 0.0);
            let (t1, t2) = ell.ray_roots(o, d).unwrap();
            let c = o.vec().dot(o.vec()) - 1.0;
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
