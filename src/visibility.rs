//! Rays against the ellipsoid: where a line meets the surface, line of sight and horizon distance.
//! The surface is smooth — no terrain, no atmospheric refraction.

use crate::ellipsoid::Ellipsoid;
use crate::position::Ecef;
use crate::vector::Vec3;

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
    /// becomes the unit sphere and solves the ray–sphere quadratic.
    pub(crate) fn ray_roots(self, origin: Ecef, dir: Vec3) -> Option<(f64, f64)> {
        let s = Vec3::new(1.0 / self.a(), 1.0 / self.a(), 1.0 / self.b());
        let o = Vec3::new(origin.x * s.x, origin.y * s.y, origin.z * s.z);
        let d = Vec3::new(dir.x * s.x, dir.y * s.y, dir.z * s.z);
        let (qa, qb, qc) = (d.dot(d), 2.0 * o.dot(d), o.dot(o) - 1.0);
        let disc = qb * qb - 4.0 * qa * qc;
        if disc.is_nan() || disc < 0.0 || qa == 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        Some(((-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa)))
    }

    /// Whether the straight line between two points clears the ellipsoid — can `a` see `b`? Both
    /// points should be at or above the surface; a point on the surface sees another only if the
    /// line does not dip below it.
    pub fn line_of_sight(self, a: Ecef, b: Ecef) -> bool {
        let d = b.vec() - a.vec();
        let dd = d.dot(d);
        if dd == 0.0 {
            return true;
        }
        // Reuse the ray solve but reject a hit that is just numerical noise at `a` itself.
        match self.ray_intersect(a, d) {
            None => true,
            Some(hit) => {
                let t = (hit.vec() - a.vec()).dot(d) / dd;
                t >= 1.0 - 1e-12 || t <= 1e-9
            }
        }
    }

    /// How far you can see from `height_m` above the surface before the Earth curves away, metres
    /// (on a sphere of the mean radius). Negative heights count as zero.
    #[inline]
    pub fn horizon_distance_m(self, height_m: f64) -> f64 {
        let r = self.mean_radius();
        let h = height_m.max(0.0);
        (2.0 * r * h + h * h).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use crate::test_util::close;
    use crate::{Ecef, Geodetic, Vec3, WGS84};

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
        // From the surface, a point 1 km up 5 km away is visible.
        let ground = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 0.0));
        assert!(WGS84.line_of_sight(ground, WGS84.to_ecef(Geodetic::new(0.0, 0.045, 1000.0))));
    }
}
