//! Reference ellipsoids and the geodetic ↔ ECEF conversions on them.

use crate::types::{Ecef, Geodetic, Vec3};

/// A reference ellipsoid of revolution: semi-major axis `a` (metres) and flattening `f`.
///
/// Use the constants ([`WGS84`], [`GRS80`], …) or build your own with [`Ellipsoid::new`].
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ellipsoid {
    a: f64,
    f: f64,
}

/// WGS 84 — GPS, and the frame CesiumJS renders in.
pub const WGS84: Ellipsoid = Ellipsoid::new(6_378_137.0, 1.0 / 298.257_223_563);
/// GRS 80 — the ellipsoid under most modern national datums (NAD83, GDA94/2020, ETRS89).
pub const GRS80: Ellipsoid = Ellipsoid::new(6_378_137.0, 1.0 / 298.257_222_101);
/// WGS 72.
pub const WGS72: Ellipsoid = Ellipsoid::new(6_378_135.0, 1.0 / 298.26);
/// International 1924 (Hayford).
pub const INTERNATIONAL_1924: Ellipsoid = Ellipsoid::new(6_378_388.0, 1.0 / 297.0);
/// Airy 1830 (OSGB36).
pub const AIRY_1830: Ellipsoid = Ellipsoid::new(6_377_563.396, 1.0 / 299.324_964_6);
/// Bessel 1841.
pub const BESSEL_1841: Ellipsoid = Ellipsoid::new(6_377_397.155, 1.0 / 299.152_812_8);
/// PZ-90 (GLONASS).
pub const PZ90: Ellipsoid = Ellipsoid::new(6_378_136.0, 1.0 / 298.257_839_303);
/// Clarke 1866 (NAD27).
pub const CLARKE_1866: Ellipsoid = Ellipsoid::new(6_378_206.4, 1.0 / 294.978_698_2);

impl Default for Ellipsoid {
    /// [`WGS84`].
    fn default() -> Self {
        WGS84
    }
}

impl Ellipsoid {
    /// An ellipsoid from its semi-major axis (metres) and flattening.
    #[inline]
    pub const fn new(a: f64, f: f64) -> Self {
        Ellipsoid { a, f }
    }

    /// Semi-major (equatorial) axis, metres.
    #[inline]
    pub const fn a(self) -> f64 {
        self.a
    }

    /// Flattening `(a - b) / a`.
    #[inline]
    pub const fn f(self) -> f64 {
        self.f
    }

    /// Semi-minor (polar) axis, metres.
    #[inline]
    pub const fn b(self) -> f64 {
        self.a * (1.0 - self.f)
    }

    /// First eccentricity squared.
    #[inline]
    pub const fn e2(self) -> f64 {
        self.f * (2.0 - self.f)
    }

    /// Second eccentricity squared.
    #[inline]
    pub const fn ep2(self) -> f64 {
        let b = self.b();
        (self.a * self.a - b * b) / (b * b)
    }

    /// Prime-vertical radius of curvature `N` at a geodetic latitude (radians), metres.
    #[inline]
    pub fn prime_vertical_radius(self, lat_rad: f64) -> f64 {
        let s = lat_rad.sin();
        self.a / (1.0 - self.e2() * s * s).sqrt()
    }

    /// Geodetic → ECEF, the exact closed form.
    #[inline]
    pub fn to_ecef(self, g: Geodetic) -> Ecef {
        let (slat, clat) = g.lat_rad().sin_cos();
        let (slon, clon) = g.lon_rad().sin_cos();
        let n = self.a / (1.0 - self.e2() * slat * slat).sqrt();
        Ecef::new(
            (n + g.height_m) * clat * clon,
            (n + g.height_m) * clat * slon,
            (n * (1.0 - self.e2()) + g.height_m) * slat,
        )
    }

    /// ECEF → geodetic by Bowring's method, iterated to convergence (two or three passes; the first
    /// pass alone is only sub-millimetre near the surface, not at orbital heights). The inverse of
    /// [`Ellipsoid::to_ecef`]: round-trips to better than 1e-8° and 1e-3 m from below the surface to
    /// GNSS orbit.
    pub fn to_geodetic(self, e: Ecef) -> Geodetic {
        let (a, b, e2, ep2) = (self.a, self.b(), self.e2(), self.ep2());
        let p = (e.x * e.x + e.y * e.y).sqrt();
        let lon = e.y.atan2(e.x);
        if p < 1e-9 {
            // On the polar axis: latitude is ±90°, height is the distance from the pole.
            let lat = core::f64::consts::FRAC_PI_2.copysign(e.z);
            return Geodetic::new(lat.to_degrees(), lon.to_degrees(), e.z.abs() - b);
        }
        // Bowring: a reduced latitude θ gives a geodetic latitude φ; re-deriving θ from φ and
        // repeating converges quadratically. Five passes is a hard cap never reached in practice.
        let mut theta = (e.z * a).atan2(p * b);
        let mut lat = 0.0;
        for _ in 0..5 {
            let (st, ct) = theta.sin_cos();
            let next = (e.z + ep2 * b * st * st * st).atan2(p - e2 * a * ct * ct * ct);
            let converged = (next - lat).abs() < 1e-15;
            lat = next;
            if converged {
                break;
            }
            theta = (b * lat.tan()).atan2(a);
        }
        let (slat, clat) = lat.sin_cos();
        let n = a / (1.0 - e2 * slat * slat).sqrt();
        // Height from whichever of cos/sin is better conditioned (near the poles cos φ → 0).
        let height = if clat.abs() > 0.1 {
            p / clat - n
        } else {
            e.z / slat - n * (1.0 - e2)
        };
        Geodetic::new(lat.to_degrees(), lon.to_degrees(), height)
    }

    /// Nearest intersection of the ray `origin + t·dir` (t ≥ 0) with the ellipsoid surface, or
    /// `None` if the ray misses (points into space) or the surface is behind it. `dir` need not be
    /// unit length. Scales space to the unit sphere and solves the ray–sphere quadratic.
    pub fn ray_intersect(self, origin: Ecef, dir: Vec3) -> Option<Ecef> {
        let s = Vec3::new(1.0 / self.a, 1.0 / self.a, 1.0 / self.b());
        let o = Vec3::new(origin.x * s.x, origin.y * s.y, origin.z * s.z);
        let d = Vec3::new(dir.x * s.x, dir.y * s.y, dir.z * s.z);
        let (qa, qb, qc) = (d.dot(d), 2.0 * o.dot(d), o.dot(o) - 1.0);
        let disc = qb * qb - 4.0 * qa * qc;
        if disc < 0.0 || qa == 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        let (t1, t2) = ((-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa));
        let t = if t1 >= 0.0 {
            t1
        } else if t2 >= 0.0 {
            t2
        } else {
            return None;
        };
        Some(origin.offset(dir * t))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn derived_constants() {
        assert!(close(WGS84.b(), 6_356_752.314_245, 1e-6));
        assert!(close(WGS84.e2(), 0.006_694_379_990_14, 1e-12));
        assert_eq!(Ellipsoid::default(), WGS84);
    }

    #[test]
    fn ecef_hits_the_known_axis_points() {
        let e = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 0.0));
        assert!(close(e.x, WGS84.a(), 1e-3) && close(e.y, 0.0, 1e-6) && close(e.z, 0.0, 1e-6));
        let n = WGS84.to_ecef(Geodetic::new(90.0, 0.0, 0.0));
        assert!(close(n.x, 0.0, 1e-3) && close(n.y, 0.0, 1e-3) && close(n.z, WGS84.b(), 1e-6));
        let up = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 1000.0));
        assert!(close(up.x, WGS84.a() + 1000.0, 1e-3));
    }

    #[test]
    fn geodetic_round_trips_on_every_ellipsoid() {
        let pts = [
            Geodetic::new(0.0, 0.0, 0.0),
            Geodetic::new(34.05, -118.24, 100.0),
            Geodetic::new(-33.87, 151.21, 2000.0),
            Geodetic::new(89.9, 10.0, 500.0),
            Geodetic::new(-89.99, -179.99, 0.0),
            Geodetic::new(45.0, 90.0, 20_200_000.0), // GPS orbit height
            Geodetic::new(10.0, 20.0, -50.0),        // below the ellipsoid
        ];
        for ell in [
            WGS84,
            GRS80,
            WGS72,
            INTERNATIONAL_1924,
            AIRY_1830,
            BESSEL_1841,
            PZ90,
            CLARKE_1866,
        ] {
            for g in pts {
                let r = ell.to_geodetic(ell.to_ecef(g));
                assert!(
                    close(r.lat_deg, g.lat_deg, 1e-8)
                        && close(r.lon_deg, g.lon_deg, 1e-8)
                        && close(r.height_m, g.height_m, 1e-3),
                    "{ell:?} {g:?} -> {r:?}"
                );
            }
        }
    }

    #[test]
    fn poles_are_handled_exactly() {
        let n = WGS84.to_geodetic(Ecef::new(0.0, 0.0, WGS84.b() + 10.0));
        assert!(close(n.lat_deg, 90.0, 1e-12) && close(n.height_m, 10.0, 1e-6));
        let s = WGS84.to_geodetic(Ecef::new(0.0, 0.0, -WGS84.b()));
        assert!(close(s.lat_deg, -90.0, 1e-12) && close(s.height_m, 0.0, 1e-6));
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
}
