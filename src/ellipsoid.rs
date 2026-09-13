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

// ── 0.2 additions: radii, latitudes, normals, line of sight, great circles ───────────────────────────

impl Ellipsoid {
    /// Meridional radius of curvature `M` at a geodetic latitude (radians), metres — the radius of
    /// the north–south curve.
    #[inline]
    pub fn meridional_radius(self, lat_rad: f64) -> f64 {
        let s = lat_rad.sin();
        let w = 1.0 - self.e2() * s * s;
        self.a * (1.0 - self.e2()) / (w * w.sqrt())
    }

    /// The mean radius `(2a + b) / 3`, metres — the sphere used by the great-circle helpers.
    #[inline]
    pub const fn mean_radius(self) -> f64 {
        (2.0 * self.a + self.b()) / 3.0
    }

    /// Distance from the centre to the surface at a geodetic latitude (radians), metres.
    pub fn geocentric_radius(self, lat_rad: f64) -> f64 {
        let (s, c) = lat_rad.sin_cos();
        let (a, b) = (self.a, self.b());
        let (ac, bs) = (a * c, b * s);
        ((a * ac * a * ac + b * bs * b * bs) / (ac * ac + bs * bs)).sqrt()
    }

    /// Geocentric latitude (angle from the equatorial plane to the centre–point line) for a geodetic
    /// latitude on the surface; both radians.
    #[inline]
    pub fn geocentric_latitude(self, geodetic_lat_rad: f64) -> f64 {
        ((1.0 - self.e2()) * geodetic_lat_rad.tan()).atan()
    }

    /// The outward unit normal to the ellipsoid at a geodetic position (the local "up").
    #[inline]
    pub fn surface_normal(self, g: Geodetic) -> Vec3 {
        let (slat, clat) = g.lat_rad().sin_cos();
        let (slon, clon) = g.lon_rad().sin_cos();
        Vec3::new(clat * clon, clat * slon, slat)
    }

    /// Geometric horizon distance from a height above the surface, metres — how far a straight line
    /// of sight reaches over a smooth ellipsoid (mean radius), ignoring refraction and terrain.
    #[inline]
    pub fn horizon_distance_m(self, height_m: f64) -> f64 {
        let r = self.mean_radius();
        let h = height_m.max(0.0);
        (2.0 * r * h + h * h).sqrt()
    }

    /// Whether the straight segment between two ECEF points clears the ellipsoid surface — a
    /// geometric line-of-sight test (no terrain, no refraction). Both points must be at or above
    /// the surface; a point on the surface sees another only if the segment does not dip below it.
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

    /// Great-circle (spherical, mean radius) distance between two positions, metres, ignoring height.
    /// For sub-metre accuracy over long lines use a geodesic library; this is the haversine.
    pub fn great_circle_distance_m(self, a: Geodetic, b: Geodetic) -> f64 {
        let (p1, p2) = (a.lat_rad(), b.lat_rad());
        let dp = p2 - p1;
        let dl = b.lon_rad() - a.lon_rad();
        let h = (dp / 2.0).sin().powi(2) + p1.cos() * p2.cos() * (dl / 2.0).sin().powi(2);
        2.0 * self.mean_radius() * h.sqrt().asin()
    }

    /// Initial bearing (forward azimuth) from `a` to `b` along the great circle, degrees from north
    /// clockwise in `[0, 360)`.
    pub fn initial_bearing_deg(self, a: Geodetic, b: Geodetic) -> f64 {
        let (p1, p2) = (a.lat_rad(), b.lat_rad());
        let dl = b.lon_rad() - a.lon_rad();
        let y = dl.sin() * p2.cos();
        let x = p1.cos() * p2.sin() - p1.sin() * p2.cos() * dl.cos();
        y.atan2(x).to_degrees().rem_euclid(360.0)
    }

    /// The position `distance_m` along the great circle from `from` on an initial `bearing_deg`;
    /// height is carried over unchanged.
    pub fn great_circle_destination(self, from: Geodetic, bearing_deg: f64, distance_m: f64) -> Geodetic {
        let delta = distance_m / self.mean_radius();
        let (sd, cd) = delta.sin_cos();
        let (s1, c1) = from.lat_rad().sin_cos();
        let (sb, cb) = bearing_deg.to_radians().sin_cos();
        let s2 = s1 * cd + c1 * sd * cb;
        let lat2 = s2.asin();
        let lon2 = from.lon_rad() + (sb * sd * c1).atan2(cd - s1 * s2);
        Geodetic::new(lat2.to_degrees(), lon2.to_degrees(), from.height_m).normalized()
    }

    /// The midpoint of the great-circle arc between two positions; height is the mean of the two.
    pub fn great_circle_midpoint(self, a: Geodetic, b: Geodetic) -> Geodetic {
        let (p1, p2) = (a.lat_rad(), b.lat_rad());
        let dl = b.lon_rad() - a.lon_rad();
        let bx = p2.cos() * dl.cos();
        let by = p2.cos() * dl.sin();
        let lat = (p1.sin() + p2.sin()).atan2(((p1.cos() + bx).powi(2) + by * by).sqrt());
        let lon = a.lon_rad() + by.atan2(p1.cos() + bx);
        Geodetic::new(lat.to_degrees(), lon.to_degrees(), (a.height_m + b.height_m) / 2.0).normalized()
    }

    /// Straight-line (ECEF) interpolation between two positions, `t = 0` at `a` and `1` at `b` —
    /// dead reckoning between two fixes a short time apart. Over long arcs the chord dips below the
    /// surface; use [`great_circle_destination`](Self::great_circle_destination) for that.
    pub fn interpolate(self, a: Geodetic, b: Geodetic, t: f64) -> Geodetic {
        let (pa, pb) = (self.to_ecef(a).vec(), self.to_ecef(b).vec());
        self.to_geodetic(Ecef::from(pa + (pb - pa) * t))
    }
}

#[cfg(test)]
mod tests_0_2 {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn radii_and_latitudes() {
        assert!(close(WGS84.mean_radius(), 6_371_008.771, 1e-3));
        // M at the equator = a(1-e²), at the pole = a/sqrt(1-e²).
        assert!(close(
            WGS84.meridional_radius(0.0),
            WGS84.a() * (1.0 - WGS84.e2()),
            1e-6
        ));
        assert!(close(
            WGS84.meridional_radius(core::f64::consts::FRAC_PI_2),
            WGS84.a() / (1.0 - WGS84.e2()).sqrt(),
            1e-6
        ));
        // N at the pole equals M at the pole.
        assert!(close(
            WGS84.prime_vertical_radius(core::f64::consts::FRAC_PI_2),
            WGS84.meridional_radius(core::f64::consts::FRAC_PI_2),
            1e-6
        ));
        assert!(close(WGS84.geocentric_radius(0.0), WGS84.a(), 1e-9));
        assert!(close(
            WGS84.geocentric_radius(core::f64::consts::FRAC_PI_2),
            WGS84.b(),
            1e-6
        ));
        // Geocentric latitude is smaller in magnitude, equal at 0 and 90.
        let g45 = WGS84.geocentric_latitude(45f64.to_radians());
        assert!(g45 < 45f64.to_radians() && g45 > 44.8f64.to_radians());
        assert_eq!(WGS84.geocentric_latitude(0.0), 0.0);
        let n = WGS84.surface_normal(Geodetic::new(90.0, 0.0, 0.0));
        assert!(close(n.z, 1.0, 1e-15) && close(n.norm(), 1.0, 1e-15));
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

    #[test]
    fn great_circle_helpers() {
        let lax = Geodetic::new(33.9425, -118.4081, 0.0);
        let jfk = Geodetic::new(40.6413, -73.7781, 0.0);
        let d = WGS84.great_circle_distance_m(lax, jfk);
        assert!(close(d, 3_974_000.0, 10_000.0), "{d}");
        let brg = WGS84.initial_bearing_deg(lax, jfk);
        assert!(close(brg, 65.9, 0.5), "{brg}");
        // Destination along that bearing for that distance lands on JFK (spherical model: ~km).
        let dest = WGS84.great_circle_destination(lax, brg, d);
        assert!(
            close(dest.lat_deg, jfk.lat_deg, 0.05) && close(dest.lon_deg, jfk.lon_deg, 0.05),
            "{dest:?}"
        );
        // Midpoint is equidistant.
        let m = WGS84.great_circle_midpoint(lax, jfk);
        let (d1, d2) = (
            WGS84.great_circle_distance_m(lax, m),
            WGS84.great_circle_distance_m(m, jfk),
        );
        assert!(close(d1, d2, 1.0) && close(d1 + d2, d, 1.0));
        // Degenerate + wrap cases.
        assert_eq!(WGS84.great_circle_distance_m(lax, lax), 0.0);
        assert_eq!(
            WGS84.initial_bearing_deg(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(1.0, 0.0, 0.0)),
            0.0
        );
        assert!(close(
            WGS84.initial_bearing_deg(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, 1.0, 0.0)),
            90.0,
            1e-9
        ));
        let across = WGS84.great_circle_destination(Geodetic::new(0.0, 179.5, 7.0), 90.0, 111_000.0);
        assert!(across.lon_deg < -179.0 && across.height_m == 7.0);
    }

    #[test]
    fn interpolation_is_linear_in_ecef() {
        let a = Geodetic::new(0.0, 0.0, 1000.0);
        let b = Geodetic::new(0.0, 0.01, 3000.0);
        let m = WGS84.interpolate(a, b, 0.5);
        // With unequal heights the chord midpoint's longitude leans towards the higher endpoint
        // (~4e-7° here); with equal heights it is exact.
        assert!(close(m.lon_deg, 0.005, 1e-6));
        let level = WGS84.interpolate(a, b.with_height(1000.0), 0.5);
        assert!(close(level.lon_deg, 0.005, 1e-12), "{}", level.lon_deg);
        // The chord dips below the mean height by d²/8R ≈ 2.4 cm over 1.1 km — the documented caveat.
        assert!(m.height_m < 2000.0 && close(m.height_m, 2000.0, 0.1), "{}", m.height_m);
        assert!(close(WGS84.interpolate(a, b, 0.0).height_m, 1000.0, 1e-6));
        assert!(close(WGS84.interpolate(a, b, 1.0).lon_deg, 0.01, 1e-9));
    }
}
