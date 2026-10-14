//! The Earth's shape — a reference [`Ellipsoid`] — and the conversions between latitude/longitude/height
//! and x/y/z on it.

use crate::angle::atan2d;
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::{Ecef, Geodetic};
use crate::vector::{fma, Vec3};

/// A reference ellipsoid: the smooth, slightly flattened sphere that latitude, longitude and height
/// are measured on. Defined by its equatorial radius `a` (metres) and flattening `f`.
///
/// Use a named one ([`WGS84`] — what GPS uses — or [`GRS80`], …) or build your own with
/// [`Ellipsoid::new`]. All the calculations in this crate are methods on it.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ellipsoid {
    a: f64,
    f: f64,
}

/// WGS 84 — GPS, and the frame CesiumJS renders in. The default.
pub const WGS84: Ellipsoid = Ellipsoid::new(6_378_137.0, 1.0 / 298.257_223_563);
/// GRS 80 — under most modern national datums (NAD83, GDA94/2020, ETRS89).
pub const GRS80: Ellipsoid = Ellipsoid::new(6_378_137.0, 1.0 / 298.257_222_101);
/// WGS 72.
pub const WGS72: Ellipsoid = Ellipsoid::new(6_378_135.0, 1.0 / 298.26);
/// International 1924 (Hayford).
pub const INTERNATIONAL_1924: Ellipsoid = Ellipsoid::new(6_378_388.0, 1.0 / 297.0);
/// Airy 1830 (OSGB36, Great Britain).
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
    /// An ellipsoid from its equatorial radius (metres) and flattening.
    #[inline]
    pub const fn new(a: f64, f: f64) -> Self {
        Ellipsoid { a, f }
    }

    /// Equatorial radius (semi-major axis), metres.
    #[inline]
    pub const fn a(self) -> f64 {
        self.a
    }

    /// Flattening `(a - b) / a`.
    #[inline]
    pub const fn f(self) -> f64 {
        self.f
    }

    /// Polar radius (semi-minor axis), metres.
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

    /// Latitude, longitude and height → x, y, z. Exact.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// let e = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 0.0));
    /// assert_eq!((e.x, e.y, e.z), (6_378_137.0, 0.0, 0.0)); // on the equator at longitude 0
    /// ```
    #[doc(alias = "geodetic2ecef")]
    #[doc(alias = "lla2ecef")]
    #[inline]
    pub fn to_ecef(self, g: Geodetic) -> Ecef {
        let (slat, clat) = g.lat_rad().sin_cos();
        let (slon, clon) = g.lon_rad().sin_cos();
        let e2 = self.e2();
        let n = self.a / fma(e2, -(slat * slat), 1.0).sqrt();
        let nh = n + g.height_m;
        Ecef::new(nh * clat * clon, nh * clat * slon, fma(n, 1.0 - e2, g.height_m) * slat)
    }

    /// x, y, z → latitude, longitude and height. The exact inverse of [`to_ecef`](Self::to_ecef).
    ///
    /// The answer is within 3 nm of the true position from deep inside the Earth up to 1 000 km,
    /// and within 20 nm out to geostationary height. That is the limit of 64-bit floating point.
    /// It handles the polar axis, the equatorial plane and the Earth's centre exactly. The centre
    /// gives latitude 90° and height −b.
    ///
    /// It uses a closed form with no iteration (Vermeille 2002, refined by Karney for
    /// GeographicLib), followed by one correction step. The step takes out the last nanometre of
    /// rounding.
    ///
    /// ```
    /// use geo3d::{Ecef, Geodetic, WGS84};
    /// let g = WGS84.to_geodetic(WGS84.to_ecef(Geodetic::new(-33.87, 151.21, 20.0)));
    /// assert!((g.lat_deg + 33.87).abs() < 1e-12 && (g.height_m - 20.0).abs() < 1e-8);
    ///
    /// let centre = WGS84.to_geodetic(Ecef::new(0.0, 0.0, 0.0));
    /// assert_eq!((centre.lat_deg, centre.height_m), (90.0, -WGS84.b()));
    /// ```
    #[doc(alias = "ecef2geodetic")]
    #[doc(alias = "ecef2lla")]
    pub fn to_geodetic(self, e: Ecef) -> Geodetic {
        let (a, f) = (self.a, self.f);
        let e2 = self.e2();
        let e2m = (1.0 - f) * (1.0 - f);
        let (e2a, e4a) = (e2.abs(), e2 * e2);
        let (x, y, z) = (e.x, e.y, e.z);
        // √(x² + y²) rather than hypot: as accurate for any distance in metres and much cheaper.
        // Only coordinates beyond 1e154 m overflow, and those take the far-away branch below.
        let r = fma(x, x, y * y).sqrt();
        let (slam, clam) = if r != 0.0 { (y / r, x / r) } else { (0.0, 1.0) };
        let mut h = fma(r, r, z * z).sqrt(); // distance from the centre
        let (sphi, cphi);
        if h.is_nan() || h > 2.0 * a / f64::EPSILON {
            // Absurdly far away: the Earth is a point and `h` is already the height.
            let rr = (x / 2.0).hypot(y / 2.0);
            let hh = (z / 2.0).hypot(rr);
            (sphi, cphi) = (z / 2.0 / hh, rr / hh);
        } else if e4a == 0.0 {
            // A sphere. The centre maps to the north pole, as for an ellipsoid.
            let zz = if h == 0.0 { 1.0 } else { z };
            let hh = zz.hypot(r);
            (sphi, cphi) = (zz / hh, r / hh);
            h -= a;
        } else {
            // Solve Vermeille's quartic for k; a prolate ellipsoid swaps the roles of R and Z.
            let (mut p, mut q) = ((r / a) * (r / a), e2m * (z / a) * (z / a));
            let rr = (p + q - e4a) / 6.0;
            if f < 0.0 {
                core::mem::swap(&mut p, &mut q);
            }
            if !(e4a * q == 0.0 && rr <= 0.0) {
                let s = e4a * p * q / 4.0;
                let r2 = rr * rr;
                let r3 = rr * r2;
                let disc = s * (2.0 * r3 + s);
                let mut u = rr;
                if disc >= 0.0 {
                    // Pick the sign of the square root that avoids cancellation.
                    let t3 = s + r3 + disc.sqrt().copysign(s + r3);
                    let t = t3.cbrt();
                    u += t + if t != 0.0 { r2 / t } else { 0.0 };
                } else {
                    let ang = (-disc).sqrt().atan2(-(s + r3));
                    u += 2.0 * rr * (ang / 3.0).cos();
                }
                let v = (u * u + e4a * q).sqrt();
                let uv = if u < 0.0 { e4a * q / (v - u) } else { u + v };
                let w = (e2a * (uv - q) / (2.0 * v)).max(0.0);
                let k = uv / ((uv + w * w).sqrt() + w);
                let (k1, k2) = if f >= 0.0 { (k, k + e2) } else { (k - e2, k) };
                let d = k1 * r / k2;
                // atan2 needs no normalising, and skipping it saves two roundings.
                (sphi, cphi) = (z / k1, r / k2);
                h = (1.0 - e2m / k1) * fma(d, d, z * z).sqrt();
            } else {
                // In the equatorial plane, inside the evolute (near the centre): the general
                // formulas give 0/0 here, so take their limits.
                let zz = ((if f >= 0.0 { e4a - p } else { p }) / e2m).sqrt();
                let xx = (if f < 0.0 { e4a - p } else { p }).sqrt();
                let hh = zz.hypot(xx);
                (sphi, cphi) = ((zz / hh).copysign(if z < 0.0 { -1.0 } else { 1.0 }), xx / hh);
                h = -a * (if f >= 0.0 { e2m } else { 1.0 }) * hh / e2a;
            }
        }
        // Longitude straight from x and y: normalising them first would add two roundings.
        let lon = if r != 0.0 { atan2d(y, x) } else { atan2d(slam, clam) };
        let mut lat = atan2d(sphi, cphi);
        if h > -1.0e6 && h <= 2.0 * a / f64::EPSILON {
            // One correction step takes the last nanometre of rounding out. Go back to the
            // meridian plane, measure the leftover (dp, dz), and correct the latitude and height
            // along the local north and up. Skipped deep inside the Earth, where the closed form
            // is already exact and the step could divide by nearly zero.
            let (s, c) = lat.to_radians().sin_cos();
            let w2 = 1.0 - e2 * s * s;
            let n = a / w2.sqrt();
            let (dp, dz) = (r - (n + h) * c, z - (n * (1.0 - e2) + h) * s);
            let m = n * (1.0 - e2) / w2; // meridional radius of curvature
            lat += ((c * dz - s * dp) / (m + h)).to_degrees();
            h += c * dp + s * dz;
        }
        Geodetic::new(lat, lon, h)
    }

    /// Straight-line interpolation between two positions through space, `t = 0` at `a` and `1` at
    /// `b` — for dead reckoning between two fixes a short time apart. Over long distances the
    /// straight line dips below the surface; use
    /// [`great_circle_destination`](Self::great_circle_destination) for those.
    pub fn interpolate(self, a: Geodetic, b: Geodetic, t: f64) -> Geodetic {
        let (pa, pb) = (self.to_ecef(a).vec(), self.to_ecef(b).vec());
        self.to_geodetic(Ecef::from(pa + (pb - pa) * t))
    }

    /// The outward unit normal to the surface at a position — the local "up".
    #[inline]
    pub fn surface_normal(self, g: Geodetic) -> Vec3 {
        let (slat, clat) = g.lat_rad().sin_cos();
        let (slon, clon) = g.lon_rad().sin_cos();
        Vec3::new(clat * clon, clat * slon, slat)
    }

    /// The mean radius `(2a + b) / 3`, metres — the sphere the great-circle helpers use.
    #[inline]
    pub const fn mean_radius(self) -> f64 {
        (2.0 * self.a + self.b()) / 3.0
    }

    /// Distance from the Earth's centre to the surface at a latitude (radians), metres.
    pub fn geocentric_radius(self, lat_rad: f64) -> f64 {
        let (s, c) = lat_rad.sin_cos();
        let (a, b) = (self.a, self.b());
        let (ac, bs) = (a * c, b * s);
        ((a * ac * a * ac + b * bs * b * bs) / (ac * ac + bs * bs)).sqrt()
    }

    /// East–west radius of curvature `N` at a latitude (radians), metres.
    #[inline]
    pub fn prime_vertical_radius(self, lat_rad: f64) -> f64 {
        let s = lat_rad.sin();
        self.a / (1.0 - self.e2() * s * s).sqrt()
    }

    /// North–south radius of curvature `M` at a latitude (radians), metres.
    #[inline]
    pub fn meridional_radius(self, lat_rad: f64) -> f64 {
        let s = lat_rad.sin();
        let w = 1.0 - self.e2() * s * s;
        self.a * (1.0 - self.e2()) / (w * w.sqrt())
    }

    /// Geocentric latitude (the angle at the Earth's centre) for a geodetic latitude on the
    /// surface; both radians.
    #[inline]
    pub fn geocentric_latitude(self, geodetic_lat_rad: f64) -> f64 {
        ((1.0 - self.e2()) * geodetic_lat_rad.tan()).atan()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::close;

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
                    close(r.lat_deg, g.lat_deg, 1e-11)
                        && close(r.lon_deg, g.lon_deg, 1e-11)
                        && close(r.height_m, g.height_m, 1e-6),
                    "{ell:?} {g:?} -> {r:?}"
                );
            }
        }
    }

    #[test]
    fn dense_global_grid_round_trips_to_machine_precision() {
        let mut lat = -89.5;
        while lat <= 89.5 {
            let mut lon = -180.0;
            while lon < 180.0 {
                for h in [-11_000.0, 0.0, 8_848.0, 12_000.0, 400_000.0, 20_200_000.0, 35_786_000.0] {
                    let g = Geodetic::new(lat, lon, h);
                    let r = WGS84.to_geodetic(WGS84.to_ecef(g));
                    assert!(
                        close(r.lat_deg, lat, 1e-11) && close(r.lon_deg, lon, 1e-11) && close(r.height_m, h, 1e-5),
                        "{g:?} -> {r:?}"
                    );
                }
                lon += 15.0;
            }
            lat += 7.0;
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
    fn deep_inside_and_far_out_round_trip() {
        // The iterated Bowring method used before 0.3.1 was up to 53 m out deep inside the Earth.
        for lat in [-89.0, -60.0, -30.0, -1.0, 0.0, 1.0, 30.0, 60.0, 89.0] {
            for h in [-6.0e6, -3.0e6, -1.0e6, 1.0e8, 1.0e9] {
                let g = Geodetic::new(lat, 40.0, h);
                let r = WGS84.to_geodetic(WGS84.to_ecef(g));
                // Tolerances scale with the distance from the centre: the rounding of x, y, z.
                let tol_m = 4.0 * f64::EPSILON * (WGS84.a() + h.abs());
                let tol_deg = (tol_m / (WGS84.a() + h).abs().max(1.0e5)).to_degrees();
                assert!(
                    close(r.lat_deg, lat, tol_deg) && close(r.height_m, h, tol_m),
                    "{g:?} -> {r:?}"
                );
            }
        }
    }

    #[test]
    fn inside_the_evolute_matches_the_exact_answer() {
        // Near the centre, in the equatorial plane, the nearest surface point is not on the
        // equator. A point at R = a·e²·cos φ / W (W = √(1 − e² sin² φ)) has latitude φ and height
        // −a(1 − e²)/W, exactly.
        for lat in [10.0f64, 30.0, 60.0, 85.0] {
            let (s, c) = lat.to_radians().sin_cos();
            let w = (1.0 - WGS84.e2() * s * s).sqrt();
            let r = WGS84.a() * WGS84.e2() * c / w;
            let g = WGS84.to_geodetic(Ecef::new(r, 0.0, 0.0));
            assert!(close(g.lat_deg, lat, 1e-10), "{lat}: {g:?}");
            assert!(
                close(g.height_m, -WGS84.a() * (1.0 - WGS84.e2()) / w, 1e-6),
                "{lat}: {g:?}"
            );
        }
        // The centre itself: latitude 90°, height −b.
        let c = WGS84.to_geodetic(Ecef::new(0.0, 0.0, 0.0));
        assert_eq!((c.lat_deg, c.lon_deg, c.height_m), (90.0, 0.0, -WGS84.b()));
    }

    #[test]
    fn spheres_prolate_ellipsoids_and_extremes() {
        let sphere = Ellipsoid::new(6_371_000.0, 0.0);
        let prolate = Ellipsoid::new(6_378_137.0, -1.0 / 150.0);
        for ell in [sphere, prolate] {
            for (lat, h) in [(0.0, 0.0), (45.0, 1_000.0), (-89.5, 20_000.0), (60.0, -500_000.0)] {
                let g = Geodetic::new(lat, -120.0, h);
                let r = ell.to_geodetic(ell.to_ecef(g));
                assert!(
                    close(r.lat_deg, lat, 1e-11) && close(r.height_m, h, 1e-6),
                    "{ell:?} {g:?} -> {r:?}"
                );
            }
        }
        // Absurdly far away, the Earth is a point: the direction is kept, and the height is the
        // distance from the centre.
        let far = WGS84.to_geodetic(Ecef::new(3.0e30, 0.0, 4.0e30));
        assert!(close(far.lat_deg, 53.130_102_354_155_98, 1e-12) && close(far.height_m, 5.0e30, 1e15));
        // NaN in gives NaN out.
        assert!(WGS84.to_geodetic(Ecef::new(f64::NAN, 0.0, 0.0)).lat_deg.is_nan());
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
