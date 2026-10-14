//! Quick spherical estimates: great-circle distance, bearing, destination and midpoint on a sphere
//! of the ellipsoid's mean radius. Within 0.6 % of the true surface distance on WGS84. Plus path
//! tools built on n-vectors (the surface normal): cross-track and along-track distance, where two
//! paths cross, the centre of a set of positions, and moving along a path. Like the n-vector
//! method (Gade 2010), they have no special cases at the poles or the antimeridian.

use crate::angle::{ang_diff, atan2d};
use crate::ellipsoid::Ellipsoid;
use crate::local::wrap_360;
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;
use crate::vector::Vec3;

/// `v` scaled to unit length, or `None` if it is too short to have a direction: shorter than the
/// rounding noise in a product of vectors of length `scale`.
fn direction(v: Vec3, scale: f64) -> Option<Vec3> {
    let len = v.norm();
    (len > 4.0 * f64::EPSILON * scale).then(|| v / len)
}

/// `n(b) − n(a)`, the difference of two surface normals, built from the latitude and longitude
/// differences with half-angle identities. Subtracting the normals themselves would leave only
/// rounding noise for nearby points; this keeps full precision at any separation.
fn normal_difference(a: Geodetic, b: Geodetic) -> Vec3 {
    let (d, e) = ang_diff(a.lon_deg, b.lon_deg);
    let (dlat, dlon) = ((b.lat_deg - a.lat_deg).to_radians(), (d + e).to_radians());
    let (slat_m, clat_m) = (a.lat_rad() + dlat / 2.0).sin_cos();
    let (slon_m, clon_m) = (a.lon_rad() + dlon / 2.0).sin_cos();
    let (sh_lat, sh_lon) = ((dlat / 2.0).sin(), (dlon / 2.0).sin());
    let (slon2, clon2) = b.lon_rad().sin_cos();
    let clat1 = a.lat_rad().cos();
    // cos φ₂ − cos φ₁ = −2 sin φₘ sin(Δφ/2); (cos λ₂ − cos λ₁, sin λ₂ − sin λ₁) = 2 sin(Δλ/2)(−sin λₘ, cos λₘ)
    let dclat = -2.0 * slat_m * sh_lat;
    Vec3::new(
        dclat * clon2 - clat1 * 2.0 * sh_lon * slon_m,
        dclat * slon2 + clat1 * 2.0 * sh_lon * clon_m,
        2.0 * clat_m * sh_lat,
    )
}

/// The position whose surface normal (n-vector) is `n`, at `height_m`.
fn from_normal(n: Vec3, height_m: f64) -> Geodetic {
    Geodetic::new(atan2d(n.z, (n.x * n.x + n.y * n.y).sqrt()), atan2d(n.y, n.x), height_m)
}

impl Ellipsoid {
    /// Great-circle distance between two positions on a sphere of the mean radius, metres, ignoring
    /// height. Fast and within 0.6 % of the true surface distance on WGS84; use
    /// [`geodesic_inverse`](Self::geodesic_inverse) for full accuracy.
    ///
    /// Computed as the angle between the two positions' unit vectors, `atan2(|a×b|, a·b)`, which is
    /// accurate to about a nanometre at every distance. (The haversine formula is off by up to a
    /// metre near the antipode.)
    #[doc(alias = "haversine")]
    pub fn great_circle_distance_m(self, a: Geodetic, b: Geodetic) -> f64 {
        self.mean_radius() * self.surface_normal(a).angle_to(self.surface_normal(b))
    }

    /// Initial bearing from `a` to `b` along the great circle, degrees clockwise from north in
    /// `[0, 360)`.
    #[doc(alias = "azimuth")]
    #[doc(alias = "heading")]
    pub fn initial_bearing_deg(self, a: Geodetic, b: Geodetic) -> f64 {
        let (p1, p2) = (a.lat_rad(), b.lat_rad());
        let dl = b.lon_rad() - a.lon_rad();
        let y = dl.sin() * p2.cos();
        let x = p1.cos() * p2.sin() - p1.sin() * p2.cos() * dl.cos();
        wrap_360(y.atan2(x).to_degrees())
    }

    /// The position `distance_m` along the great circle from `from` on an initial `bearing_deg`;
    /// height is carried over.
    ///
    /// Worked with n-vectors: the start's normal, turned towards the bearing by the distance's
    /// angle. The latitude comes from `atan2`, so it stays exact right up to the poles (an `asin`
    /// form loses up to 11 cm there).
    pub fn great_circle_destination(self, from: Geodetic, bearing_deg: f64, distance_m: f64) -> Geodetic {
        let (sd, cd) = (distance_m / self.mean_radius()).sin_cos();
        let (sb, cb) = bearing_deg.to_radians().sin_cos();
        let (slat, clat) = from.lat_rad().sin_cos();
        let (slon, clon) = from.lon_rad().sin_cos();
        let (up, north, east) = (
            Vec3::new(clat * clon, clat * slon, slat),
            Vec3::new(-slat * clon, -slat * slon, clat),
            Vec3::new(-slon, clon, 0.0),
        );
        from_normal(up * cd + (north * cb + east * sb) * sd, from.height_m)
    }

    /// The halfway point along the great circle between two positions; height is the mean of the two.
    pub fn great_circle_midpoint(self, a: Geodetic, b: Geodetic) -> Geodetic {
        let (p1, p2) = (a.lat_rad(), b.lat_rad());
        let dl = b.lon_rad() - a.lon_rad();
        let bx = p1.cos() + p2.cos() * dl.cos();
        let by = p2.cos() * dl.sin();
        let lat = (p1.sin() + p2.sin()).atan2((bx * bx + by * by).sqrt());
        let lon = a.lon_rad() + by.atan2(bx);
        Geodetic::new(lat.to_degrees(), lon.to_degrees(), (a.height_m + b.height_m) / 2.0).normalized()
    }

    /// The unit normal to the plane of the great circle from `start` towards `end`, or `None` if
    /// they are the same or opposite points, which don't fix a circle. Computed as
    /// `n₁ × (n₂ − n₁)`, which equals `n₁ × n₂` but keeps full precision for short paths.
    fn path_normal(self, start: Geodetic, end: Geodetic) -> Option<Vec3> {
        let d = normal_difference(start, end);
        direction(self.surface_normal(start).cross(d), d.norm())
    }

    /// How far `p` is to the side of the great-circle path from `start` towards `end`, metres:
    /// positive to the right as you travel, negative to the left. The path is the whole great
    /// circle, so this is the distance to the nearest point on it. NaN if `start` and `end` are
    /// the same or opposite points.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// // A path due north along the prime meridian, and a point just east of it.
    /// let (start, end) = (Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(10.0, 0.0, 0.0));
    /// let off = WGS84.cross_track_distance_m(start, end, Geodetic::new(1.0, 0.1, 0.0));
    /// assert!((off - 11_117.8).abs() < 0.1);   // 11.1 km to the right
    /// ```
    #[doc(alias = "cross_track")]
    #[doc(alias = "xtd")]
    pub fn cross_track_distance_m(self, start: Geodetic, end: Geodetic, p: Geodetic) -> f64 {
        let Some(c) = self.path_normal(start, end) else {
            return f64::NAN;
        };
        // c · n_p, as c · (n_p − n_start): c is at right angles to n_start, so this is the same
        // value without the rounding noise of the larger vector.
        let side = c.dot(normal_difference(start, p));
        self.mean_radius() * (-side).atan2(c.cross(self.surface_normal(p)).norm())
    }

    /// How far along the great-circle path from `start` towards `end` the point nearest `p` is,
    /// metres; negative if it's behind `start`. NaN if `start` and `end` are the same or opposite
    /// points.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// let (start, end) = (Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(10.0, 0.0, 0.0));
    /// let along = WGS84.along_track_distance_m(start, end, Geodetic::new(1.0, 0.1, 0.0));
    /// assert!((along - 111_195.2).abs() < 0.1);   // level with 1° north
    /// ```
    #[doc(alias = "along_track")]
    pub fn along_track_distance_m(self, start: Geodetic, end: Geodetic, p: Geodetic) -> f64 {
        let Some(c) = self.path_normal(start, end) else {
            return f64::NAN;
        };
        // The angle from start to p dropped onto the path's plane, written with d = n_p − n_start
        // so nothing cancels: sin ∝ (n_start × d)·c and cos ∝ 1 + n_start·d.
        let (a, d) = (self.surface_normal(start), normal_difference(start, p));
        self.mean_radius() * a.cross(d).dot(c).atan2(1.0 + a.dot(d))
    }

    /// Where two great-circle paths cross, on the surface (height 0). Two great circles cross
    /// twice, on opposite sides of the Earth; this gives the crossing nearer the four endpoints.
    ///
    /// `None` if either path's ends are the same or opposite points, or both paths run along the
    /// same great circle.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// let g = |lat, lon| Geodetic::new(lat, lon, 0.0);
    /// let x = WGS84.great_circle_intersection(g(10.0, 20.0), g(30.0, 40.0), g(50.0, 60.0), g(70.0, 80.0)).unwrap();
    /// assert!((x.lat_deg - 40.318_643_07).abs() < 1e-8 && (x.lon_deg - 55.901_867_88).abs() < 1e-8);
    /// ```
    pub fn great_circle_intersection(self, a1: Geodetic, a2: Geodetic, b1: Geodetic, b2: Geodetic) -> Option<Geodetic> {
        let x = direction(self.path_normal(a1, a2)?.cross(self.path_normal(b1, b2)?), 1.0)?;
        let n = |g| self.surface_normal(g);
        let near = n(a1) + n(a2) + n(b1) + n(b2);
        Some(from_normal(if x.dot(near) < 0.0 { -x } else { x }, 0.0))
    }

    /// The centre of a set of positions: the direction of the sum of their surface normals
    /// (n-vectors), so it works across the poles and the antimeridian. Height is the average.
    ///
    /// `None` for no positions, or ones that balance out exactly (like two opposite points).
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// // Either side of the antimeridian: the centre is on it, not on the prime meridian.
    /// let m = WGS84.mean_position(&[Geodetic::new(0.0, 179.0, 0.0), Geodetic::new(0.0, -179.0, 0.0)]).unwrap();
    /// assert!((m.lon_deg.abs() - 180.0).abs() < 1e-12);
    /// ```
    #[doc(alias = "centroid")]
    #[doc(alias = "center")]
    pub fn mean_position(self, points: &[Geodetic]) -> Option<Geodetic> {
        let sum = points.iter().fold(Vec3::ZERO, |s, &p| s + self.surface_normal(p));
        let height = points.iter().map(|p| p.height_m).sum::<f64>() / points.len() as f64;
        Some(from_normal(direction(sum / points.len() as f64, 1.0)?, height))
    }

    /// The position a fraction `t` of the way along the great circle from `a` to `b`: `t = 0` is
    /// `a`, `1` is `b`, and it moves at a constant speed. Height changes in a straight line.
    /// Values of `t` outside 0 to 1 carry on along the same circle.
    ///
    /// Exactly opposite points have no single great circle between them, so they give NaN.
    ///
    /// ```
    /// use geo3d::{Geodetic, WGS84};
    /// // From 89° N on the prime meridian, over the pole, to 89° N on the antimeridian.
    /// let (a, b) = (Geodetic::new(89.0, 0.0, 0.0), Geodetic::new(89.0, 180.0, 0.0));
    /// let p = WGS84.great_circle_interpolate(a, b, 0.6);
    /// assert!((p.lat_deg - 89.8).abs() < 1e-9);   // 0.2° past the pole
    /// ```
    #[doc(alias = "slerp")]
    pub fn great_circle_interpolate(self, a: Geodetic, b: Geodetic, t: f64) -> Geodetic {
        let na = self.surface_normal(a);
        let height = a.height_m + (b.height_m - a.height_m) * t;
        // The unit vector at right angles to `a`, towards `b`, in the circle's plane:
        // nb − na(na·nb), written with d = nb − na so it keeps its precision for nearby points.
        let d = normal_difference(a, b);
        let Some(w) = direction(d - na * na.dot(d), d.norm()) else {
            // `b` is `a` (stay put), or opposite it (no single circle).
            return if d.norm() < 1.0 {
                Geodetic::new(a.lat_deg, a.lon_deg, height)
            } else {
                Geodetic::new(f64::NAN, f64::NAN, height)
            };
        };
        let angle = na.cross(d).norm().atan2(1.0 + na.dot(d)); // the angle from a to b
        let (s, c) = (t * angle).sin_cos();
        from_normal(na * c + w * s, height)
    }
}

#[cfg(test)]
mod tests {
    use crate::test_util::close;
    use crate::{Geodetic, WGS84};

    #[test]
    fn haversine_holds_up_at_and_near_the_antipode() {
        let half = core::f64::consts::PI * WGS84.mean_radius();
        // Exactly opposite points: rounding can push the haversine term past 1 (NaN with asin).
        for i in 0..=180 {
            let lat = -90.0 + i as f64;
            let d = WGS84.great_circle_distance_m(Geodetic::new(lat, 10.0, 0.0), Geodetic::new(-lat, -170.0, 0.0));
            assert!(close(d, half, 1e-6), "lat {lat}: {d}");
        }
        // 1 m short of half way round along the equator.
        let r = WGS84.mean_radius();
        let lon = 180.0 - (1.0 / r).to_degrees();
        let d = WGS84.great_circle_distance_m(Geodetic::new(0.0, 0.0, 0.0), Geodetic::new(0.0, lon, 0.0));
        assert!(close(d, half - 1.0, 1e-6), "{}", d - half);
        // And a millimetre apart.
        let mm = (0.001 / r).to_degrees();
        let d = WGS84.great_circle_distance_m(Geodetic::new(45.0, 7.0, 0.0), Geodetic::new(45.0 + mm, 7.0, 0.0));
        assert!(close(d, 0.001, 1e-9), "{d}");
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

    /// Short paths keep their precision. A path due north along a meridian has an exact answer:
    /// a point Δλ east of it at latitude φ is R·asin(cos φ·sin Δλ) to the right. (The textbook
    /// n₁ × n₂ form was 1.2 m out on a 1 mm path.)
    #[test]
    fn cross_and_along_track_stay_precise_on_short_paths() {
        let r = WGS84.mean_radius();
        for (lat, lon) in [(37.3, 123.4), (-62.1, -71.7), (5.5, 179.9)] {
            for len in [1e3, 1.0, 1e-3, 1e-6] {
                let dlat = (len / r).to_degrees();
                let (s, e) = (Geodetic::new(lat, lon, 0.0), Geodetic::new(lat + dlat, lon, 0.0));
                for off in [1.0, 100.0] {
                    let (plat, dlon) = (
                        lat + dlat / 2.0,
                        (off / r).to_degrees() / (lat + dlat / 2.0).to_radians().cos(),
                    );
                    let p = Geodetic::new(plat, lon + dlon, 0.0);
                    let want = r * (plat.to_radians().cos() * dlon.to_radians().sin()).asin();
                    let got = WGS84.cross_track_distance_m(s, e, p);
                    assert!((got - want).abs() < 1e-8, "{len} m path at {lat}: {got} vs {want}");
                    // The nearest point on the meridian is a little poleward of p:
                    // tan φ_foot = tan φ_p / cos Δλ.
                    let (sp, cp) = plat.to_radians().sin_cos();
                    let want = r * (sp.atan2(cp * dlon.to_radians().cos()) - lat.to_radians());
                    let along = WGS84.along_track_distance_m(s, e, p);
                    // Allow for how well the path's ends can be stored in degrees.
                    let tol = 1e-8 + off * 1e-9 / len;
                    assert!(
                        (along - want).abs() < tol,
                        "{len} m path, {off} m off: {along} vs {want}"
                    );
                }
            }
        }
    }

    /// Due north along a meridian, the latitude grows by exactly the angle travelled, right up to
    /// the pole. (The old asin form was up to 11 cm out arriving near a pole.)
    #[test]
    fn destination_is_exact_near_the_poles() {
        let r = WGS84.mean_radius();
        for i in 0..2_000 {
            let start = 80.0 + 10.0 * (i as f64 / 2_000.0);
            let end = 90.0 - 10f64.powi(-(i % 13));
            if end <= start {
                continue;
            }
            for (lat0, sign) in [(start, 1.0), (-start, -1.0)] {
                let bearing = if sign > 0.0 { 0.0 } else { 180.0 };
                let got = WGS84.great_circle_destination(
                    Geodetic::new(lat0, 33.0, 0.0),
                    bearing,
                    (end - start).to_radians() * r,
                );
                let err = (got.lat_deg - sign * end).abs().to_radians() * r;
                assert!(err < 1e-8, "from {lat0} to {}: {} ({err} m)", sign * end, got.lat_deg);
            }
        }
    }

    #[test]
    fn path_tools_signs_and_degenerate_cases() {
        let g = |lat, lon| Geodetic::new(lat, lon, 0.0);
        let r = WGS84.mean_radius();
        // Metres per degree on the sphere.
        let deg = r * core::f64::consts::PI / 180.0;
        // Due east along the equator: north of it is on the left.
        let (s, e) = (g(0.0, 10.0), g(0.0, 20.0));
        assert!(close(WGS84.cross_track_distance_m(s, e, g(1.0, 15.0)), -deg, 1e-6));
        assert!(close(WGS84.cross_track_distance_m(s, e, g(-1.0, 15.0)), deg, 1e-6));
        // Along-track: level with 5° along, and negative behind the start.
        assert!(close(WGS84.along_track_distance_m(s, e, g(1.0, 15.0)), 5.0 * deg, 1e-6));
        assert!(WGS84.along_track_distance_m(s, e, g(0.5, 9.0)) < 0.0);
        // On the path: zero across, and the distance along it.
        let on = WGS84.great_circle_destination(s, 90.0, 300_000.0);
        assert!(WGS84.cross_track_distance_m(s, e, on).abs() < 1e-6);
        assert!(close(WGS84.along_track_distance_m(s, e, on), 300_000.0, 1e-6));
        // A path over the pole, and one across the antimeridian.
        let x = WGS84.cross_track_distance_m(g(80.0, 0.0), g(80.0, 180.0), g(89.0, 90.0));
        assert!(close(x.abs(), deg, 1e-6), "{x}");
        let x = WGS84.cross_track_distance_m(g(0.0, 179.0), g(0.0, -179.0), g(1.0, 180.0));
        assert!(close(x, -deg, 1e-6), "{x}");
        // A path needs two distinct, non-opposite ends.
        assert!(WGS84.cross_track_distance_m(s, s, g(1.0, 1.0)).is_nan());
        assert!(WGS84.along_track_distance_m(s, g(0.0, -170.0), g(1.0, 1.0)).is_nan());
        // Crossings: the near one, and none for one circle or a pointless path.
        let c = WGS84
            .great_circle_intersection(g(0.0, -10.0), g(0.0, 10.0), g(-10.0, 5.0), g(10.0, 5.0))
            .unwrap();
        assert!(c.lat_deg.abs() < 1e-12 && close(c.lon_deg, 5.0, 1e-12), "{c:?}");
        assert!(WGS84
            .great_circle_intersection(g(0.0, 0.0), g(0.0, 10.0), g(0.0, 20.0), g(0.0, 30.0))
            .is_none());
        assert!(WGS84
            .great_circle_intersection(s, s, g(1.0, 1.0), g(2.0, 2.0))
            .is_none());
        // Centre: heights average; nothing, or two opposite points, has no centre.
        let m = WGS84
            .mean_position(&[Geodetic::new(10.0, 0.0, 100.0), Geodetic::new(-10.0, 0.0, 300.0)])
            .unwrap();
        assert!(
            m.lat_deg.abs() < 1e-12 && m.lon_deg == 0.0 && m.height_m == 200.0,
            "{m:?}"
        );
        assert!(WGS84.mean_position(&[]).is_none());
        assert!(WGS84.mean_position(&[g(0.0, 0.0), g(0.0, 180.0)]).is_none());
        let pole = WGS84
            .mean_position(&[g(80.0, 0.0), g(80.0, 120.0), g(80.0, -120.0)])
            .unwrap();
        assert!(close(pole.lat_deg, 90.0, 1e-12), "{pole:?}");
        // Interpolation: the ends, the middle, heights, and the antipodes.
        let (a, b) = (Geodetic::new(10.0, 20.0, 0.0), Geodetic::new(-30.0, 100.0, 1_000.0));
        let p0 = WGS84.great_circle_interpolate(a, b, 0.0);
        let p1 = WGS84.great_circle_interpolate(a, b, 1.0);
        assert!(close(p0.lat_deg, 10.0, 1e-12) && close(p0.lon_deg, 20.0, 1e-12) && p0.height_m == 0.0);
        assert!(close(p1.lat_deg, -30.0, 1e-12) && close(p1.lon_deg, 100.0, 1e-12) && p1.height_m == 1_000.0);
        let mid = WGS84.great_circle_interpolate(a, b, 0.5);
        let gm = WGS84.great_circle_midpoint(a, b);
        assert!(close(mid.lat_deg, gm.lat_deg, 1e-12) && close(mid.lon_deg, gm.lon_deg, 1e-12));
        let d = WGS84.great_circle_distance_m(a, b);
        assert!(close(
            WGS84.great_circle_distance_m(a, WGS84.great_circle_interpolate(a, b, 0.3)),
            0.3 * d,
            1e-6
        ));
        assert!(WGS84.great_circle_interpolate(a, a, 0.5).lat_deg == 10.0);
        assert!(WGS84
            .great_circle_interpolate(g(0.0, 0.0), g(0.0, 180.0), 0.5)
            .lat_deg
            .is_nan());
    }
}
