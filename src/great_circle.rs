//! Quick spherical estimates: great-circle distance, bearing, destination and midpoint on a sphere
//! of the ellipsoid's mean radius. Within 0.6 % of the true surface distance on WGS84.

use crate::ellipsoid::Ellipsoid;
use crate::position::Geodetic;

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
        y.atan2(x).to_degrees().rem_euclid(360.0)
    }

    /// The position `distance_m` along the great circle from `from` on an initial `bearing_deg`;
    /// height is carried over.
    pub fn great_circle_destination(self, from: Geodetic, bearing_deg: f64, distance_m: f64) -> Geodetic {
        let delta = distance_m / self.mean_radius();
        let (sd, cd) = delta.sin_cos();
        let (s1, c1) = from.lat_rad().sin_cos();
        let (sb, cb) = bearing_deg.to_radians().sin_cos();
        let s2 = (s1 * cd + c1 * sd * cb).clamp(-1.0, 1.0);
        let lat2 = s2.asin();
        let lon2 = from.lon_rad() + (sb * sd * c1).atan2(cd - s1 * s2);
        Geodetic::new(lat2.to_degrees(), lon2.to_degrees(), from.height_m).normalized()
    }

    /// The halfway point along the great circle between two positions; height is the mean of the two.
    pub fn great_circle_midpoint(self, a: Geodetic, b: Geodetic) -> Geodetic {
        let (p1, p2) = (a.lat_rad(), b.lat_rad());
        let dl = b.lon_rad() - a.lon_rad();
        let bx = p2.cos() * dl.cos();
        let by = p2.cos() * dl.sin();
        let lat = (p1.sin() + p2.sin()).atan2(((p1.cos() + bx).powi(2) + by * by).sqrt());
        let lon = a.lon_rad() + by.atan2(p1.cos() + bx);
        Geodetic::new(lat.to_degrees(), lon.to_degrees(), (a.height_m + b.height_m) / 2.0).normalized()
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
}
