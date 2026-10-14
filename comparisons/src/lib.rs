//! One adapter module per library, all with the same shape: degrees and metres in, degrees and
//! metres out. Each adapter is what a user of that library would write, unit and type conversions
//! included, so the benchmarks time them and the accuracy run tests exactly the same code.
//!
//! Every adapter returns `Option`: `None` means the library itself refused the input (a checked
//! constructor said no). geo3d never refuses.

/// Latitude (degrees), longitude (degrees), height above the ellipsoid (metres).
pub type Lla = [f64; 3];
/// ECEF x, y, z (metres).
pub type Xyz = [f64; 3];
/// East, north, up (metres).
pub type EnuV = [f64; 3];
/// Azimuth (degrees, `[0, 360)`), elevation (degrees), range (metres).
pub type AerV = [f64; 3];
/// Latitude, longitude (degrees).
pub type LatLon = [f64; 2];

/// A tiny deterministic generator (xorshift64*), so every run tests the same points.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed.max(1))
    }

    /// Uniform in `[0, 1)`.
    pub fn f64(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[lo, hi)`.
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.f64()
    }

    /// A point uniformly distributed over the sphere.
    pub fn lat_lon(&mut self) -> LatLon {
        [self.range(-1.0, 1.0).asin().to_degrees(), self.range(-180.0, 180.0)]
    }
}

/// A star-shaped (so never self-crossing) polygon of `n` vertices around `centre`, about
/// `size_deg` across, counter-clockwise, with longitudes in `[-180, 180)` as GeoJSON has them.
pub fn star(rng: &mut Rng, centre: LatLon, size_deg: f64, n: usize) -> Vec<LatLon> {
    // Bearings spread round the circle with gaps under 180°, so the centre is inside; taken in
    // decreasing order, which is counter-clockwise.
    let start = rng.range(0.0, 360.0);
    let bearings: Vec<f64> = (0..n)
        .rev()
        .map(|i| start + 360.0 * (i as f64 + rng.range(0.25, 0.75)) / n as f64)
        .collect();
    bearings
        .iter()
        .map(|&b| {
            let r = size_deg / 2.0 * rng.range(0.3, 1.0);
            let lat = centre[0] + r * b.to_radians().cos();
            let lon = centre[1] + r * b.to_radians().sin() / centre[0].to_radians().cos();
            [lat, (lon + 180.0).rem_euclid(360.0) - 180.0]
        })
        .collect()
}

pub mod with_geo3d {
    use super::*;
    use geo3d::{Aer, Ecef, Geodetic, LocalFrame, WGS84};

    #[inline]
    fn geo(g: Lla) -> Geodetic {
        Geodetic::new(g[0], g[1], g[2])
    }

    #[inline]
    pub fn to_ecef(g: Lla) -> Option<Xyz> {
        Some(WGS84.to_ecef(geo(g)).to_array())
    }

    #[inline]
    pub fn to_geodetic(e: Xyz) -> Option<Lla> {
        let g = WGS84.to_geodetic(Ecef::new(e[0], e[1], e[2]));
        Some([g.lat_deg, g.lon_deg, g.height_m])
    }

    pub struct Frame(LocalFrame);

    impl Frame {
        #[inline]
        pub fn new(obs: Lla) -> Option<Self> {
            Some(Frame(LocalFrame::new(geo(obs))))
        }

        #[inline]
        pub fn enu(&self, t: Lla) -> Option<EnuV> {
            Some(self.0.enu_of(geo(t)).vec().to_array())
        }

        #[inline]
        pub fn aer(&self, t: Lla) -> Option<AerV> {
            let Aer {
                azimuth_deg,
                elevation_deg,
                range_m,
            } = self.0.look_angles(geo(t));
            Some([azimuth_deg, elevation_deg, range_m])
        }
    }

    /// Distance, initial and final azimuth (degrees). Vincenty for now: `None` near antipodes.
    #[inline]
    pub fn inverse(a: LatLon, b: LatLon) -> Option<[f64; 3]> {
        let g = WGS84.geodesic_inverse(Geodetic::new(a[0], a[1], 0.0), Geodetic::new(b[0], b[1], 0.0))?;
        Some([g.distance_m, g.initial_bearing_deg, g.final_bearing_deg])
    }

    /// Arrival latitude, longitude and final azimuth.
    #[inline]
    pub fn direct(a: LatLon, azimuth_deg: f64, distance_m: f64) -> Option<[f64; 3]> {
        let (p, azi2) = WGS84.geodesic_direct(Geodetic::new(a[0], a[1], 0.0), azimuth_deg, distance_m);
        Some([p.lat_deg, p.lon_deg, azi2])
    }

    #[inline]
    pub fn great_circle_m(a: LatLon, b: LatLon) -> Option<f64> {
        Some(WGS84.great_circle_distance_m(Geodetic::new(a[0], a[1], 0.0), Geodetic::new(b[0], b[1], 0.0)))
    }

    /// Signed area in m², counter-clockwise positive.
    #[inline]
    pub fn polygon_area_m2(vertices: &[LatLon]) -> Option<f64> {
        let v: Vec<Geodetic> = vertices.iter().map(|p| Geodetic::new(p[0], p[1], 0.0)).collect();
        WGS84.polygon_area(&v).map(|p| p.area_m2)
    }

    /// Distance from `p` to the great circle through `a` and `b`, metres, on a sphere of radius
    /// 6 371 008.771 m (unsigned, for comparison).
    #[inline]
    pub fn cross_track_m(a: LatLon, b: LatLon, p: LatLon) -> Option<f64> {
        let g = |x: LatLon| Geodetic::new(x[0], x[1], 0.0);
        Some(WGS84.cross_track_distance_m(g(a), g(b), g(p)).abs())
    }
}

pub mod with_map_3d {
    use super::*;
    use map_3d::Ellipsoid::WGS84;

    #[inline]
    pub fn to_ecef(g: Lla) -> Option<Xyz> {
        let (x, y, z) = map_3d::geodetic2ecef(g[0].to_radians(), g[1].to_radians(), g[2], WGS84);
        Some([x, y, z])
    }

    #[inline]
    pub fn to_geodetic(e: Xyz) -> Option<Lla> {
        let (lat, lon, h) = map_3d::ecef2geodetic(e[0], e[1], e[2], WGS84);
        Some([lat.to_degrees(), lon.to_degrees(), h])
    }

    /// map_3d has no frame object: the observer is passed on every call. Only its conversion to
    /// radians can be done once.
    pub struct Frame([f64; 3]);

    impl Frame {
        #[inline]
        pub fn new(obs: Lla) -> Option<Self> {
            Some(Frame([obs[0].to_radians(), obs[1].to_radians(), obs[2]]))
        }

        #[inline]
        pub fn enu(&self, t: Lla) -> Option<EnuV> {
            let [lat0, lon0, h0] = self.0;
            let (e, n, u) = map_3d::geodetic2enu(t[0].to_radians(), t[1].to_radians(), t[2], lat0, lon0, h0, WGS84);
            Some([e, n, u])
        }

        #[inline]
        pub fn aer(&self, t: Lla) -> Option<AerV> {
            let [lat0, lon0, h0] = self.0;
            let (az, el, r) = map_3d::geodetic2aer(t[0].to_radians(), t[1].to_radians(), t[2], lat0, lon0, h0, WGS84);
            Some([az.to_degrees(), el.to_degrees(), r])
        }
    }

    /// Haversine on a sphere of radius 6 371 000 m; takes degrees.
    #[inline]
    pub fn great_circle_m(a: LatLon, b: LatLon) -> Option<f64> {
        Some(map_3d::distance((a[0], a[1]), (b[0], b[1])))
    }
}

pub mod with_nav_types {
    use super::*;
    use nav_types::{ECEF, WGS84};

    /// The checked constructor (`from_degrees_and_meters` panics instead on out-of-range input).
    #[inline]
    fn wgs(g: Lla) -> Option<WGS84<f64>> {
        WGS84::try_from_degrees_and_meters(g[0], g[1], g[2])
    }

    #[inline]
    pub fn to_ecef(g: Lla) -> Option<Xyz> {
        let e = ECEF::from(wgs(g)?);
        Some([e.x(), e.y(), e.z()])
    }

    #[inline]
    pub fn to_geodetic(e: Xyz) -> Option<Lla> {
        let g = WGS84::from(ECEF::new(e[0], e[1], e[2]));
        Some([g.latitude_degrees(), g.longitude_degrees(), g.altitude()])
    }

    /// nav-types has no frame object; the observer's ECEF position is the part that can be kept.
    pub struct Frame(ECEF<f64>);

    impl Frame {
        #[inline]
        pub fn new(obs: Lla) -> Option<Self> {
            Some(Frame(ECEF::from(wgs(obs)?)))
        }

        #[inline]
        pub fn enu(&self, t: Lla) -> Option<EnuV> {
            let v = ECEF::from(wgs(t)?) - self.0;
            Some([v.east(), v.north(), v.up()])
        }
    }

    /// Haversine on a sphere of radius a = 6 378 137 m, plus the height difference.
    #[inline]
    pub fn great_circle_m(a: LatLon, b: LatLon) -> Option<f64> {
        Some(wgs([a[0], a[1], 0.0])?.distance(&wgs([b[0], b[1], 0.0])?))
    }
}

pub mod with_sguaba {
    use super::*;
    use sguaba::math::RigidBodyTransform;
    use sguaba::systems::{Ecef, Wgs84};
    use sguaba::{system, Coordinate};
    use uom::si::angle::degree;
    use uom::si::f64::{Angle, Length};
    use uom::si::length::meter;

    system!(pub struct ObserverEnu using ENU);

    /// The checked builder: `None` for a latitude outside ±90° (or NaN).
    #[inline]
    fn wgs(g: Lla) -> Option<Wgs84> {
        Some(
            Wgs84::builder()
                .latitude(Angle::new::<degree>(g[0]))?
                .longitude(Angle::new::<degree>(g[1]))
                .altitude(Length::new::<meter>(g[2]))
                .build(),
        )
    }

    #[inline]
    fn metres(c: [Length; 3]) -> [f64; 3] {
        c.map(|l| l.get::<meter>())
    }

    #[inline]
    pub fn to_ecef(g: Lla) -> Option<Xyz> {
        Some(metres(Coordinate::<Ecef>::from_wgs84(&wgs(g)?).to_cartesian()))
    }

    #[inline]
    pub fn to_geodetic(e: Xyz) -> Option<Lla> {
        let c = sguaba::coordinate!(x = Length::new::<meter>(e[0]), y = Length::new::<meter>(e[1]), z = Length::new::<meter>(e[2]); in Ecef);
        let g = c.to_wgs84();
        Some([
            g.latitude().get::<degree>(),
            g.longitude().get::<degree>(),
            g.altitude().get::<meter>(),
        ])
    }

    pub struct Frame(RigidBodyTransform<Ecef, ObserverEnu>);

    impl Frame {
        #[inline]
        pub fn new(obs: Lla) -> Option<Self> {
            // SAFETY: `ObserverEnu` is only ever used for the frame at `obs` (sguaba uses `unsafe`
            // to mark frame-of-reference assumptions, not memory safety).
            Some(Frame(unsafe { RigidBodyTransform::ecef_to_enu_at(&wgs(obs)?) }))
        }

        #[inline]
        fn local(&self, t: Lla) -> Option<Coordinate<ObserverEnu>> {
            Some(self.0.transform(Coordinate::<Ecef>::from_wgs84(&wgs(t)?)))
        }

        #[inline]
        pub fn enu(&self, t: Lla) -> Option<EnuV> {
            Some(metres(self.local(t)?.to_cartesian()))
        }

        #[inline]
        pub fn aer(&self, t: Lla) -> Option<AerV> {
            let c = self.local(t)?;
            let b = c.bearing_from_origin()?;
            Some([
                b.azimuth().get::<degree>().rem_euclid(360.0),
                b.elevation().get::<degree>(),
                c.distance_from_origin().get::<meter>(),
            ])
        }
    }

    /// Haversine on a sphere of radius a = 6 378 137 m.
    #[inline]
    pub fn great_circle_m(a: LatLon, b: LatLon) -> Option<f64> {
        Some(
            wgs([a[0], a[1], 0.0])?
                .haversine_distance_on_surface(&wgs([b[0], b[1], 0.0])?)
                .get::<meter>(),
        )
    }
}

pub mod with_jord {
    use super::*;
    use jord::ellipsoidal::Ellipsoid;
    use jord::local::EnuFrame;
    use jord::spherical::{GreatCircle, Loop, Sphere};
    use jord::{GeocentricPosition, GeodeticPosition, LatLong, Length, NVector, PositionVector, Surface};

    #[inline]
    fn geo(g: Lla) -> GeodeticPosition {
        GeodeticPosition::new(NVector::from_lat_long_degrees(g[0], g[1]), Length::from_metres(g[2]))
    }

    #[inline]
    fn ecef(g: Lla) -> GeocentricPosition {
        Ellipsoid::WGS84.geodetic_to_geocentric_position(geo(g))
    }

    #[inline]
    pub fn to_ecef(g: Lla) -> Option<Xyz> {
        let v = ecef(g).as_metres();
        Some([v.x(), v.y(), v.z()])
    }

    #[inline]
    pub fn to_geodetic(e: Xyz) -> Option<Lla> {
        let g = Ellipsoid::WGS84.geocentric_to_geodetic_position(GeocentricPosition::from_metres(e[0], e[1], e[2]));
        let ll = LatLong::from_nvector(g.horizontal_position());
        Some([
            ll.latitude().as_degrees(),
            ll.longitude().as_degrees(),
            g.height().as_metres(),
        ])
    }

    pub struct Frame(EnuFrame);

    impl Frame {
        #[inline]
        pub fn new(obs: Lla) -> Option<Self> {
            Some(Frame(EnuFrame::from_geodetic(geo(obs), Ellipsoid::WGS84)))
        }

        #[inline]
        pub fn enu(&self, t: Lla) -> Option<EnuV> {
            let v = self.0.local_vector_to(ecef(t)).as_metres();
            Some([v.x(), v.y(), v.z()])
        }

        #[inline]
        pub fn aer(&self, t: Lla) -> Option<AerV> {
            let v = self.0.local_vector_to(ecef(t));
            Some([
                v.bearing().as_degrees(),
                v.elevation().as_degrees(),
                v.slant_range().as_metres(),
            ])
        }
    }

    /// Great-circle distance from n-vectors (an atan2 form, not haversine), radius 6 371 000.8 m.
    #[inline]
    pub fn great_circle_m(a: LatLon, b: LatLon) -> Option<f64> {
        let (p, q) = (
            NVector::from_lat_long_degrees(a[0], a[1]),
            NVector::from_lat_long_degrees(b[0], b[1]),
        );
        Some(Sphere::EARTH.distance(p, q).as_metres())
    }

    /// Distance from `p` to the great circle through `a` and `b`, metres, on the sphere of radius
    /// 6 371 000.8 m (unsigned, for comparison).
    #[inline]
    pub fn cross_track_m(a: LatLon, b: LatLon, p: LatLon) -> Option<f64> {
        let n = |x: LatLon| NVector::from_lat_long_degrees(x[0], x[1]);
        Some(
            Sphere::EARTH
                .cross_track_distance(n(p), GreatCircle::new(n(a), n(b)))
                .as_metres()
                .abs(),
        )
    }

    /// Area in m² on the sphere of radius 6 371 000.8 m, from the spherical excess (unsigned).
    #[inline]
    pub fn polygon_area_m2(vertices: &[LatLon]) -> Option<f64> {
        let vs: Vec<NVector> = vertices
            .iter()
            .map(|p| NVector::from_lat_long_degrees(p[0], p[1]))
            .collect();
        let r = Sphere::EARTH.radius().as_metres();
        Some(Loop::new(&vs).spherical_excess().as_radians() * r * r)
    }
}

pub mod with_geographiclib {
    use super::*;
    use geographiclib_rs::{DirectGeodesic, Geodesic, InverseGeodesic, PolygonArea, Winding};

    /// Distance, initial and final azimuth (degrees, wrapped into `[0, 360)` like geo3d's).
    #[inline]
    pub fn inverse(a: LatLon, b: LatLon) -> Option<[f64; 3]> {
        let (s12, azi1, azi2, _a12): (f64, f64, f64, f64) = Geodesic::wgs84().inverse(a[0], a[1], b[0], b[1]);
        Some([s12, azi1.rem_euclid(360.0), azi2.rem_euclid(360.0)])
    }

    /// Arrival latitude, longitude and final azimuth.
    #[inline]
    pub fn direct(a: LatLon, azimuth_deg: f64, distance_m: f64) -> Option<[f64; 3]> {
        let (lat, lon, azi2): (f64, f64, f64) = Geodesic::wgs84().direct(a[0], a[1], azimuth_deg, distance_m);
        Some([lat, lon, azi2])
    }

    /// Signed area in m², counter-clockwise positive.
    #[inline]
    pub fn polygon_area_m2(vertices: &[LatLon]) -> Option<f64> {
        let g = Geodesic::wgs84();
        let mut pa = PolygonArea::new(&g, Winding::CounterClockwise);
        for p in vertices {
            pa.add_point(p[0], p[1]);
        }
        Some(pa.compute(true).1)
    }
}

/// geo's `Geodesic` metric space calls geographiclib-rs (distance only, or position only).
pub mod with_geo {
    use super::*;
    use geo::{
        ChamberlainDuquetteArea, Coord, CrossTrackDistance, Destination, Distance, Geodesic, GeodesicArea, Haversine,
        LineString, Point, Polygon,
    };

    #[inline]
    fn point(p: LatLon) -> Point {
        Point::new(p[1], p[0]) // x = longitude, y = latitude
    }

    #[inline]
    pub fn geodesic_distance_m(a: LatLon, b: LatLon) -> Option<f64> {
        Some(Geodesic.distance(point(a), point(b)))
    }

    #[inline]
    pub fn geodesic_destination(a: LatLon, azimuth_deg: f64, distance_m: f64) -> Option<LatLon> {
        let p = Geodesic.destination(point(a), azimuth_deg, distance_m);
        Some([p.y(), p.x()])
    }

    /// Haversine on a sphere of radius 6 371 008.8 m.
    #[inline]
    pub fn great_circle_m(a: LatLon, b: LatLon) -> Option<f64> {
        Some(Haversine.distance(point(a), point(b)))
    }

    #[inline]
    fn polygon(vertices: &[LatLon]) -> Polygon {
        let ring: Vec<Coord> = vertices.iter().map(|p| Coord { x: p[1], y: p[0] }).collect();
        Polygon::new(LineString::from(ring), vec![])
    }

    /// `CrossTrackDistance` (haversine bearings), radius 6 371 008.8 m; unsigned.
    #[inline]
    pub fn cross_track_m(a: LatLon, b: LatLon, p: LatLon) -> Option<f64> {
        Some(point(p).cross_track_distance(&point(a), &point(b)))
    }

    /// `GeodesicArea`, which calls geographiclib-rs: signed m², counter-clockwise positive.
    #[inline]
    pub fn geodesic_area_m2(vertices: &[LatLon]) -> Option<f64> {
        Some(polygon(vertices).geodesic_area_signed())
    }

    /// `ChamberlainDuquetteArea`, on a sphere of radius 6 378 137 m: signed m².
    #[inline]
    pub fn spherical_area_m2(vertices: &[LatLon]) -> Option<f64> {
        Some(polygon(vertices).chamberlain_duquette_signed_area())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every adapter answers the same simple question the same way, so a unit or argument-order
    /// slip in one of them shows up here rather than as a fake accuracy result.
    /// A 1° square: the geodesic adapters agree exactly; the spherical ones within 1 %.
    #[test]
    fn area_adapters_agree_on_a_square() {
        let square = [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]];
        let a = with_geo3d::polygon_area_m2(&square).unwrap();
        assert!((a - 12_308_778_361.469).abs() < 0.01, "{a}");
        for (name, x) in [
            ("geographiclib-rs", with_geographiclib::polygon_area_m2(&square)),
            ("geo geodesic", with_geo::geodesic_area_m2(&square)),
        ] {
            assert!((x.unwrap() - a).abs() < 0.01, "{name}: {x:?}");
        }
        for (name, x) in [
            ("geo spherical", with_geo::spherical_area_m2(&square)),
            ("jord", with_jord::polygon_area_m2(&square)),
        ] {
            assert!((x.unwrap() / a - 1.0).abs() < 0.01, "{name}: {x:?}");
        }
    }

    #[test]
    fn adapters_agree_on_a_simple_case() {
        let obs = [-33.94, 151.18, 21.0];
        let t = [-33.80, 151.30, 2_500.0];
        let e = with_geo3d::to_ecef(t).unwrap();
        for (name, x) in [
            ("map_3d", with_map_3d::to_ecef(t)),
            ("nav-types", with_nav_types::to_ecef(t)),
            ("sguaba", with_sguaba::to_ecef(t)),
            ("jord", with_jord::to_ecef(t)),
        ] {
            let x = x.unwrap();
            assert!(
                (0..3).all(|i| (x[i] - e[i]).abs() < 1e-6),
                "{name} to_ecef {x:?} vs {e:?}"
            );
        }
        let enu = with_geo3d::Frame::new(obs).unwrap().enu(t).unwrap();
        for (name, v) in [
            ("map_3d", with_map_3d::Frame::new(obs).unwrap().enu(t)),
            ("nav-types", with_nav_types::Frame::new(obs).unwrap().enu(t)),
            ("sguaba", with_sguaba::Frame::new(obs).unwrap().enu(t)),
            ("jord", with_jord::Frame::new(obs).unwrap().enu(t)),
        ] {
            let v = v.unwrap();
            assert!(
                (0..3).all(|i| (v[i] - enu[i]).abs() < 1e-6),
                "{name} enu {v:?} vs {enu:?}"
            );
        }
        let aer = with_geo3d::Frame::new(obs).unwrap().aer(t).unwrap();
        for (name, v) in [
            ("map_3d", with_map_3d::Frame::new(obs).unwrap().aer(t)),
            ("sguaba", with_sguaba::Frame::new(obs).unwrap().aer(t)),
            ("jord", with_jord::Frame::new(obs).unwrap().aer(t)),
        ] {
            let v = v.unwrap();
            assert!(
                (0..3).all(|i| (v[i] - aer[i]).abs() < 1e-6),
                "{name} aer {v:?} vs {aer:?}"
            );
        }
        let (a, b) = ([33.9425, -118.4081], [40.6413, -73.7781]);
        let g = with_geographiclib::inverse(a, b).unwrap();
        let ours = with_geo3d::inverse(a, b).unwrap();
        assert!((g[0] - ours[0]).abs() < 1e-3 && (g[1] - ours[1]).abs() < 1e-8 && (g[2] - ours[2]).abs() < 1e-8);
        assert!((with_geo::geodesic_distance_m(a, b).unwrap() - g[0]).abs() < 1e-6);
        let d = with_geographiclib::direct(a, g[1], g[0]).unwrap();
        let p = with_geo::geodesic_destination(a, g[1], g[0]).unwrap();
        assert!((d[0] - b[0]).abs() < 1e-9 && (d[1] - b[1]).abs() < 1e-9 && (p[0] - b[0]).abs() < 1e-9);
        for (name, s) in [
            ("geo3d", with_geo3d::great_circle_m(a, b)),
            ("map_3d", with_map_3d::great_circle_m(a, b)),
            ("nav-types", with_nav_types::great_circle_m(a, b)),
            ("sguaba", with_sguaba::great_circle_m(a, b)),
            ("jord", with_jord::great_circle_m(a, b)),
            ("geo", with_geo::great_circle_m(a, b)),
        ] {
            let s = s.unwrap();
            assert!((s - g[0]).abs() < 0.01 * g[0], "{name} great circle {s} vs {}", g[0]);
        }
    }
}
