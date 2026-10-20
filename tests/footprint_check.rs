//! The prepared footprint against `geo`'s `coordinate_position`, the reference: every point placed
//! the same way, inside, on the boundary or outside.
#![cfg(feature = "footprint")]

use geo::algorithm::coordinate_position::CoordPos;
use geo::{Coord, CoordinatePosition, Geometry, LineString, MultiPolygon, Polygon};
use geo3d::{AltitudeBand, Containment, Footprint, Geodetic, Prism, Volume};

fn ring(pts: &[(f64, f64)]) -> LineString<f64> {
    LineString::from(pts.to_vec())
}

/// A wavy circle of `n` vertices around (cx, cy): a coastline-like ring.
fn wavy(n: usize, cx: f64, cy: f64, r: f64) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            let a = std::f64::consts::TAU * k as f64 / n as f64;
            let rr = r + 0.05 * r * (a * 40.0).sin();
            (cx + rr * a.cos(), cy + rr * a.sin())
        })
        .collect();
    v.push(v[0]);
    v
}

/// A square, a square with a hole, a 5,000-vertex ring around an 800-vertex hole, a spike, a flat
/// ring, an empty polygon, and three polygons that share edges and a corner.
fn shapes() -> Vec<Geometry<f64>> {
    let square = ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]);
    let hole = ring(&[(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0), (3.0, 3.0)]);
    let spike = ring(&[(0.0, 0.0), (10.0, 0.0), (5.0, 9.0), (5.0, 1.0), (4.0, 9.0), (0.0, 0.0)]);
    let flat = ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 0.0), (0.0, 0.0)]);
    vec![
        Polygon::new(square.clone(), vec![]).into(),
        Polygon::new(square, vec![hole]).into(),
        Polygon::new(ring(&wavy(5000, 5.0, 5.0, 4.0)), vec![ring(&wavy(800, 5.0, 5.0, 1.5))]).into(),
        Polygon::new(spike, vec![]).into(),
        Polygon::new(flat, vec![]).into(),
        Polygon::new(LineString::new(vec![]), vec![]).into(),
        MultiPolygon::new(vec![
            Polygon::new(
                ring(&[(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (0.0, 5.0), (0.0, 0.0)]),
                vec![],
            ),
            Polygon::new(
                ring(&[(5.0, 5.0), (10.0, 5.0), (10.0, 10.0), (5.0, 10.0), (5.0, 5.0)]),
                vec![],
            ),
            Polygon::new(
                ring(&[(5.0, 0.0), (10.0, 0.0), (10.0, 5.0), (5.0, 5.0), (5.0, 0.0)]),
                vec![],
            ),
        ])
        .into(),
    ]
}

fn footprint(g: &Geometry<f64>) -> Footprint {
    match g {
        Geometry::Polygon(p) => Footprint::from(p),
        Geometry::MultiPolygon(m) => Footprint::from(m),
        _ => unreachable!("a polygonal shape"),
    }
}

fn expected(pos: CoordPos) -> Containment {
    match pos {
        CoordPos::Inside => Containment::Inside,
        CoordPos::OnBoundary => Containment::OnBoundary,
        CoordPos::Outside => Containment::Outside,
    }
}

/// Every point placed as geo places it: on a grid that lands on edges and vertices, on every vertex
/// and edge midpoint, and at random.
#[test]
fn a_footprint_places_every_point_exactly_as_geo_does() {
    let mut seed = 7u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut checked = 0;
    for g in shapes() {
        let fp = footprint(&g);
        let mut pts: Vec<Coord<f64>> = Vec::new();
        for i in 0..=48 {
            for j in 0..=48 {
                pts.push(Coord {
                    x: f64::from(i) / 4.0 - 1.0,
                    y: f64::from(j) / 4.0 - 1.0,
                });
            }
        }
        let rings: Vec<&LineString<f64>> = match &g {
            Geometry::Polygon(p) => std::iter::once(p.exterior()).chain(p.interiors()).collect(),
            Geometry::MultiPolygon(m) => {
                m.0.iter()
                    .flat_map(|p| std::iter::once(p.exterior()).chain(p.interiors()))
                    .collect()
            }
            _ => unreachable!(),
        };
        for r in rings {
            for l in r.lines() {
                pts.push(l.start);
                pts.push(Coord {
                    x: (l.start.x + l.end.x) / 2.0,
                    y: (l.start.y + l.end.y) / 2.0,
                });
            }
        }
        for _ in 0..2000 {
            pts.push(Coord {
                x: next() * 12.0 - 1.0,
                y: next() * 12.0 - 1.0,
            });
        }
        for c in pts {
            let got = fp.containment(Geodetic::from_lon_lat(c.x, c.y, 0.0));
            assert_eq!(got, expected(g.coordinate_position(&c)), "{c:?} against {g:?}");
            checked += 1;
        }
    }
    assert!(checked > 40_000, "{checked} points checked");
}

/// Building from geo3d's own positions gives the same footprint as from `geo`, closed or not, and in
/// either winding direction.
#[test]
fn a_footprint_from_positions_closes_its_rings_and_takes_either_direction() {
    let at = |lon, lat| Geodetic::from_lon_lat(lon, lat, 123.0);
    let open = [at(0.0, 0.0), at(10.0, 0.0), at(10.0, 10.0), at(0.0, 10.0)];
    let mut closed = open.to_vec();
    closed.push(open[0]);
    let mut reversed = closed.clone();
    reversed.reverse();
    let reference = Footprint::from(&Polygon::new(
        ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]),
        vec![],
    ));
    for (name, fp) in [
        ("open", Footprint::polygon(&open, &[])),
        ("closed", Footprint::polygon(&closed, &[])),
        ("reversed", Footprint::polygon(&reversed, &[])),
    ] {
        for (lon, lat) in [(5.0, 5.0), (0.0, 5.0), (10.0, 10.0), (11.0, 5.0), (5.0, -0.1)] {
            let p = at(lon, lat);
            assert_eq!(fp.containment(p), reference.containment(p), "{name} at ({lon}, {lat})");
        }
    }
    assert_eq!(
        Footprint::polygon(&open, &[]).bounds(),
        Some((
            Geodetic::from_lon_lat(0.0, 0.0, 0.0),
            Geodetic::from_lon_lat(10.0, 10.0, 0.0)
        ))
    );
    assert_eq!(Footprint::polygon(&[], &[]).bounds(), None);
    assert_eq!(
        Footprint::polygon(&[], &[]).containment(at(0.0, 0.0)),
        Containment::Outside
    );
}

/// Longitude is x and latitude is y: an asymmetric zone tells a swap apart.
#[test]
fn longitude_is_x_and_latitude_is_y() {
    // a thin strip, wide in longitude and narrow in latitude
    let at = |lon, lat| Geodetic::from_lon_lat(lon, lat, 0.0);
    let strip = Footprint::polygon(
        &[at(150.0, -34.0), at(152.0, -34.0), at(152.0, -33.9), at(150.0, -33.9)],
        &[],
    );
    assert_eq!(strip.containment(at(151.0, -33.95)), Containment::Inside);
    assert_eq!(
        strip.containment(Geodetic::new(151.0, -33.95, 0.0)),
        Containment::Outside,
        "swapped"
    );
    assert_eq!(
        Geodetic::from_lon_lat(151.0, -33.95, 7.0),
        Geodetic::new(-33.95, 151.0, 7.0)
    );
}

/// A prism is inside the band and inside the footprint or on its boundary.
#[test]
fn a_prism_is_its_footprint_through_its_band() {
    let at = |lon, lat, h| Geodetic::from_lon_lat(lon, lat, h);
    let square = [
        at(0.0, 0.0, 0.0),
        at(10.0, 0.0, 0.0),
        at(10.0, 10.0, 0.0),
        at(0.0, 10.0, 0.0),
    ];
    let prism = Prism::new(
        Footprint::polygon(&square, &[]),
        AltitudeBand::new(Some(100.0), Some(200.0)),
    );
    assert!(prism.contains(at(5.0, 5.0, 150.0)));
    assert!(
        prism.contains(at(10.0, 5.0, 100.0)),
        "the boundary counts, as for every volume"
    );
    assert!(!prism.contains(at(5.0, 5.0, 250.0)), "above");
    assert!(!prism.contains(at(5.0, 5.0, 50.0)), "below");
    assert!(!prism.contains(at(11.0, 5.0, 150.0)), "beside");
}

/// A point that is not a number, or is infinite, is outside, and never panics: input from a
/// stream can hold anything.
#[test]
fn a_point_that_is_not_finite_is_outside() {
    let at = |lon, lat| Geodetic::from_lon_lat(lon, lat, 0.0);
    let square = Footprint::polygon(&[at(0.0, 0.0), at(10.0, 0.0), at(10.0, 10.0), at(0.0, 10.0)], &[]);
    for (lon, lat) in [
        (f64::NAN, 5.0),
        (5.0, f64::NAN),
        (f64::INFINITY, 5.0),
        (5.0, f64::NEG_INFINITY),
        (f64::NAN, f64::NAN),
    ] {
        assert_eq!(square.containment(at(lon, lat)), Containment::Outside, "({lon}, {lat})");
    }
    // a footprint whose vertices are not finite answers without panicking
    let broken = Footprint::polygon(&[at(0.0, 0.0), at(f64::NAN, 0.0), at(10.0, 10.0)], &[]);
    let _ = broken.containment(at(5.0, 5.0));
}
