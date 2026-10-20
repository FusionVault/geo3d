//! A ground footprint, prepared once so that "is this point inside?" is fast: one polygon or
//! several, each with holes, as an airspace sector, a port, a delivery zone or a country.
//!
//! Each ring's edges are indexed by latitude, so a point meets only the few edges that cross its
//! latitude instead of every edge. A coastline or a border has thousands. A point outside the
//! footprint's bounding box meets none.
//!
//! **Edges are straight in longitude and latitude** (planar, in degrees), as in GeoJSON and most
//! 2D GIS. Rings that cross the antimeridian are not supported: split them at ±180° first.
//!
//! **The answer is exact.** A ring's winding number counts only the edges whose latitude range
//! holds the point's, each with Shewchuk's adaptive orientation test (the `robust` crate). So a
//! point on an edge or a vertex is always on the boundary, whatever the rounding, and neither
//! answer depends on the order of the edges. It is the answer `geo`'s `coordinate_position` gives.

use alloc::vec;
use alloc::vec::Vec;

use crate::position::Geodetic;

/// Where a point lies against a [`Footprint`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Containment {
    /// In the interior.
    Inside,
    /// On an edge or a vertex.
    OnBoundary,
    /// Outside, or in a hole.
    Outside,
}

/// Buckets never outnumber this, so a large footprint's index stays small.
const MAX_BUCKETS: usize = 1024;
/// Rebuild with fewer buckets when edges spanning many buckets would be stored this many times over.
const MAX_COPIES: usize = 8;

/// A vertex: `x` is the longitude and `y` the latitude, degrees.
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
struct Xy {
    x: f64,
    y: f64,
}

impl From<Xy> for robust::Coord<f64> {
    #[inline]
    fn from(p: Xy) -> Self {
        robust::Coord { x: p.x, y: p.y }
    }
}

/// One ring: its edges, and for each bucket of its latitude range the edges that overlap it.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
struct Ring {
    edges: Vec<(Xy, Xy)>,
    /// A ring of one point has no edge: a point equal to it is on its boundary.
    single: Option<Xy>,
    y0: f64,
    y1: f64,
    dy: f64,
    buckets: Vec<Vec<u32>>,
}

impl Ring {
    /// The ring through `points`, closed whether or not the last repeats the first.
    fn new(points: impl IntoIterator<Item = Xy>) -> Ring {
        let points: Vec<Xy> = points.into_iter().collect();
        let mut edges: Vec<(Xy, Xy)> = points.windows(2).map(|w| (w[0], w[1])).collect();
        if let (Some(&first), Some(&last)) = (points.first(), points.last()) {
            if points.len() > 1 && first != last {
                edges.push((last, first));
            }
        }
        let single = (points.len() == 1).then(|| points[0]);
        let (mut y0, mut y1) = (f64::INFINITY, f64::NEG_INFINITY);
        for (s, e) in &edges {
            y0 = y0.min(s.y.min(e.y));
            y1 = y1.max(s.y.max(e.y));
        }
        let mut ring = Ring {
            edges,
            single,
            y0,
            y1,
            dy: 1.0,
            buckets: Vec::new(),
        };
        let mut n = ring.edges.len().clamp(1, MAX_BUCKETS);
        loop {
            ring.dy = if y1 > y0 { (y1 - y0) / n as f64 } else { 1.0 };
            ring.buckets = vec![Vec::new(); n];
            let mut copies = 0;
            for (i, (s, e)) in ring.edges.iter().enumerate() {
                let (lo, hi) = (ring.bucket(s.y.min(e.y)), ring.bucket(s.y.max(e.y)));
                copies += hi - lo + 1;
                for b in &mut ring.buckets[lo..=hi] {
                    b.push(i as u32);
                }
            }
            if n == 1 || copies <= MAX_COPIES * ring.edges.len() {
                return ring;
            }
            n /= 2;
        }
    }

    /// The bucket of latitude `y`. It grows with `y`, so an edge's buckets run from its lowest
    /// point's to its highest point's.
    #[inline]
    fn bucket(&self, y: f64) -> usize {
        (((y - self.y0) / self.dy) as usize).min(self.buckets.len() - 1)
    }

    /// Where `c` lies against this ring.
    fn containment(&self, c: Xy) -> Containment {
        if self.edges.is_empty() {
            return match self.single {
                Some(p) if p == c => Containment::OnBoundary,
                _ => Containment::Outside,
            };
        }
        // no edge's latitude range holds `c`'s: no edge takes part
        if !(c.y >= self.y0 && c.y <= self.y1) {
            return Containment::Outside;
        }
        let mut winding = 0i32;
        for &i in &self.buckets[self.bucket(c.y)] {
            let (s, e) = self.edges[i as usize];
            // an upward edge includes its start and excludes its end, a downward edge the reverse,
            // a horizontal edge counts only for the boundary, and a crossing is strictly to the right
            if s.y <= c.y {
                if e.y >= c.y {
                    let o = robust::orient2d(s.into(), e.into(), c.into());
                    if o > 0.0 && e.y != c.y {
                        winding += 1;
                    } else if o == 0.0 && between(c.x, s.x, e.x) {
                        return Containment::OnBoundary;
                    }
                }
            } else if e.y <= c.y {
                let o = robust::orient2d(s.into(), e.into(), c.into());
                if o < 0.0 {
                    winding -= 1;
                } else if o == 0.0 && between(c.x, s.x, e.x) {
                    return Containment::OnBoundary;
                }
            }
        }
        if winding == 0 {
            Containment::Outside
        } else {
            Containment::Inside
        }
    }
}

/// Whether `v` lies between `a` and `b`, ends included, in either order.
#[inline]
fn between(v: f64, a: f64, b: f64) -> bool {
    if a < b {
        v >= a && v <= b
    } else {
        v >= b && v <= a
    }
}

/// A ground footprint, prepared for fast point tests: one polygon or several, each with holes.
///
/// Build it once, then ask [`containment`](Footprint::containment) for as many points as you like.
/// Heights are ignored. Extrude it through an [`AltitudeBand`](crate::AltitudeBand) with
/// [`Prism`](crate::Prism) for an airspace volume.
///
/// ```
/// use geo3d::{Containment, Footprint, Geodetic};
/// // a 1° square near Sydney, with a hole in the middle
/// let at = |lon, lat| Geodetic::from_lon_lat(lon, lat, 0.0);
/// let square = [at(151.0, -34.0), at(152.0, -34.0), at(152.0, -33.0), at(151.0, -33.0)];
/// let hole = [at(151.4, -33.6), at(151.6, -33.6), at(151.6, -33.4), at(151.4, -33.4)];
/// let zone = Footprint::polygon(&square, &[&hole]);
/// assert_eq!(zone.containment(at(151.2, -33.8)), Containment::Inside);
/// assert_eq!(zone.containment(at(151.5, -33.5)), Containment::Outside); // in the hole
/// assert_eq!(zone.containment(at(152.0, -33.5)), Containment::OnBoundary);
/// ```
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Footprint {
    /// Each polygon's exterior ring and holes.
    polygons: Vec<(Ring, Vec<Ring>)>,
    /// The south-west and north-east corners of every vertex, or `None` when there is none.
    bounds: Option<(Xy, Xy)>,
}

impl Footprint {
    /// One polygon: its exterior ring and its holes.
    ///
    /// A ring is its vertices in order. It is closed whether or not the last vertex repeats the
    /// first, and either winding direction works. Heights are ignored.
    pub fn polygon(exterior: &[Geodetic], holes: &[&[Geodetic]]) -> Self {
        Footprint::multi_polygon(&[(exterior, holes)])
    }

    /// Several polygons, each its exterior ring and its holes, as in
    /// [`polygon`](Footprint::polygon).
    ///
    /// A point on the boundaries of an odd number of them is on the boundary, and on an even number
    /// it is inside: two polygons that share an edge cover it.
    pub fn multi_polygon(polygons: &[(&[Geodetic], &[&[Geodetic]])]) -> Self {
        let xy = |ring: &[Geodetic]| {
            ring.iter()
                .map(|g| Xy {
                    x: g.lon_deg,
                    y: g.lat_deg,
                })
                .collect::<Vec<_>>()
        };
        Footprint::from_rings(
            polygons
                .iter()
                .map(|(exterior, holes)| (xy(exterior), holes.iter().map(|h| xy(h)).collect())),
        )
    }

    /// The footprint of rings already in longitude, latitude order.
    fn from_rings(polygons: impl IntoIterator<Item = (Vec<Xy>, Vec<Vec<Xy>>)>) -> Self {
        let mut bounds: Option<(Xy, Xy)> = None;
        let polygons = polygons
            .into_iter()
            // an empty polygon takes no part
            .filter(|(exterior, _)| !exterior.is_empty())
            .map(|(exterior, holes)| {
                for p in &exterior {
                    bounds = Some(match bounds {
                        None => (*p, *p),
                        Some((lo, hi)) => (
                            Xy {
                                x: lo.x.min(p.x),
                                y: lo.y.min(p.y),
                            },
                            Xy {
                                x: hi.x.max(p.x),
                                y: hi.y.max(p.y),
                            },
                        ),
                    });
                }
                (Ring::new(exterior), holes.into_iter().map(Ring::new).collect())
            })
            .collect();
        Footprint { polygons, bounds }
    }

    /// Where `point` lies: inside, on the boundary, or outside (holes included). Its height is
    /// ignored.
    pub fn containment(&self, point: Geodetic) -> Containment {
        let c = Xy {
            x: point.lon_deg,
            y: point.lat_deg,
        };
        // the holes lie inside the exteriors, so the exteriors' box bounds everything
        match self.bounds {
            Some((lo, hi)) if c.x >= lo.x && c.x <= hi.x && c.y >= lo.y && c.y <= hi.y => {}
            _ => return Containment::Outside,
        }
        let (mut inside, mut boundaries) = (false, 0usize);
        for (exterior, holes) in &self.polygons {
            match exterior.containment(c) {
                Containment::Outside => {}
                Containment::OnBoundary => boundaries += 1,
                Containment::Inside => {
                    let mut in_polygon = true;
                    for hole in holes {
                        match hole.containment(c) {
                            Containment::Outside => {}
                            Containment::OnBoundary => {
                                boundaries += 1;
                                in_polygon = false;
                                break;
                            }
                            Containment::Inside => {
                                in_polygon = false;
                                break;
                            }
                        }
                    }
                    inside |= in_polygon;
                }
            }
        }
        if boundaries % 2 == 1 {
            Containment::OnBoundary
        } else if inside {
            Containment::Inside
        } else {
            Containment::Outside
        }
    }

    /// The south-west and north-east corners of the footprint's bounding box, or `None` for an
    /// empty footprint. Heights are 0.
    pub fn bounds(&self) -> Option<(Geodetic, Geodetic)> {
        self.bounds.map(|(lo, hi)| {
            (
                Geodetic::from_lon_lat(lo.x, lo.y, 0.0),
                Geodetic::from_lon_lat(hi.x, hi.y, 0.0),
            )
        })
    }
}

#[cfg(feature = "geo-types")]
mod geo_types_impls {
    use super::{Footprint, Xy};
    use alloc::vec::Vec;

    fn ring(ls: &geo_types::LineString<f64>) -> Vec<Xy> {
        ls.0.iter().map(|c| Xy { x: c.x, y: c.y }).collect()
    }

    fn polygon(p: &geo_types::Polygon<f64>) -> (Vec<Xy>, Vec<Vec<Xy>>) {
        (ring(p.exterior()), p.interiors().iter().map(ring).collect())
    }

    impl From<&geo_types::Polygon<f64>> for super::Footprint {
        fn from(p: &geo_types::Polygon<f64>) -> Self {
            Footprint::from_rings([polygon(p)])
        }
    }

    impl From<&geo_types::MultiPolygon<f64>> for super::Footprint {
        fn from(m: &geo_types::MultiPolygon<f64>) -> Self {
            Footprint::from_rings(m.0.iter().map(polygon))
        }
    }
}
