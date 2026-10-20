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

/// Adds what the edge from `$s` to `$e` contributes to `$c`'s winding number to `$winding`, or
/// returns [`Containment::OnBoundary`] from the enclosing function when `$c` is on the edge. An
/// upward edge includes its start and excludes its end, a downward edge the reverse, a horizontal
/// edge counts only for the boundary, and a crossing is strictly to the right. The one edge rule,
/// for an indexed ring and a ring used once alike; a macro so each loop keeps its plain branches.
macro_rules! edge_winding {
    ($s:expr, $e:expr, $c:expr, $winding:ident) => {{
        let (s, e, c): (Xy, Xy, Xy) = ($s, $e, $c);
        if s.y <= c.y {
            if e.y >= c.y {
                let o = robust::orient2d(s.into(), e.into(), c.into());
                if o > 0.0 && e.y != c.y {
                    $winding += 1;
                } else if o == 0.0 && between(c.x, s.x, e.x) {
                    return Containment::OnBoundary;
                }
            }
        } else if e.y <= c.y {
            let o = robust::orient2d(s.into(), e.into(), c.into());
            if o < 0.0 {
                $winding -= 1;
            } else if o == 0.0 && between(c.x, s.x, e.x) {
                return Containment::OnBoundary;
            }
        }
    }};
}

/// Tallies one polygon into `$inside` and `$boundaries`, its exterior placing the point at
/// `$exterior` and each `$h` of `$holes` at `$place`: inside a polygon is inside its exterior and in
/// none of its holes, asked only while it can still be. The one rule for combining rings, for a
/// footprint and a polygon used once alike; see [`tallied`] for the answer. A macro so each loop
/// keeps its plain branches.
macro_rules! tally_polygon {
    ($exterior:expr, $holes:expr, |$h:ident| $place:expr, $inside:ident, $boundaries:ident) => {
        match $exterior {
            Containment::Outside => {}
            Containment::OnBoundary => $boundaries += 1,
            Containment::Inside => {
                let mut in_polygon = true;
                for $h in $holes {
                    match $place {
                        Containment::Outside => {}
                        Containment::OnBoundary => {
                            $boundaries += 1;
                            in_polygon = false;
                            break;
                        }
                        Containment::Inside => {
                            in_polygon = false;
                            break;
                        }
                    }
                }
                $inside |= in_polygon;
            }
        }
    };
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
    #[inline]
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
            edge_winding!(s, e, c, winding);
        }
        by_winding(winding)
    }
}

#[inline]
fn by_winding(winding: i32) -> Containment {
    if winding == 0 {
        Containment::Outside
    } else {
        Containment::Inside
    }
}

/// Where `c` lies against the ring through `vertices`, edge by edge with no index. It is closed
/// as a [`Ring`] is: whether or not the last vertex repeats the first, and a ring of one vertex is
/// that point.
#[inline]
fn ring_once(vertices: impl IntoIterator<Item = Xy>, c: Xy) -> Containment {
    let mut vertices = vertices.into_iter();
    let Some(first) = vertices.next() else {
        return Containment::Outside;
    };
    let (mut last, mut winding, mut single) = (first, 0i32, true);
    for v in vertices {
        edge_winding!(last, v, c, winding);
        (last, single) = (v, false);
    }
    if single {
        return if first == c {
            Containment::OnBoundary
        } else {
            Containment::Outside
        };
    }
    if last != first {
        edge_winding!(last, first, c, winding);
    }
    by_winding(winding)
}

/// What the boundaries and insides tallied over a footprint's polygons make: a point on the
/// boundaries of an odd number of polygons is on the boundary.
#[inline]
fn tallied(inside: bool, boundaries: usize) -> Containment {
    if boundaries % 2 == 1 {
        Containment::OnBoundary
    } else if inside {
        Containment::Inside
    } else {
        Containment::Outside
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
    /// A point on the boundaries of an odd number of them is on the boundary. On an even number it
    /// is outside, unless it is inside another of them. This is `geo`'s rule. So a point on an edge
    /// that two of them share is outside: merge the two into one polygon when it should be inside.
    ///
    /// ```
    /// use geo3d::{Containment, Footprint, Geodetic};
    /// let at = |lon, lat| Geodetic::from_lon_lat(lon, lat, 0.0);
    /// let west = [at(0.0, 0.0), at(1.0, 0.0), at(1.0, 1.0), at(0.0, 1.0)];
    /// let east = [at(1.0, 0.0), at(2.0, 0.0), at(2.0, 1.0), at(1.0, 1.0)];
    /// let both = Footprint::multi_polygon(&[(&west, &[]), (&east, &[])]);
    /// assert_eq!(both.containment(at(0.5, 0.5)), Containment::Inside);
    /// assert_eq!(both.containment(at(1.0, 0.5)), Containment::Outside); // on the edge both share
    /// assert_eq!(both.containment(at(0.0, 0.5)), Containment::OnBoundary); // on one boundary
    /// ```
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
    #[inline]
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
            tally_polygon!(
                exterior.containment(c),
                holes,
                |hole| hole.containment(c),
                inside,
                boundaries
            );
        }
        tallied(inside, boundaries)
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

/// Where `point` lies against one polygon used once, tested edge by edge with no index: its rings
/// are its exterior and then its holes, each its vertices in order.
///
/// It gives a [`Footprint`]'s answer for the same rings, by the same exact rule, valid or not. It
/// allocates nothing, so a ring can be read straight from where it lies, such as a row of WKB. A
/// polygon tested against many points is faster as a [`Footprint`]: its edges are indexed once.
///
/// ```
/// use geo3d::{polygon_containment, Containment, Geodetic};
/// let at = |lon, lat| Geodetic::from_lon_lat(lon, lat, 0.0);
/// let square = [at(151.0, -34.0), at(152.0, -34.0), at(152.0, -33.0), at(151.0, -33.0)];
/// let hole = [at(151.4, -33.6), at(151.6, -33.6), at(151.6, -33.4), at(151.4, -33.4)];
/// assert_eq!(polygon_containment([&square[..], &hole], at(151.2, -33.8)), Containment::Inside);
/// assert_eq!(polygon_containment([&square[..], &hole], at(151.5, -33.5)), Containment::Outside);
/// ```
#[inline]
pub fn polygon_containment<R>(rings: impl IntoIterator<Item = R>, point: Geodetic) -> Containment
where
    R: IntoIterator,
    R::Item: core::borrow::Borrow<Geodetic>,
{
    multi_polygon_containment([rings], point)
}

/// Where `point` lies against several polygons used once, each its rings as in
/// [`polygon_containment`]. A point on the boundaries of an odd number of them is on the boundary,
/// as in a [`Footprint`].
#[inline]
pub fn multi_polygon_containment<P, R>(polygons: impl IntoIterator<Item = P>, point: Geodetic) -> Containment
where
    P: IntoIterator<Item = R>,
    R: IntoIterator,
    R::Item: core::borrow::Borrow<Geodetic>,
{
    let c = Xy {
        x: point.lon_deg,
        y: point.lat_deg,
    };
    // a footprint's bounding box puts a point that is not a number outside; so does this
    if !(c.x.is_finite() && c.y.is_finite()) {
        return Containment::Outside;
    }
    let place = |ring: R| {
        ring_once(
            ring.into_iter().map(|g| {
                let g = core::borrow::Borrow::<Geodetic>::borrow(&g);
                Xy {
                    x: g.lon_deg,
                    y: g.lat_deg,
                }
            }),
            c,
        )
    };
    let (mut inside, mut boundaries) = (false, 0usize);
    for rings in polygons {
        let mut rings = rings.into_iter();
        // a polygon with no rings takes no part
        if let Some(exterior) = rings.next() {
            tally_polygon!(place(exterior), rings, |hole| place(hole), inside, boundaries);
        }
    }
    tallied(inside, boundaries)
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
