//! UTM and UPS: the world grid of 60 transverse Mercator zones plus two polar stereographic caps,
//! with the zone rules and limits of NGA.STND.0037 as GeographicLib applies them.

use super::{PolarStereographic, TransverseMercator};
use crate::angle::{ang_diff, ang_normalize};
use crate::ellipsoid::{Ellipsoid, WGS84};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;

/// A position on the UTM/UPS grid: a zone, a hemisphere, and an easting and northing in metres.
///
/// - **UTM** covers 80° S to 84° N in 60 zones, each 6° of longitude wide (with the Norway and
///   Svalbard exceptions). `zone` is 1 to 60. Eastings are 500 000 m on the zone's central
///   meridian; northings count from the equator, plus 10 000 000 m in the south.
/// - **UPS** covers the polar caps. `zone` is 0. Both coordinates are 2 000 000 m at the pole.
///
/// Height isn't part of a grid position: set it again with [`Geodetic::with_height`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct UtmUps {
    /// The UTM zone, 1 to 60, or 0 for UPS.
    pub zone: u8,
    /// Northern hemisphere.
    pub north: bool,
    /// Easting, metres.
    pub easting_m: f64,
    /// Northing, metres.
    pub northing_m: f64,
}

/// 100 km: the grid limits are whole multiples of it, and positions up to one tile beyond them are
/// still accepted (GeographicLib's "slop").
const TILE: f64 = 100_000.0;

impl UtmUps {
    /// Build from a zone (0 for UPS), hemisphere, easting and northing.
    #[inline]
    pub const fn new(zone: u8, north: bool, easting_m: f64, northing_m: f64) -> Self {
        UtmUps {
            zone,
            north,
            easting_m,
            northing_m,
        }
    }

    /// The zone's central meridian, degrees (NaN for UPS).
    pub fn central_meridian_deg(self) -> f64 {
        if self.zone == 0 {
            f64::NAN
        } else {
            6.0 * self.zone as f64 - 183.0
        }
    }

    /// The zone UTM/UPS uses for a position: 0 (UPS) beyond 84° N and 80° S, otherwise the UTM
    /// zone, including the Norway (32V) and Svalbard (31X, 33X, 35X, 37X) exceptions. `None` for
    /// a NaN or infinite latitude or longitude, or a latitude beyond ±90°.
    pub fn standard_zone(g: Geodetic) -> Option<u8> {
        let (lat, lon) = (g.lat_deg, g.lon_deg);
        if !(lat.abs() <= 90.0 && lon.is_finite()) {
            return None;
        }
        if !(-80.0..84.0).contains(&lat) {
            return Some(0);
        }
        let mut ilon = ang_normalize(lon).floor() as i32;
        if ilon == 180 {
            ilon = -180; // now in [-180, 180)
        }
        let mut zone = (ilon + 186) / 6;
        let band = ((lat.floor() as i32 + 80) / 8 - 10).clamp(-10, 9);
        if band == 7 && zone == 31 && ilon >= 3 {
            zone = 32; // Norway
        } else if band == 9 && (0..42).contains(&ilon) {
            zone = 2 * ((ilon + 183) / 12) + 1; // Svalbard
        }
        Some(zone as u8)
    }

    /// Whether an easting and northing are within the grid's limits for this zone kind and
    /// hemisphere, give or take 100 km.
    fn in_range(utm: bool, north: bool, easting_m: f64, northing_m: f64) -> bool {
        // (min easting, max easting, min northing, max northing) in 100 km tiles.
        let (e0, e1, n0, n1) = match (utm, north) {
            (false, false) => (8.0, 32.0, 8.0, 32.0),
            (false, true) => (13.0, 27.0, 13.0, 27.0),
            (true, false) => (1.0, 9.0, 10.0, 100.0 + 95.0),
            (true, true) => (1.0, 9.0, 10.0 - 100.0, 95.0),
        };
        // Written so that a NaN fails.
        (e0 - 1.0) * TILE <= easting_m
            && easting_m <= (e1 + 1.0) * TILE
            && (n0 - 1.0) * TILE <= northing_m
            && northing_m <= (n1 + 1.0) * TILE
    }
}

impl Ellipsoid {
    /// The UTM/UPS grid position of `g` (height ignored), in its standard zone. Accurate to
    /// about 5 nm.
    ///
    /// `None` for a NaN or infinite latitude or longitude, or a latitude beyond ±90°.
    ///
    /// ```
    /// use geo3d::{Geodetic, UtmUps, WGS84};
    /// let u = WGS84.to_utm_ups(Geodetic::new(-33.87, 151.21, 20.0)).unwrap();   // Sydney
    /// assert_eq!((u.zone, u.north), (56, false));
    /// assert!((u.easting_m - 334_435.706_141_829).abs() < 1e-8);
    /// assert!((u.northing_m - 6_250_816.397_804_656).abs() < 1e-8);
    /// ```
    #[doc(alias = "utm")]
    #[doc(alias = "ups")]
    pub fn to_utm_ups(self, g: Geodetic) -> Option<UtmUps> {
        self.to_utm_ups_zone(g, UtmUps::standard_zone(g)?)
    }

    /// The UTM/UPS grid position of `g` in a chosen zone (1 to 60, or 0 for UPS), for keeping
    /// nearby data in one zone. The hemisphere follows the latitude.
    ///
    /// `None` if the position is too far from the zone: its easting or northing would be more than
    /// 100 km beyond the grid's limits (for UTM, an easting from 0 to 1 000 km). Also `None` for a
    /// zone above 60, and for the inputs [`to_utm_ups`](Self::to_utm_ups) refuses.
    pub fn to_utm_ups_zone(self, g: Geodetic, zone: u8) -> Option<UtmUps> {
        if !(g.lat_deg.abs() <= 90.0 && g.lon_deg.is_finite()) || zone > 60 {
            return None;
        }
        let north = !g.lat_deg.is_sign_negative();
        let (easting_m, northing_m) = if zone == 0 {
            if g.lat_deg.abs() < 70.0 {
                return None;
            }
            let p = self.ups().forward(north, g);
            (p.x_m + 20.0 * TILE, p.y_m + 20.0 * TILE)
        } else {
            let lon0 = 6.0 * zone as f64 - 183.0;
            let dlon = ang_diff(lon0, g.lon_deg).0;
            if dlon.is_nan() || dlon.abs() > 60.0 {
                return None;
            }
            let p = self.utm().project::<false>(lon0, g);
            (p.x_m + 5.0 * TILE, p.y_m + if north { 0.0 } else { 100.0 * TILE })
        };
        UtmUps::in_range(zone != 0, north, easting_m, northing_m).then_some(UtmUps {
            zone,
            north,
            easting_m,
            northing_m,
        })
    }

    /// The position (height 0) of a UTM/UPS grid position: the reverse of
    /// [`to_utm_ups`](Self::to_utm_ups).
    ///
    /// `None` for a zone above 60, a NaN easting or northing, or one more than 100 km beyond the
    /// grid's limits.
    ///
    /// ```
    /// use geo3d::{UtmUps, WGS84};
    /// let g = WGS84.from_utm_ups(UtmUps::new(56, false, 334_786.0, 6_252_080.0)).unwrap();
    /// assert!((g.lat_deg + 33.858_663_916_529_91).abs() < 1e-12);
    /// assert!((g.lon_deg - 151.214_022_860_895_08).abs() < 1e-12);
    /// ```
    pub fn from_utm_ups(self, u: UtmUps) -> Option<Geodetic> {
        if u.zone > 60 || !UtmUps::in_range(u.zone != 0, u.north, u.easting_m, u.northing_m) {
            return None;
        }
        Some(if u.zone == 0 {
            self.ups()
                .reverse(u.north, u.easting_m - 20.0 * TILE, u.northing_m - 20.0 * TILE)
        } else {
            let y = u.northing_m - if u.north { 0.0 } else { 100.0 * TILE };
            self.utm()
                .reverse(u.central_meridian_deg(), u.easting_m - 5.0 * TILE, y)
        })
    }

    /// UTM's transverse Mercator on this ellipsoid (WGS84's is built at compile time).
    fn utm(self) -> TransverseMercator {
        if self == WGS84 {
            TransverseMercator::UTM
        } else {
            TransverseMercator::new(self, 0.9996)
        }
    }

    fn ups(self) -> PolarStereographic {
        PolarStereographic::new(self, 0.994)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zone_edges_and_exceptions() {
        let z = |lat, lon| UtmUps::standard_zone(Geodetic::new(lat, lon, 0.0)).unwrap();
        // Ordinary zones: 6° wide, zone 1 starting at 180° W; 180° E is zone 1 too.
        assert_eq!((z(10.0, -180.0), z(10.0, -174.000_000_1), z(10.0, -174.0)), (1, 1, 2));
        assert_eq!((z(10.0, 179.999), z(10.0, 180.0), z(10.0, 0.0)), (60, 1, 31));
        // Norway: band V (56° to 64° N) from 3° E is zone 32.
        assert_eq!(
            (
                z(56.0, 3.0),
                z(55.999, 3.0),
                z(63.999, 3.0),
                z(64.0, 3.0),
                z(60.0, 2.999)
            ),
            (32, 31, 32, 31, 31)
        );
        // Svalbard: band X (72° to 84° N), 0° to 42° E, has only the odd zones 31, 33, 35 and 37.
        assert_eq!(
            (
                z(72.0, 8.999),
                z(72.0, 9.0),
                z(72.0, 21.0),
                z(72.0, 33.0),
                z(72.0, 42.0)
            ),
            (31, 33, 35, 37, 38)
        );
        assert_eq!(z(71.999, 9.0), 32);
        // UPS beyond 84° N and 80° S.
        assert_eq!(
            (z(83.999, 0.0), z(84.0, 0.0), z(-80.0, 0.0), z(-80.000_1, 0.0)),
            (31, 0, 31, 0)
        );
    }
}
