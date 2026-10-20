//! Positions: [`Geodetic`] (latitude, longitude, height), [`Ecef`] (x, y, z turning with the Earth)
//! and [`Eci`] (x, y, z fixed in space).

use crate::vector::Vec3;

/// A position as latitude and longitude in **degrees** and height in **metres** — what a GPS
/// receiver reports.
///
/// The height is above the ellipsoid, not above mean sea level; the two differ by up to about
/// 100 m depending on where you are. A [`Geoid`](crate::Geoid) converts between them.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Geodetic {
    /// Latitude, degrees, north positive.
    pub lat_deg: f64,
    /// Longitude, degrees, east positive.
    pub lon_deg: f64,
    /// Height above the ellipsoid, metres.
    pub height_m: f64,
}

impl Geodetic {
    /// Build from latitude and longitude in degrees and height in metres. Latitude comes first.
    #[inline]
    pub const fn new(lat_deg: f64, lon_deg: f64, height_m: f64) -> Self {
        Geodetic {
            lat_deg,
            lon_deg,
            height_m,
        }
    }

    /// Build from longitude first, then latitude, in degrees, and height in metres: the order of
    /// GeoJSON, WKT and GeoArrow, whose `x` is the longitude. Use it at those boundaries, so the two
    /// cannot be swapped by accident.
    ///
    /// ```
    /// use geo3d::Geodetic;
    /// assert_eq!(Geodetic::from_lon_lat(151.2, -33.9, 10.0), Geodetic::new(-33.9, 151.2, 10.0));
    /// ```
    #[inline]
    pub const fn from_lon_lat(lon_deg: f64, lat_deg: f64, height_m: f64) -> Self {
        Geodetic::new(lat_deg, lon_deg, height_m)
    }

    /// Latitude in radians.
    #[inline]
    pub fn lat_rad(self) -> f64 {
        self.lat_deg.to_radians()
    }

    /// Longitude in radians.
    #[inline]
    pub fn lon_rad(self) -> f64 {
        self.lon_deg.to_radians()
    }

    /// The same position at a different height.
    #[inline]
    pub const fn with_height(self, height_m: f64) -> Self {
        Geodetic { height_m, ..self }
    }

    /// Longitude wrapped into `(-180, 180]` and latitude clamped into `[-90, 90]`; height unchanged.
    pub fn normalized(self) -> Geodetic {
        // Exact: a longitude already in range comes back unchanged.
        let lon = crate::angle::ang_normalize(self.lon_deg);
        Geodetic::new(self.lat_deg.clamp(-90.0, 90.0), lon, self.height_m)
    }

    /// The point on the opposite side of the Earth, at the same height.
    pub fn antipode(self) -> Geodetic {
        Geodetic::new(-self.lat_deg, self.lon_deg + 180.0, self.height_m).normalized()
    }
}

/// A position in **Earth-Centred, Earth-Fixed** coordinates, metres. The origin is the Earth's
/// centre, `x` points to latitude 0° longitude 0°, `y` to longitude 90° east and `z` to the North
/// Pole. The axes turn with the Earth, so a parked car keeps the same `Ecef`.
///
/// Straight lines, distances and angles work here like in any 3D space. It is also the frame
/// CesiumJS renders in.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ecef {
    /// X, metres.
    pub x: f64,
    /// Y, metres.
    pub y: f64,
    /// Z, metres.
    pub z: f64,
}

impl Ecef {
    /// Build from components in metres.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Ecef { x, y, z }
    }

    /// As a plain vector.
    #[inline]
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    /// Straight-line distance to another point, metres.
    #[inline]
    pub fn distance_to(self, o: Ecef) -> f64 {
        (self.vec() - o.vec()).norm()
    }

    /// Move by a vector (metres).
    #[inline]
    pub fn offset(self, d: Vec3) -> Ecef {
        Ecef::from(self.vec() + d)
    }

    /// As `[x, y, z]`.
    #[inline]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl From<Vec3> for Ecef {
    #[inline]
    fn from(v: Vec3) -> Self {
        Ecef::new(v.x, v.y, v.z)
    }
}
impl From<[f64; 3]> for Ecef {
    #[inline]
    fn from([x, y, z]: [f64; 3]) -> Self {
        Ecef::new(x, y, z)
    }
}
impl From<Ecef> for [f64; 3] {
    #[inline]
    fn from(e: Ecef) -> Self {
        e.to_array()
    }
}

/// A position in **Earth-Centred Inertial** coordinates, metres: like [`Ecef`], but the axes do not
/// turn with the Earth. This is the TEME frame that SGP4 / TLE satellite propagators output, not
/// J2000 / GCRF. Convert with [`Ecef::to_eci`] / [`Eci::to_ecef`]; see the [`eci`](crate::eci) module.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Eci {
    /// X, metres (towards the vernal equinox).
    pub x: f64,
    /// Y, metres.
    pub y: f64,
    /// Z, metres (along the rotation axis).
    pub z: f64,
}

impl Eci {
    /// Build from components in metres.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Eci { x, y, z }
    }

    /// As a plain vector.
    #[inline]
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    /// As `[x, y, z]`.
    #[inline]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl From<Vec3> for Eci {
    #[inline]
    fn from(v: Vec3) -> Self {
        Eci::new(v.x, v.y, v.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::close;

    #[test]
    fn geodetic_normalisation_and_antipode() {
        let g = Geodetic::new(95.0, 190.0, 5.0).normalized();
        assert_eq!(g, Geodetic::new(90.0, -170.0, 5.0));
        assert_eq!(Geodetic::new(0.0, 180.0, 0.0).normalized().lon_deg, 180.0);
        assert_eq!(Geodetic::new(0.0, -180.0, 0.0).normalized().lon_deg, 180.0);
        assert_eq!(
            Geodetic::new(10.0, 20.0, 1.0).antipode(),
            Geodetic::new(-10.0, -160.0, 1.0)
        );
        assert_eq!(
            Geodetic::new(10.0, 20.0, 1.0).antipode().antipode(),
            Geodetic::new(10.0, 20.0, 1.0)
        );
        assert_eq!(Geodetic::new(1.0, 2.0, 3.0).with_height(9.0).height_m, 9.0);
        assert!(close(
            Geodetic::new(90.0, 180.0, 0.0).lat_rad(),
            core::f64::consts::FRAC_PI_2,
            1e-15
        ));
        assert!(close(
            Geodetic::new(90.0, 180.0, 0.0).lon_rad(),
            core::f64::consts::PI,
            1e-15
        ));
    }

    #[test]
    fn ecef_and_eci_helpers() {
        let e = Ecef::from([1.0, 2.0, 3.0]);
        assert_eq!(<[f64; 3]>::from(e), [1.0, 2.0, 3.0]);
        assert_eq!(e.offset(Vec3::new(1.0, 1.0, 1.0)), Ecef::new(2.0, 3.0, 4.0));
        assert!(close(e.distance_to(Ecef::new(1.0, 2.0, 0.0)), 3.0, 1e-12));
        let i = Eci::from(Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(i.vec(), Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(i.to_array(), [1.0, 2.0, 3.0]);
    }

    #[cfg(feature = "serde")]
    #[test]
    fn coordinate_types_round_trip_through_json() {
        use crate::Aer;
        let g = Geodetic::new(34.05, -118.24, 100.0);
        let json = serde_json::to_string(&g).unwrap();
        assert_eq!(json, r#"{"lat_deg":34.05,"lon_deg":-118.24,"height_m":100.0}"#);
        assert_eq!(serde_json::from_str::<Geodetic>(&json).unwrap(), g);
        let aer = Aer::new(90.0, 5.0, 1000.0);
        assert_eq!(
            serde_json::from_str::<Aer>(&serde_json::to_string(&aer).unwrap()).unwrap(),
            aer
        );
    }
}
