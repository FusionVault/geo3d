//! The coordinate types. All are `Copy`, plain `f64` fields, no invariants beyond units.

use core::ops::{Add, Mul, Neg, Sub};

/// A Cartesian 3-vector (metres, or unitless for directions).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };

    /// Build from components.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product (right-handed).
    #[inline]
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Euclidean length.
    #[inline]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// The unit vector in this direction, or `None` for the zero vector.
    #[inline]
    pub fn normalized(self) -> Option<Vec3> {
        let n = self.norm();
        (n > 0.0).then(|| self * (1.0 / n))
    }

    /// As `[x, y, z]`.
    #[inline]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl From<[f64; 3]> for Vec3 {
    #[inline]
    fn from([x, y, z]: [f64; 3]) -> Self {
        Vec3 { x, y, z }
    }
}
impl From<Vec3> for [f64; 3] {
    #[inline]
    fn from(v: Vec3) -> Self {
        v.to_array()
    }
}
impl Add for Vec3 {
    type Output = Vec3;
    #[inline]
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl Sub for Vec3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Mul<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, k: f64) -> Vec3 {
        Vec3::new(self.x * k, self.y * k, self.z * k)
    }
}
impl Neg for Vec3 {
    type Output = Vec3;
    #[inline]
    fn neg(self) -> Vec3 {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

/// A geodetic position: latitude and longitude in **degrees**, height in **metres** above the
/// ellipsoid (ellipsoidal height, not orthometric/MSL).
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
    /// Build from latitude and longitude in degrees and height in metres.
    #[inline]
    pub const fn new(lat_deg: f64, lon_deg: f64, height_m: f64) -> Self {
        Geodetic {
            lat_deg,
            lon_deg,
            height_m,
        }
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
}

/// An Earth-Centred, Earth-Fixed Cartesian position, metres. `x` through the equator/prime meridian,
/// `z` through the north pole.
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

    /// Straight-line distance to another ECEF point, metres.
    #[inline]
    pub fn distance_to(self, o: Ecef) -> f64 {
        (self.vec() - o.vec()).norm()
    }

    /// Translate by a vector (metres).
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

/// A vector in a local **East-North-Up** frame, metres (see [`LocalFrame`](crate::LocalFrame)).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Enu {
    /// East, metres.
    pub east: f64,
    /// North, metres.
    pub north: f64,
    /// Up (along the ellipsoid normal), metres.
    pub up: f64,
}

impl Enu {
    /// Build from components in metres.
    #[inline]
    pub const fn new(east: f64, north: f64, up: f64) -> Self {
        Enu { east, north, up }
    }

    /// Length, metres.
    #[inline]
    pub fn norm(self) -> f64 {
        self.vec().norm()
    }

    /// As a plain vector `(east, north, up)`.
    #[inline]
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.east, self.north, self.up)
    }

    /// The same vector in North-East-Down.
    #[inline]
    pub const fn to_ned(self) -> Ned {
        Ned::new(self.north, self.east, -self.up)
    }

    /// Azimuth / elevation / range of this vector (see [`Aer::from_enu`]).
    #[inline]
    pub fn to_aer(self) -> Aer {
        Aer::from_enu(self)
    }
}

/// A vector in a local **North-East-Down** frame, metres (the aerospace convention).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Ned {
    /// North, metres.
    pub north: f64,
    /// East, metres.
    pub east: f64,
    /// Down, metres.
    pub down: f64,
}

impl Ned {
    /// Build from components in metres.
    #[inline]
    pub const fn new(north: f64, east: f64, down: f64) -> Self {
        Ned { north, east, down }
    }

    /// The same vector in East-North-Up.
    #[inline]
    pub const fn to_enu(self) -> Enu {
        Enu::new(self.east, self.north, -self.down)
    }

    /// Length, metres.
    #[inline]
    pub fn norm(self) -> f64 {
        self.to_enu().norm()
    }
}

/// Azimuth-Elevation-Range: a pointing from an observer. Azimuth in degrees from north, clockwise,
/// normalised to `[0, 360)`; elevation in degrees above the local horizon; range in metres.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Aer {
    /// Azimuth, degrees clockwise from north.
    pub azimuth_deg: f64,
    /// Elevation, degrees above the horizon.
    pub elevation_deg: f64,
    /// Slant range, metres.
    pub range_m: f64,
}

impl Aer {
    /// Build from azimuth and elevation in degrees and range in metres.
    #[inline]
    pub const fn new(azimuth_deg: f64, elevation_deg: f64, range_m: f64) -> Self {
        Aer {
            azimuth_deg,
            elevation_deg,
            range_m,
        }
    }

    /// The unit direction of this pointing in the local ENU frame (range ignored).
    #[inline]
    pub fn direction_enu(self) -> Enu {
        let (saz, caz) = self.azimuth_deg.to_radians().sin_cos();
        let (sel, cel) = self.elevation_deg.to_radians().sin_cos();
        Enu::new(saz * cel, caz * cel, sel)
    }

    /// The ENU vector of this pointing at its range.
    #[inline]
    pub fn to_enu(self) -> Enu {
        let d = self.direction_enu();
        Enu::new(d.east * self.range_m, d.north * self.range_m, d.up * self.range_m)
    }

    /// The pointing that reaches an ENU vector. A zero vector has azimuth 0 and elevation 0.
    #[inline]
    pub fn from_enu(v: Enu) -> Aer {
        let range = v.norm();
        let azimuth_deg = v.east.atan2(v.north).to_degrees().rem_euclid(360.0);
        let elevation_deg = if range == 0.0 {
            0.0
        } else {
            (v.up / range).asin().to_degrees()
        };
        Aer {
            azimuth_deg,
            elevation_deg,
            range_m: range,
        }
    }
}

#[cfg(all(test, feature = "serde"))]
mod serde_tests {
    use super::*;

    #[test]
    fn coordinate_types_round_trip_through_json() {
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
