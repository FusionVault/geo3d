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

// ── 0.2 additions ────────────────────────────────────────────────────────────────────────────────────

impl Vec3 {
    /// Distance between two points, metres.
    #[inline]
    pub fn distance_to(self, o: Vec3) -> f64 {
        (self - o).norm()
    }

    /// The angle between two vectors, radians in `[0, π]`; `0` if either is the zero vector.
    pub fn angle_to(self, o: Vec3) -> f64 {
        let d = self.norm() * o.norm();
        if d == 0.0 {
            return 0.0;
        }
        (self.dot(o) / d).clamp(-1.0, 1.0).acos()
    }

    /// The projection of this vector onto `onto` (the zero vector if `onto` is zero).
    #[inline]
    pub fn project_onto(self, onto: Vec3) -> Vec3 {
        let d = onto.dot(onto);
        if d == 0.0 {
            Vec3::ZERO
        } else {
            onto * (self.dot(onto) / d)
        }
    }
}

impl core::ops::Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn div(self, k: f64) -> Vec3 {
        Vec3::new(self.x / k, self.y / k, self.z / k)
    }
}

macro_rules! local_vector_ops {
    ($t:ident { $a:ident, $b:ident, $c:ident }) => {
        impl Add for $t {
            type Output = $t;
            #[inline]
            fn add(self, o: $t) -> $t {
                $t::new(self.$a + o.$a, self.$b + o.$b, self.$c + o.$c)
            }
        }
        impl Sub for $t {
            type Output = $t;
            #[inline]
            fn sub(self, o: $t) -> $t {
                $t::new(self.$a - o.$a, self.$b - o.$b, self.$c - o.$c)
            }
        }
        impl Mul<f64> for $t {
            type Output = $t;
            #[inline]
            fn mul(self, k: f64) -> $t {
                $t::new(self.$a * k, self.$b * k, self.$c * k)
            }
        }
        impl Neg for $t {
            type Output = $t;
            #[inline]
            fn neg(self) -> $t {
                $t::new(-self.$a, -self.$b, -self.$c)
            }
        }
        impl From<[f64; 3]> for $t {
            #[inline]
            fn from([a, b, c]: [f64; 3]) -> Self {
                $t::new(a, b, c)
            }
        }
        impl From<$t> for [f64; 3] {
            #[inline]
            fn from(v: $t) -> Self {
                [v.$a, v.$b, v.$c]
            }
        }
    };
}
local_vector_ops!(Enu { east, north, up });
local_vector_ops!(Ned { north, east, down });

impl Ned {
    /// As a plain vector `(north, east, down)`.
    #[inline]
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.north, self.east, self.down)
    }

    /// Azimuth / elevation / range of this vector.
    #[inline]
    pub fn to_aer(self) -> Aer {
        Aer::from_enu(self.to_enu())
    }
}

impl Aer {
    /// The pointing that reaches a NED vector.
    #[inline]
    pub fn from_ned(v: Ned) -> Aer {
        Aer::from_enu(v.to_enu())
    }

    /// The NED vector of this pointing at its range.
    #[inline]
    pub fn to_ned(self) -> Ned {
        self.to_enu().to_ned()
    }

    /// Whether the pointing is above the local horizon (elevation > 0).
    #[inline]
    pub fn is_above_horizon(self) -> bool {
        self.elevation_deg > 0.0
    }

    /// The same pointing with azimuth wrapped into `[0, 360)`.
    #[inline]
    pub fn normalized(self) -> Aer {
        Aer {
            azimuth_deg: self.azimuth_deg.rem_euclid(360.0),
            ..self
        }
    }
}

impl Geodetic {
    /// Longitude wrapped into `(-180, 180]` and latitude clamped into `[-90, 90]`; height unchanged.
    pub fn normalized(self) -> Geodetic {
        let mut lon = self.lon_deg.rem_euclid(360.0);
        if lon > 180.0 {
            lon -= 360.0;
        }
        Geodetic::new(self.lat_deg.clamp(-90.0, 90.0), lon, self.height_m)
    }

    /// The diametrically opposite point at the same height.
    pub fn antipode(self) -> Geodetic {
        Geodetic::new(-self.lat_deg, self.lon_deg + 180.0, self.height_m).normalized()
    }
}

/// An Earth-Centred Inertial position (metres): ECEF rotated by the sidereal angle, see [`crate::eci`].
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
mod tests_0_2 {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn vec3_extras() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 2.0, 0.0);
        assert!(close(a.angle_to(b), core::f64::consts::FRAC_PI_2, 1e-12));
        assert_eq!(a.angle_to(Vec3::ZERO), 0.0);
        assert_eq!(Vec3::new(3.0, 4.0, 0.0).project_onto(a), Vec3::new(3.0, 0.0, 0.0));
        assert_eq!(Vec3::new(1.0, 1.0, 1.0).project_onto(Vec3::ZERO), Vec3::ZERO);
        assert_eq!(Vec3::new(2.0, 4.0, 6.0) / 2.0, Vec3::new(1.0, 2.0, 3.0));
        assert!(close(Vec3::new(3.0, 4.0, 0.0).distance_to(Vec3::ZERO), 5.0, 1e-12));
        assert_eq!(Vec3::new(1.0, 2.0, 2.0).normalized().unwrap().norm(), 1.0);
        assert!(Vec3::ZERO.normalized().is_none());
        assert_eq!(<[f64; 3]>::from(Vec3::new(1.0, 2.0, 3.0)), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn local_vector_operators() {
        let e = Enu::new(1.0, 2.0, 3.0);
        assert_eq!(e + e, Enu::new(2.0, 4.0, 6.0));
        assert_eq!(e - e, Enu::default());
        assert_eq!(e * 2.0, Enu::new(2.0, 4.0, 6.0));
        assert_eq!(-e, Enu::new(-1.0, -2.0, -3.0));
        assert_eq!(Enu::from([1.0, 2.0, 3.0]), e);
        assert_eq!(<[f64; 3]>::from(e), [1.0, 2.0, 3.0]);
        let n = e.to_ned();
        assert_eq!(n, Ned::new(2.0, 1.0, -3.0));
        assert_eq!(n + n, Ned::new(4.0, 2.0, -6.0));
        assert_eq!(n.vec(), Vec3::new(2.0, 1.0, -3.0));
        assert_eq!(n.to_enu(), e);
        assert!(close(n.norm(), e.norm(), 1e-12));
    }

    #[test]
    fn aer_ned_and_horizon() {
        let aer = Aer::new(45.0, 10.0, 1000.0);
        let back = Aer::from_ned(aer.to_ned());
        assert!(
            close(back.azimuth_deg, 45.0, 1e-9)
                && close(back.elevation_deg, 10.0, 1e-9)
                && close(back.range_m, 1000.0, 1e-9)
        );
        assert_eq!(aer.to_ned().to_aer().normalized().azimuth_deg.round(), 45.0);
        assert!(aer.is_above_horizon() && !Aer::new(0.0, -1.0, 1.0).is_above_horizon());
        assert_eq!(Aer::new(-90.0, 0.0, 1.0).normalized().azimuth_deg, 270.0);
        assert_eq!(Aer::from_enu(Enu::default()), Aer::new(0.0, 0.0, 0.0));
    }

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
        let i = Eci::from(Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(i.vec(), Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(i.to_array(), [1.0, 2.0, 3.0]);
    }
}
