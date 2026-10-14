//! Offsets and pointings from an observer: [`Enu`] (east, north, up), [`Ned`] (north, east, down)
//! and [`Aer`] (azimuth, elevation, range).

use crate::vector::{vector_ops, Vec3};

/// An offset in a local **East-North-Up** frame, metres: how far east, north and up something is
/// from an observer. Get one from a [`LocalFrame`](crate::LocalFrame). Supports `+`, `-`, `* f64`
/// and unary `-`.
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

    /// The same offset in North-East-Down.
    #[inline]
    pub const fn to_ned(self) -> Ned {
        Ned::new(self.north, self.east, -self.up)
    }

    /// Azimuth, elevation and range of this offset (see [`Aer::from_enu`]).
    #[inline]
    pub fn to_aer(self) -> Aer {
        Aer::from_enu(self)
    }
}

/// An offset in a local **North-East-Down** frame, metres — the aerospace convention. The same as
/// [`Enu`] with the axes reordered and up flipped to down. Supports `+`, `-`, `* f64` and unary `-`.
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

    /// Length, metres.
    #[inline]
    pub fn norm(self) -> f64 {
        self.to_enu().norm()
    }

    /// As a plain vector `(north, east, down)`.
    #[inline]
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.north, self.east, self.down)
    }

    /// The same offset in East-North-Up.
    #[inline]
    pub const fn to_enu(self) -> Enu {
        Enu::new(self.east, self.north, -self.down)
    }

    /// Azimuth, elevation and range of this offset.
    #[inline]
    pub fn to_aer(self) -> Aer {
        Aer::from_enu(self.to_enu())
    }
}

vector_ops!(Enu { east, north, up });
vector_ops!(Ned { north, east, down });

/// Where to point from an observer: **azimuth** in degrees clockwise from north, **elevation** in
/// degrees above the horizon, and **range** in metres. The look angles a radar, antenna or
/// telescope uses.
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

    /// The pointing that reaches an ENU offset; azimuth in `[0, 360)`. A zero offset has azimuth 0
    /// and elevation 0. Elevation keeps full precision right up to the zenith.
    #[doc(alias = "enu2aer")]
    #[inline]
    pub fn from_enu(v: Enu) -> Aer {
        let azimuth_deg = v.east.atan2(v.north).to_degrees().rem_euclid(360.0);
        // atan2 rather than asin(up / range), which loses precision near ±90°. sqrt(e² + n²) is as
        // accurate as hypot at any distance in metres (no overflow below 1e154) and much cheaper.
        let elevation_deg = v.up.atan2((v.east * v.east + v.north * v.north).sqrt()).to_degrees();
        Aer {
            azimuth_deg,
            elevation_deg,
            range_m: v.norm(),
        }
    }

    /// The pointing that reaches a NED offset.
    #[inline]
    pub fn from_ned(v: Ned) -> Aer {
        Aer::from_enu(v.to_enu())
    }

    /// The ENU offset this pointing reaches at its range.
    #[doc(alias = "aer2enu")]
    #[inline]
    pub fn to_enu(self) -> Enu {
        let d = self.direction_enu();
        Enu::new(d.east * self.range_m, d.north * self.range_m, d.up * self.range_m)
    }

    /// The NED offset this pointing reaches at its range.
    #[inline]
    pub fn to_ned(self) -> Ned {
        self.to_enu().to_ned()
    }

    /// The unit direction of this pointing in ENU (range ignored).
    #[inline]
    pub fn direction_enu(self) -> Enu {
        let (saz, caz) = self.azimuth_deg.to_radians().sin_cos();
        let (sel, cel) = self.elevation_deg.to_radians().sin_cos();
        Enu::new(saz * cel, caz * cel, sel)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::close;

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
    fn elevation_keeps_its_precision_near_the_zenith() {
        // 1 nm off vertical at 1 m range: asin(up / range) rounds to exactly 90°.
        let el = Enu::new(1e-9, 0.0, 1.0).to_aer().elevation_deg;
        assert!(close(el, 90.0 - 1e-9_f64.to_degrees(), 1e-13), "{el}");
        assert_eq!(Enu::new(0.0, 0.0, -5.0).to_aer().elevation_deg, -90.0);
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
        assert!(close(Aer::new(90.0, 0.0, 1_000.0).to_enu().east, 1_000.0, 1e-9));
    }
}
