//! The polar stereographic projection, exact (no series), after GeographicLib.

use super::MapPoint;
use crate::angle::{ang_normalize, atan2d, atand, lat_fix, sincosd, tand};
use crate::auxlat::{eatanhe, signed_eccentricity, tauf, taupf};
use crate::ellipsoid::{Ellipsoid, WGS84};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;

/// The polar stereographic projection on an ellipsoid, centred on a pole with scale `k0` there.
/// It's conformal and exact everywhere in its hemisphere. UPS, the polar half of UTM/UPS, is this
/// projection with `k0 = 0.994`, as [`PolarStereographic::UPS`].
///
/// It gives plain x and y from the pole, without the false easting and northing a grid adds. In
/// the north, y points away from 180° E; in the south, towards it.
///
/// ```
/// use geo3d::{Geodetic, PolarStereographic};
/// let p = PolarStereographic::UPS.forward(true, Geodetic::new(85.0, 45.0, 0.0));
/// let back = PolarStereographic::UPS.reverse(true, p.x_m, p.y_m);
/// assert!((back.lat_deg - 85.0).abs() < 1e-13 && (back.lon_deg - 45.0).abs() < 1e-13);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarStereographic {
    a: f64,
    f: f64,
    k0: f64,
    e2: f64,
    e2m: f64,
}

impl PolarStereographic {
    /// UPS's projection: WGS84 with `k0 = 0.994`.
    pub const UPS: PolarStereographic = PolarStereographic::new(WGS84, 0.994);

    /// The projection on `ellipsoid` with scale `k0` at the pole (positive).
    pub const fn new(ellipsoid: Ellipsoid, k0: f64) -> Self {
        let f = ellipsoid.f();
        let e2 = f * (2.0 - f);
        PolarStereographic {
            a: ellipsoid.a(),
            f,
            k0,
            e2,
            e2m: 1.0 - e2,
        }
    }

    /// The scale at the pole.
    #[inline]
    pub const fn k0(&self) -> f64 {
        self.k0
    }

    fn es(&self) -> f64 {
        signed_eccentricity(self.e2, self.f)
    }

    /// `2·k0·a / c`, where `c = (1 − f)·exp(e·atanh e)`: the distance from the pole per unit of
    /// the projection's radial function.
    fn scale_length(&self) -> f64 {
        2.0 * self.k0 * self.a / ((1.0 - self.f) * eatanhe(1.0, self.es()).exp())
    }

    /// Project `g` (height ignored) about the north pole (`north`) or the south pole.
    pub fn forward(&self, north: bool, g: Geodetic) -> MapPoint {
        let lat = lat_fix(g.lat_deg) * if north { 1.0 } else { -1.0 };
        let tau = tand(lat);
        let secphi = 1f64.hypot(tau);
        let taup = taupf(tau, self.es());
        let rho = 1f64.hypot(taup) + taup.abs();
        // Distance from the pole, in units of scale_length: exactly 0 at the pole.
        let rho = if taup >= 0.0 {
            if lat != 90.0 {
                1.0 / rho
            } else {
                0.0
            }
        } else {
            rho
        } * self.scale_length();
        let k = if lat != 90.0 {
            (rho / self.a) * secphi * (self.e2m + self.e2 / (secphi * secphi)).sqrt()
        } else {
            self.k0
        };
        let (slon, clon) = sincosd(g.lon_deg);
        MapPoint {
            x_m: slon * rho,
            y_m: clon * if north { -rho } else { rho },
            convergence_deg: ang_normalize(if north { g.lon_deg } else { -g.lon_deg }),
            scale: k,
        }
    }

    /// The position (height 0) at `x_m`, `y_m` from the north pole (`north`) or the south pole:
    /// the reverse of [`forward`](Self::forward).
    pub fn reverse(&self, north: bool, x_m: f64, y_m: f64) -> Geodetic {
        self.reverse_full(north, x_m, y_m).0
    }

    /// [`reverse`](Self::reverse), with the grid convergence (degrees) and scale there.
    pub(crate) fn reverse_full(&self, north: bool, x_m: f64, y_m: f64) -> (Geodetic, f64, f64) {
        let rho = x_m.hypot(y_m);
        let t = if rho != 0.0 {
            rho / self.scale_length()
        } else {
            f64::EPSILON * f64::EPSILON
        };
        let taup = (1.0 / t - t) / 2.0;
        let tau = tauf(taup, self.es());
        let secphi = 1f64.hypot(tau);
        let k = if rho != 0.0 {
            (rho / self.a) * secphi * (self.e2m + self.e2 / (secphi * secphi)).sqrt()
        } else {
            self.k0
        };
        let lat = if north { 1.0 } else { -1.0 } * atand(tau);
        let lon = atan2d(x_m, if north { -y_m } else { y_m });
        let gamma = ang_normalize(if north { lon } else { -lon });
        (Geodetic::new(lat, lon, 0.0), gamma, k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_over_both_caps() {
        let ps = PolarStereographic::UPS;
        for i in 0..=200 {
            for j in 0..36 {
                let lat = 60.0 + 30.0 * i as f64 / 200.0;
                let lon = -175.0 + 10.0 * j as f64;
                for north in [true, false] {
                    let g = Geodetic::new(if north { lat } else { -lat }, lon, 0.0);
                    let p = ps.forward(north, g);
                    let (back, gamma, k) = ps.reverse_full(north, p.x_m, p.y_m);
                    assert!((back.lat_deg - g.lat_deg).abs() < 1e-13, "{g:?} -> {back:?}");
                    if lat < 90.0 {
                        assert!((back.lon_deg - lon).abs() < 1e-11, "{g:?} -> {back:?}");
                        assert!((gamma - p.convergence_deg).abs() < 1e-11 && (k - p.scale).abs() < 1e-14);
                    }
                }
            }
        }
    }

    #[test]
    fn the_pole_and_the_scale() {
        let ps = PolarStereographic::UPS;
        let p = ps.forward(true, Geodetic::new(90.0, 123.0, 0.0));
        assert_eq!((p.x_m, p.y_m, p.scale), (0.0, -0.0, 0.994));
        assert!(ps.reverse(false, 0.0, 0.0).lat_deg == -90.0);
        // Away from the pole the scale grows: about 1.0 by 81° (UPS's design) and more beyond.
        let k81 = ps.forward(true, Geodetic::new(81.07, 0.0, 0.0)).scale;
        assert!((k81 - 1.0).abs() < 1e-3, "{k81}");
        // North: y points away from 180° E, so the 180° meridian is +y; south, the opposite.
        assert!(ps.forward(true, Geodetic::new(85.0, 180.0, 0.0)).y_m > 0.0);
        assert!(ps.forward(false, Geodetic::new(-85.0, 180.0, 0.0)).y_m < 0.0);
    }
}
