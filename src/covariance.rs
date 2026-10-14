//! [`Cov3`]: how uncertain a 3D position is, carried between frames.

use crate::ellipsoid::Ellipsoid;
use crate::frame::LocalFrame;
use crate::local::{Aer, Enu};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;

const DEG: f64 = core::f64::consts::PI / 180.0;

/// A 3×3 covariance: how uncertain a 3D position (or velocity) is, and how its errors are linked.
///
/// Rows and columns follow the coordinates' own order: x, y, z for ECEF; east, north, up for ENU;
/// north, east, down for NED; azimuth, elevation, range for AER; latitude, longitude, height for
/// geodetic. Units are those of the coordinates, squared: m² for metres, deg² for degrees.
///
/// Make one from standard deviations, move it between frames, and read the standard deviations
/// back:
///
/// | From → to | Use | |
/// |---|---|---|
/// | ENU ↔ ECEF | [`LocalFrame::enu_cov_to_ecef`], [`LocalFrame::ecef_cov_to_enu`] | exact |
/// | ENU ↔ NED | [`Cov3::enu_to_ned`], [`Cov3::ned_to_enu`] | exact |
/// | AER ↔ ENU | [`Aer::cov_to_enu`], [`Enu::cov_to_aer`] | first order |
/// | geodetic ↔ ECEF | [`Ellipsoid::geodetic_cov_to_ecef`], [`Ellipsoid::ecef_cov_to_geodetic`] | first order |
///
/// "First order" means the conversion is treated as straight (linear) across the spread of the
/// errors, as Orekit and Stone Soup do. That's accurate while the errors are small next to the
/// distances involved. At long range with large angular noise (several degrees), the true spread
/// curves into a banana shape that no 3×3 covariance can describe.
///
/// ```
/// use geo3d::{Aer, Cov3};
/// // A radar return: 0.2° in azimuth and elevation, 20 m in range, at 50 km.
/// let hit = Aer::new(40.0, 10.0, 50_000.0);
/// let enu = hit.cov_to_enu(Cov3::from_std_devs([0.2, 0.2, 20.0]));
/// let [east, north, up] = enu.std_devs();   // 134 m, 114 m, 172 m
/// assert!((up - 171.9).abs() < 0.1);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cov3(pub [[f64; 3]; 3]);

impl Cov3 {
    /// Independent errors with these standard deviations: the variances `σ²` on the diagonal,
    /// zero elsewhere.
    pub const fn from_std_devs(s: [f64; 3]) -> Self {
        Cov3([
            [s[0] * s[0], 0.0, 0.0],
            [0.0, s[1] * s[1], 0.0],
            [0.0, 0.0, s[2] * s[2]],
        ])
    }

    /// The standard deviation of each coordinate: the square roots of the diagonal.
    pub fn std_devs(self) -> [f64; 3] {
        [self.0[0][0].sqrt(), self.0[1][1].sqrt(), self.0[2][2].sqrt()]
    }

    /// This ENU covariance in NED order. Exact.
    pub fn enu_to_ned(self) -> Cov3 {
        self.transformed([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
    }

    /// This NED covariance in ENU order. Exact.
    pub fn ned_to_enu(self) -> Cov3 {
        self.enu_to_ned() // the same swap both ways
    }

    /// `J · self · Jᵀ`, kept exactly symmetric.
    fn transformed(self, j: [[f64; 3]; 3]) -> Cov3 {
        let s = self.0;
        let mut js = [[0.0; 3]; 3];
        for (row, jr) in js.iter_mut().zip(&j) {
            for (c, out) in row.iter_mut().enumerate() {
                *out = jr[0] * s[0][c] + jr[1] * s[1][c] + jr[2] * s[2][c];
            }
        }
        let mut out = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in r..3 {
                out[r][c] = js[r][0] * j[c][0] + js[r][1] * j[c][1] + js[r][2] * j[c][2];
                out[c][r] = out[r][c];
            }
        }
        Cov3(out)
    }
}

impl LocalFrame {
    /// An ENU covariance at this observer (m²), as an ECEF covariance (m²). Exact: it's a rotation.
    /// For NED, use [`Cov3::ned_to_enu`] first.
    ///
    /// ```
    /// use geo3d::{Cov3, Geodetic, LocalFrame};
    /// // A GPS fix: 3 m horizontally, 5 m vertically.
    /// let fix = LocalFrame::new(Geodetic::new(-33.87, 151.21, 20.0));
    /// let ecef = fix.enu_cov_to_ecef(Cov3::from_std_devs([3.0, 3.0, 5.0]));
    /// let back = fix.ecef_cov_to_enu(ecef);
    /// assert!((back.std_devs()[2] - 5.0).abs() < 1e-12);
    /// ```
    pub fn enu_cov_to_ecef(&self, cov: Cov3) -> Cov3 {
        let (e, n, u) = self.axes();
        cov.transformed([[e.x, n.x, u.x], [e.y, n.y, u.y], [e.z, n.z, u.z]])
    }

    /// An ECEF covariance (m²) as an ENU covariance at this observer (m²). Exact: it's a rotation.
    pub fn ecef_cov_to_enu(&self, cov: Cov3) -> Cov3 {
        let (e, n, u) = self.axes();
        cov.transformed([[e.x, e.y, e.z], [n.x, n.y, n.z], [u.x, u.y, u.z]])
    }
}

impl Aer {
    /// The ENU covariance (m²) of this pointing, given its AER covariance (azimuth and elevation
    /// in deg², range in m²). First order: see [`Cov3`].
    pub fn cov_to_enu(self, cov: Cov3) -> Cov3 {
        let (saz, caz) = self.azimuth_deg.to_radians().sin_cos();
        let (sel, cel) = self.elevation_deg.to_radians().sin_cos();
        let r = self.range_m;
        // How east, north and up change with azimuth and elevation (per degree) and range.
        cov.transformed([
            [r * cel * caz * DEG, -r * sel * saz * DEG, cel * saz],
            [-r * cel * saz * DEG, -r * sel * caz * DEG, cel * caz],
            [0.0, r * cel * DEG, sel],
        ])
    }
}

impl Enu {
    /// The AER covariance (azimuth and elevation in deg², range in m²) of this offset, given its
    /// ENU covariance (m²). First order: see [`Cov3`].
    ///
    /// Straight up or down, azimuth has no meaning, so the azimuth and elevation variances come
    /// out NaN; the range variance is still right. At zero range everything is NaN.
    pub fn cov_to_aer(self, cov: Cov3) -> Cov3 {
        let (e, n, u) = (self.east, self.north, self.up);
        let h2 = e * e + n * n;
        let (h, r2) = (h2.sqrt(), h2 + u * u);
        let r = r2.sqrt();
        cov.transformed([
            [n / h2 / DEG, -e / h2 / DEG, 0.0],
            [-e * u / (r2 * h) / DEG, -n * u / (r2 * h) / DEG, h / r2 / DEG],
            [e / r, n / r, u / r],
        ])
    }
}

impl Ellipsoid {
    /// The ECEF covariance (m²) of a position `at`, given its geodetic covariance (latitude and
    /// longitude in deg², height in m²). First order: see [`Cov3`].
    pub fn geodetic_cov_to_ecef(self, at: Geodetic, cov: Cov3) -> Cov3 {
        let (e, n, u) = LocalFrame::on(self, at).axes();
        // A degree of latitude moves (M + h)·π/180 metres north; of longitude, (N + h)·cos φ·π/180 east.
        let (dlat, dlon) = self.metres_per_degree(at);
        cov.transformed([
            [n.x * dlat, e.x * dlon, u.x],
            [n.y * dlat, e.y * dlon, u.y],
            [n.z * dlat, e.z * dlon, u.z],
        ])
    }

    /// The geodetic covariance (latitude and longitude in deg², height in m²) of a position `at`,
    /// given its ECEF covariance (m²). First order: see [`Cov3`]. At the poles, longitude has no
    /// meaning, so its variance is enormous.
    ///
    /// ```
    /// use geo3d::{Cov3, Geodetic, WGS84};
    /// let at = Geodetic::new(0.0, 0.0, 0.0);
    /// // 111 m north–south at the equator is about 0.001° of latitude.
    /// let geo = WGS84.ecef_cov_to_geodetic(at, Cov3::from_std_devs([1.0, 1.0, 110.574]));
    /// assert!((geo.std_devs()[0] - 0.001).abs() < 1e-7);
    /// ```
    pub fn ecef_cov_to_geodetic(self, at: Geodetic, cov: Cov3) -> Cov3 {
        let (e, n, u) = LocalFrame::on(self, at).axes();
        let (dlat, dlon) = self.metres_per_degree(at);
        cov.transformed([
            [n.x / dlat, n.y / dlat, n.z / dlat],
            [e.x / dlon, e.y / dlon, e.z / dlon],
            [u.x, u.y, u.z],
        ])
    }

    /// Metres per degree of latitude and of longitude at `at`, along the surface at its height.
    fn metres_per_degree(self, at: Geodetic) -> (f64, f64) {
        let lat = at.lat_rad();
        let h = at.height_m;
        (
            (self.meridional_radius(lat) + h) * DEG,
            (self.prime_vertical_radius(lat) + h) * lat.cos() * DEG,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Ecef, WGS84};

    fn outer(v: [f64; 3]) -> Cov3 {
        Cov3(core::array::from_fn(|i| core::array::from_fn(|j| v[i] * v[j])))
    }

    /// Pushing a rank-one covariance `v·vᵀ` through a conversion `f` must give `(Jv)(Jv)ᵀ`, where
    /// `Jv` is the central difference of `f` along `v`. This checks every Jacobian entry.
    fn check_jacobian(
        propagate: impl Fn(Cov3) -> Cov3,
        f: impl Fn([f64; 3]) -> [f64; 3],
        x: [f64; 3],
        steps: [f64; 3],
    ) {
        for k in 0..3 {
            let mut v = [0.0; 3];
            v[k] = 1.0;
            let (mut lo, mut hi) = (x, x);
            lo[k] -= steps[k];
            hi[k] += steps[k];
            let (a, b) = (f(lo), f(hi));
            let jv: [f64; 3] = core::array::from_fn(|i| (b[i] - a[i]) / (2.0 * steps[k]));
            let got = propagate(outer(v)).0;
            let scale = jv.iter().fold(0.0f64, |m, x| m.max(x.abs()));
            for i in 0..3 {
                for j in 0..3 {
                    let want = jv[i] * jv[j];
                    assert!(
                        (got[i][j] - want).abs() <= 1e-6 * scale * scale,
                        "at {x:?}, column {k}: [{i}][{j}] = {} vs finite difference {want}",
                        got[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn jacobians_match_finite_differences() {
        for &(az, el, r) in &[(40.0, 10.0, 50_000.0), (200.0, -30.0, 1_000.0), (100.0, 80.0, 3.0e7)] {
            let aer = Aer::new(az, el, r);
            check_jacobian(
                |c| aer.cov_to_enu(c),
                |p| {
                    let v = Aer::new(p[0], p[1], p[2]).to_enu();
                    [v.east, v.north, v.up]
                },
                [az, el, r],
                [1e-6, 1e-6, r * 1e-7],
            );
            let enu = aer.to_enu();
            check_jacobian(
                |c| enu.cov_to_aer(c),
                |p| {
                    let a = Enu::new(p[0], p[1], p[2]).to_aer();
                    [a.azimuth_deg, a.elevation_deg, a.range_m]
                },
                [enu.east, enu.north, enu.up],
                [r * 1e-7; 3],
            );
        }
        for &(lat, lon, h) in &[(-33.87, 151.21, 20.0), (0.0, 10.0, -5_000.0), (89.0, -120.0, 4.0e5)] {
            let at = Geodetic::new(lat, lon, h);
            check_jacobian(
                |c| WGS84.geodetic_cov_to_ecef(at, c),
                |p| WGS84.to_ecef(Geodetic::new(p[0], p[1], p[2])).to_array(),
                [lat, lon, h],
                [1e-6, 1e-6, 1.0],
            );
            let x = WGS84.to_ecef(at);
            check_jacobian(
                |c| WGS84.ecef_cov_to_geodetic(at, c),
                |p| {
                    let g = WGS84.to_geodetic(Ecef::new(p[0], p[1], p[2]));
                    [g.lat_deg, g.lon_deg, g.height_m]
                },
                x.to_array(),
                [x.to_array().iter().map(|c| c * c).sum::<f64>().sqrt() * 1e-7; 3],
            );
        }
    }

    /// Normal deviates from a fixed seed (xorshift64* and Box–Muller).
    struct Normal(u64);
    impl Normal {
        fn unit(&mut self) -> f64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            ((self.0.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        }
        fn next(&mut self) -> f64 {
            (-2.0 * self.unit().ln()).sqrt() * (2.0 * core::f64::consts::PI * self.unit()).cos()
        }
    }

    /// The sample covariance of `f` over noisy inputs agrees with the propagated covariance.
    fn monte_carlo(x: [f64; 3], sigma: [f64; 3], f: impl Fn([f64; 3]) -> [f64; 3], propagated: Cov3) {
        let (mut rng, n) = (Normal(0x9e37_79b9_7f4a_7c15), 200_000);
        let samples: Vec<[f64; 3]> = (0..n)
            .map(|_| f(core::array::from_fn(|i| x[i] + sigma[i] * rng.next())))
            .collect();
        let mean: [f64; 3] = core::array::from_fn(|i| samples.iter().map(|s| s[i]).sum::<f64>() / n as f64);
        let p = propagated.0;
        for i in 0..3 {
            for j in 0..3 {
                let s = samples.iter().map(|v| (v[i] - mean[i]) * (v[j] - mean[j])).sum::<f64>() / (n - 1) as f64;
                // Compare as correlations, so every entry has the same scale.
                let d = (s - p[i][j]) / (p[i][i] * p[j][j]).sqrt();
                assert!(d.abs() < 0.01, "[{i}][{j}]: sample {s} vs propagated {}", p[i][j]);
            }
        }
    }

    #[test]
    fn propagated_covariances_match_monte_carlo() {
        // A radar return, and a GPS fix given in latitude, longitude and height.
        let (aer, sigma) = (Aer::new(40.0, 10.0, 50_000.0), [0.2, 0.2, 20.0]);
        monte_carlo(
            [40.0, 10.0, 50_000.0],
            sigma,
            |p| {
                let v = Aer::new(p[0], p[1], p[2]).to_enu();
                [v.east, v.north, v.up]
            },
            aer.cov_to_enu(Cov3::from_std_devs(sigma)),
        );
        let (at, sigma) = (Geodetic::new(-33.87, 151.21, 20.0), [3e-5, 3e-5, 5.0]);
        monte_carlo(
            [-33.87, 151.21, 20.0],
            sigma,
            |p| WGS84.to_ecef(Geodetic::new(p[0], p[1], p[2])).to_array(),
            WGS84.geodetic_cov_to_ecef(at, Cov3::from_std_devs(sigma)),
        );
    }

    /// Equal to rounding, entry by entry on the scale of each entry's standard deviations, and
    /// exactly symmetric.
    fn assert_same(b: Cov3, c: Cov3) {
        for i in 0..3 {
            for j in 0..3 {
                let scale = (c.0[i][i] * c.0[j][j]).sqrt();
                assert!((b.0[i][j] - c.0[i][j]).abs() < 1e-11 * scale, "{b:?} vs {c:?}");
                assert_eq!(b.0[i][j], b.0[j][i]);
            }
        }
    }

    #[test]
    fn round_trips_and_exact_cases() {
        let cov = Cov3([[4.0, 1.0, -0.5], [1.0, 9.0, 2.0], [-0.5, 2.0, 25.0]]);
        let frame = LocalFrame::new(Geodetic::new(51.5, -0.12, 30.0));
        assert_same(frame.ecef_cov_to_enu(frame.enu_cov_to_ecef(cov)), cov);
        let aer = Aer::new(123.0, 45.0, 10_000.0);
        assert_same(aer.to_enu().cov_to_aer(aer.cov_to_enu(cov)), cov);
        assert_same(cov.enu_to_ned().ned_to_enu(), cov);
        // A GPS fix in degrees: about 3 m across, 5 m up, with some correlation.
        let geo = Cov3([[9e-10, 1e-10, 2e-5], [1e-10, 9e-10, -1e-5], [2e-5, -1e-5, 25.0]]);
        let at = Geodetic::new(-33.87, 151.21, 20.0);
        assert_same(WGS84.ecef_cov_to_geodetic(at, WGS84.geodetic_cov_to_ecef(at, geo)), geo);
        // NED swaps north and east, and flips the sign of the up/down cross terms.
        let ned = cov.enu_to_ned().0;
        assert_eq!((ned[0][0], ned[1][1], ned[2][2]), (9.0, 4.0, 25.0));
        assert_eq!((ned[0][1], ned[0][2], ned[1][2]), (1.0, -2.0, 0.5));
        // A rotation keeps the total variance (the trace).
        let ecef = frame.enu_cov_to_ecef(cov).0;
        assert!((ecef[0][0] + ecef[1][1] + ecef[2][2] - 38.0).abs() < 1e-12);
        assert_eq!(Cov3::from_std_devs([3.0, 4.0, 5.0]).std_devs(), [3.0, 4.0, 5.0]);
        // Straight up, azimuth has no meaning.
        let up = Enu::new(0.0, 0.0, 1_000.0).cov_to_aer(Cov3::from_std_devs([1.0; 3]));
        assert!(!up.0[0][0].is_finite());
        assert!((up.std_devs()[2] - 1.0).abs() < 1e-12);
    }
}
