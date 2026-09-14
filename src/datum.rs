//! Datum shifts.

use crate::types::{fma, Ecef, Vec3};

/// A 7-parameter Helmert (Bursa-Wolf) transformation between ECEF frames: translations (metres),
/// rotations (**radians**) and a scale factor (unitless, e.g. `ppm × 1e-6`).
///
/// Rotation sign follows the **position-vector** convention (EPSG method 9606, the one most European
/// datums publish). For parameters published in the coordinate-frame convention (EPSG 9607, common in
/// North America and Australia) use [`Helmert7::coordinate_frame`], which flips the rotation signs.
/// Parameter values are per datum pair — look them up in EPSG; none are hardcoded here, because a
/// wrong constant is worse than none.
///
/// ```
/// use geo3d::{Ecef, Helmert7};
/// // A published position-vector set: tx ty tz (m), rx ry rz (arc-seconds), ds (ppm).
/// let h = Helmert7::position_vector(-446.448, 125.157, -542.060, -0.1502, -0.2470, -0.8421, 20.4894);
/// let p = Ecef::new(3_874_938.0, 116_218.0, 5_047_168.0);
/// let back = h.apply_inverse(h.apply(p));
/// assert!((back.x - p.x).abs() < 1e-6 && (back.y - p.y).abs() < 1e-6 && (back.z - p.z).abs() < 1e-6);
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Helmert7 {
    /// Translation along X, metres.
    pub tx: f64,
    /// Translation along Y, metres.
    pub ty: f64,
    /// Translation along Z, metres.
    pub tz: f64,
    /// Rotation about X, radians (position-vector convention).
    pub rx: f64,
    /// Rotation about Y, radians (position-vector convention).
    pub ry: f64,
    /// Rotation about Z, radians (position-vector convention).
    pub rz: f64,
    /// Scale difference, unitless (`ppm × 1e-6`).
    pub scale: f64,
}

const ARCSEC: f64 = core::f64::consts::PI / (180.0 * 3600.0);

impl Helmert7 {
    /// Raw parameters: metres, radians, unitless scale.
    #[inline]
    pub const fn new(tx: f64, ty: f64, tz: f64, rx: f64, ry: f64, rz: f64, scale: f64) -> Self {
        Helmert7 {
            tx,
            ty,
            tz,
            rx,
            ry,
            rz,
            scale,
        }
    }

    /// From parameters as EPSG publishes them in the **position-vector** convention: translations in
    /// metres, rotations in arc-seconds, scale difference in parts per million.
    #[inline]
    pub fn position_vector(
        tx: f64,
        ty: f64,
        tz: f64,
        rx_arcsec: f64,
        ry_arcsec: f64,
        rz_arcsec: f64,
        ds_ppm: f64,
    ) -> Self {
        Helmert7::new(
            tx,
            ty,
            tz,
            rx_arcsec * ARCSEC,
            ry_arcsec * ARCSEC,
            rz_arcsec * ARCSEC,
            ds_ppm * 1e-6,
        )
    }

    /// From parameters published in the **coordinate-frame** convention (rotation signs flipped).
    #[inline]
    pub fn coordinate_frame(
        tx: f64,
        ty: f64,
        tz: f64,
        rx_arcsec: f64,
        ry_arcsec: f64,
        rz_arcsec: f64,
        ds_ppm: f64,
    ) -> Self {
        Helmert7::position_vector(tx, ty, tz, -rx_arcsec, -ry_arcsec, -rz_arcsec, ds_ppm)
    }

    /// Apply the transformation to an ECEF point (source frame → target frame). Small-angle form,
    /// evaluated with FMA so each rotated component carries one rounding instead of several.
    #[inline]
    pub fn apply(&self, p: Ecef) -> Ecef {
        let m = 1.0 + self.scale;
        Ecef::new(
            fma(m, fma(self.ry, p.z, fma(self.rz, -p.y, p.x)), self.tx),
            fma(m, fma(self.rx, -p.z, fma(self.rz, p.x, p.y)), self.ty),
            fma(m, fma(self.rx, p.y, fma(self.ry, -p.x, p.z)), self.tz),
        )
    }

    /// Undo the transformation exactly (target frame → source frame): solves the linear form for
    /// the source point, so `apply_inverse(apply(p))` round-trips to floating-point precision.
    pub fn apply_inverse(&self, p: Ecef) -> Ecef {
        let m = 1.0 + self.scale;
        // q = R·s where R = [[1, -rz, ry], [rz, 1, -rx], [-ry, rx, 1]]; solve by Cramer's rule.
        let q = Vec3::new((p.x - self.tx) / m, (p.y - self.ty) / m, (p.z - self.tz) / m);
        let (rx, ry, rz) = (self.rx, self.ry, self.rz);
        let det = 1.0 + rx * rx + ry * ry + rz * rz;
        let x = (q.x * (1.0 + rx * rx) + q.y * (rz + rx * ry) + q.z * (rx * rz - ry)) / det;
        let y = (q.x * (rx * ry - rz) + q.y * (1.0 + ry * ry) + q.z * (rx + ry * rz)) / det;
        let z = (q.x * (rx * rz + ry) + q.y * (ry * rz - rx) + q.z * (1.0 + rz * rz)) / det;
        Ecef::new(x, y, z)
    }

    /// The **first-order** inverse: every parameter negated — the form in which reverse parameter
    /// sets are conventionally published. Good to about a centimetre for typical (ppm, arc-second)
    /// sets; use [`apply_inverse`](Self::apply_inverse) when you need the exact reverse.
    #[inline]
    pub fn inverse(&self) -> Self {
        Helmert7::new(-self.tx, -self.ty, -self.tz, -self.rx, -self.ry, -self.rz, -self.scale)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_and_translation() {
        let p = Ecef::new(1000.0, 2000.0, 3000.0);
        assert_eq!(Helmert7::default().apply(p), p);
        let t = Helmert7::new(10.0, -5.0, 2.0, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(t.apply(p), Ecef::new(1010.0, 1995.0, 3002.0));
    }

    #[test]
    fn exact_inverse_round_trips_and_first_order_is_close() {
        let h = Helmert7::position_vector(-446.448, 125.157, -542.060, -0.1502, -0.2470, -0.8421, 20.4894);
        let p = Ecef::new(3_874_938.0, 116_218.0, 5_047_168.0);
        let exact = h.apply_inverse(h.apply(p));
        assert!(exact.distance_to(p) < 1e-6, "{exact:?}");
        let approx = h.inverse().apply(h.apply(p));
        assert!(approx.distance_to(p) < 0.02, "{approx:?}");
    }

    #[test]
    fn conventions_differ_by_rotation_sign() {
        let a = Helmert7::position_vector(0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0);
        let b = Helmert7::coordinate_frame(0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 0.0);
        assert_eq!(a.rx, -b.rx);
        assert!((a.rz - 3.0 * ARCSEC).abs() < 1e-18);
    }
}
