//! [`Quaternion`]: a rotation as four numbers, the way flight controllers store attitude.

use crate::angle::sincosd;
use crate::attitude::Attitude;
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;

/// A rotation as a unit quaternion `w + xi + yj + zk`: the scalar part `w`, then the vector part
/// `x, y, z`.
///
/// It uses the Hamilton convention, scalar first, as PX4 and MAVLink's `ATTITUDE_QUATERNION` do.
/// For an attitude it turns body (FRD) vectors into NED, like [`Attitude::to_rotation`]. A
/// quaternion and its negative are the same rotation; geo3d always returns the one with `w > 0`
/// (or, for a half turn with `w = 0`, the first non-zero of `x`, `y`, `z` positive).
///
/// Combine rotations with `*`: `(a * b)` turns by `b`, then by `a`, matching the rotation matrices
/// `A · B`.
///
/// ```
/// use geo3d::{Attitude, Quaternion};
/// // Heading 90°: a quarter turn about down.
/// let q = Attitude::new(90.0, 0.0, 0.0).to_quaternion();
/// let half = std::f64::consts::FRAC_1_SQRT_2;
/// assert!((q.w - half).abs() < 2e-16 && (q.z - half).abs() < 2e-16);
/// assert!((Attitude::from_quaternion(q).yaw_deg - 90.0).abs() < 1e-12);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Quaternion {
    /// The scalar part, `cos(angle / 2)`.
    pub w: f64,
    /// The vector part along x, `sin(angle / 2)` times the axis.
    pub x: f64,
    /// The vector part along y.
    pub y: f64,
    /// The vector part along z.
    pub z: f64,
}

impl Default for Quaternion {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Quaternion {
    /// No rotation.
    pub const IDENTITY: Quaternion = Quaternion::new(1.0, 0.0, 0.0, 0.0);

    /// Build from the scalar part `w` and the vector part `x, y, z`.
    #[inline]
    pub const fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Quaternion { w, x, y, z }
    }

    /// The length; 1 for a rotation.
    #[inline]
    pub fn norm(self) -> f64 {
        (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Scaled to length 1, with the sign chosen as described above. A zero quaternion gives NaN.
    pub fn normalized(self) -> Quaternion {
        let s = self.sign() / self.norm();
        Quaternion::new(self.w * s, self.x * s, self.y * s, self.z * s)
    }

    /// −1 if the negative is the canonical one (see above), otherwise 1.
    fn sign(self) -> f64 {
        let lead = [self.w, self.x, self.y, self.z]
            .into_iter()
            .find(|&c| c != 0.0)
            .unwrap_or(0.0);
        if lead < 0.0 {
            -1.0
        } else {
            1.0
        }
    }

    /// The reverse rotation (for a unit quaternion).
    #[inline]
    pub const fn conjugate(self) -> Quaternion {
        Quaternion::new(self.w, -self.x, -self.y, -self.z)
    }

    /// The same rotation as a 3×3 matrix, row by row: `out = M · in`. Works for any non-zero
    /// quaternion, not just unit ones.
    pub fn to_rotation(self) -> [[f64; 3]; 3] {
        let Quaternion { w, x, y, z } = self;
        let s = 2.0 / (w * w + x * x + y * y + z * z);
        let (xx, yy, zz) = (s * x * x, s * y * y, s * z * z);
        let (xy, xz, yz) = (s * x * y, s * x * z, s * y * z);
        let (wx, wy, wz) = (s * w * x, s * w * y, s * w * z);
        [
            [1.0 - (yy + zz), xy - wz, xz + wy],
            [xy + wz, 1.0 - (xx + zz), yz - wx],
            [xz - wy, yz + wx, 1.0 - (xx + yy)],
        ]
    }

    /// The quaternion of a rotation matrix (row by row). Uses Shepperd's method: it starts from
    /// the largest of the four components, so it stays precise for every rotation, including half
    /// turns where `w` is close to 0.
    pub fn from_rotation(m: [[f64; 3]; 3]) -> Quaternion {
        let trace = m[0][0] + m[1][1] + m[2][2];
        let q = if trace >= m[0][0] && trace >= m[1][1] && trace >= m[2][2] {
            let w = (1.0 + trace).sqrt() * 2.0; // 4w
            Quaternion::new(
                w / 4.0,
                (m[2][1] - m[1][2]) / w,
                (m[0][2] - m[2][0]) / w,
                (m[1][0] - m[0][1]) / w,
            )
        } else if m[0][0] >= m[1][1] && m[0][0] >= m[2][2] {
            let x = (1.0 + m[0][0] - m[1][1] - m[2][2]).sqrt() * 2.0; // 4x
            Quaternion::new(
                (m[2][1] - m[1][2]) / x,
                x / 4.0,
                (m[0][1] + m[1][0]) / x,
                (m[0][2] + m[2][0]) / x,
            )
        } else if m[1][1] >= m[2][2] {
            let y = (1.0 - m[0][0] + m[1][1] - m[2][2]).sqrt() * 2.0; // 4y
            Quaternion::new(
                (m[0][2] - m[2][0]) / y,
                (m[0][1] + m[1][0]) / y,
                y / 4.0,
                (m[1][2] + m[2][1]) / y,
            )
        } else {
            let z = (1.0 - m[0][0] - m[1][1] + m[2][2]).sqrt() * 2.0; // 4z
            Quaternion::new(
                (m[1][0] - m[0][1]) / z,
                (m[0][2] + m[2][0]) / z,
                (m[1][2] + m[2][1]) / z,
                z / 4.0,
            )
        };
        q.normalized()
    }
}

/// The Hamilton product: `a * b` turns by `b`, then by `a`.
impl core::ops::Mul for Quaternion {
    type Output = Quaternion;
    #[inline]
    fn mul(self, b: Quaternion) -> Quaternion {
        let a = self;
        Quaternion::new(
            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        )
    }
}

impl Attitude {
    /// This attitude as a unit quaternion turning body (FRD) vectors into NED, with the canonical
    /// sign (see [`Quaternion`]). Built
    /// from the half angles directly, so it's exact at whole multiples of 180°.
    pub fn to_quaternion(self) -> Quaternion {
        let (sy, cy) = sincosd(self.yaw_deg / 2.0);
        let (sp, cp) = sincosd(self.pitch_deg / 2.0);
        let (sr, cr) = sincosd(self.roll_deg / 2.0);
        let q = Quaternion::new(
            cy * cp * cr + sy * sp * sr,
            cy * cp * sr - sy * sp * cr,
            cy * sp * cr + sy * cp * sr,
            sy * cp * cr - cy * sp * sr,
        );
        // Already unit length to rounding; only the sign is chosen.
        let s = q.sign();
        Quaternion::new(s * q.w, s * q.x, s * q.y, s * q.z)
    }

    /// The yaw, pitch and roll of a quaternion that turns body (FRD) vectors into NED. It needn't
    /// be exactly unit length. Exact pointing straight up or down, as
    /// [`from_rotation`](Self::from_rotation) is.
    pub fn from_quaternion(q: Quaternion) -> Attitude {
        Attitude::from_rotation(q.to_rotation())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(seed: u64) -> impl FnMut() -> f64 {
        let mut s = seed;
        move || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            (s.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn max_diff(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
        (0..9)
            .map(|i| (a[i / 3][i % 3] - b[i / 3][i % 3]).abs())
            .fold(0.0, f64::max)
    }

    fn mat_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        core::array::from_fn(|r| core::array::from_fn(|c| (0..3).map(|k| a[r][k] * b[k][c]).sum()))
    }

    fn close(a: Quaternion, b: Quaternion, tol: f64) -> bool {
        (a.w - b.w)
            .abs()
            .max((a.x - b.x).abs())
            .max((a.y - b.y).abs())
            .max((a.z - b.z).abs())
            <= tol
    }

    #[test]
    fn known_quaternions() {
        let h = core::f64::consts::FRAC_1_SQRT_2;
        // Quarter turns about down (yaw), right (pitch) and forward (roll). (sin 45° and cos 45° are
        // a rounding step apart in floating point.)
        assert!(close(
            Attitude::new(90.0, 0.0, 0.0).to_quaternion(),
            Quaternion::new(h, 0.0, 0.0, h),
            2e-16
        ));
        assert!(close(
            Attitude::new(0.0, 90.0, 0.0).to_quaternion(),
            Quaternion::new(h, 0.0, h, 0.0),
            2e-16
        ));
        assert!(close(
            Attitude::new(0.0, 0.0, 90.0).to_quaternion(),
            Quaternion::new(h, h, 0.0, 0.0),
            2e-16
        ));
        // A half turn in yaw is exact: w = 0, and the sign convention keeps w ≥ 0.
        assert_eq!(
            Attitude::new(180.0, 0.0, 0.0).to_quaternion(),
            Quaternion::new(0.0, 0.0, 0.0, 1.0)
        );
        assert_eq!(
            Attitude::new(-180.0, 0.0, 0.0).to_quaternion(),
            Quaternion::new(0.0, 0.0, 0.0, 1.0)
        );
        assert_eq!(Attitude::new(0.0, 0.0, 0.0).to_quaternion(), Quaternion::IDENTITY);
        assert_eq!(Quaternion::default(), Quaternion::IDENTITY);
        assert_eq!(
            Quaternion::from_rotation(Quaternion::IDENTITY.to_rotation()),
            Quaternion::IDENTITY
        );
    }

    #[test]
    fn matches_the_rotation_matrix_everywhere() {
        let mut u = rng(0x9e37_79b9_7f4a_7c15);
        for k in 0..20_000 {
            let pitch = if k % 2 == 0 {
                u() * 180.0 - 90.0
            } else {
                (90.0 - 10f64.powf(-u() * 15.0)) * if k % 4 == 1 { 1.0 } else { -1.0 }
            };
            let a = Attitude::new(u() * 360.0, pitch, u() * 360.0 - 180.0);
            let q = a.to_quaternion();
            assert!((q.norm() - 1.0).abs() < 4e-16 && q.sign() > 0.0, "{q:?}");
            // The quaternion and the angles give the same rotation, and each converts back.
            let m = a.to_rotation();
            assert!(max_diff(q.to_rotation(), m) < 2e-15, "{a:?}");
            let back = Quaternion::from_rotation(m);
            assert!(
                (back.w - q.w)
                    .abs()
                    .max((back.x - q.x).abs())
                    .max((back.y - q.y).abs())
                    .max((back.z - q.z).abs())
                    < 2e-15,
                "{a:?}: {back:?} vs {q:?}"
            );
            let b = Attitude::from_quaternion(q);
            assert!(max_diff(b.to_rotation(), m) < 4e-15, "{a:?} -> {b:?}");
            assert!((b.pitch_deg - pitch).abs() < 1e-13, "{a:?} -> {b:?}");
        }
    }

    #[test]
    fn half_turns_and_any_axis() {
        // Shepperd's method on rotations near a half turn (w near 0), about random axes.
        let mut u = rng(42);
        for _ in 0..20_000 {
            let (ax, ay, az) = (u() - 0.5, u() - 0.5, u() - 0.5);
            let len = (ax * ax + ay * ay + az * az).sqrt();
            let angle = core::f64::consts::PI - 10f64.powf(-u() * 12.0) * if u() < 0.5 { 1.0 } else { -1.0 };
            let (s, c) = (angle / 2.0).sin_cos();
            let q = Quaternion::new(c, s * ax / len, s * ay / len, s * az / len).normalized();
            let back = Quaternion::from_rotation(q.to_rotation());
            let d = (back.w - q.w)
                .abs()
                .max((back.x - q.x).abs())
                .max((back.y - q.y).abs())
                .max((back.z - q.z).abs());
            assert!(d < 4e-15, "{q:?} -> {back:?}");
        }
    }

    #[test]
    fn products_compose_like_matrices() {
        let a = Attitude::new(30.0, 10.0, -5.0).to_quaternion();
        let b = Attitude::new(200.0, -40.0, 70.0).to_quaternion();
        assert!(max_diff((a * b).to_rotation(), mat_mul(a.to_rotation(), b.to_rotation())) < 2e-15);
        let i = a * a.conjugate();
        assert!(close(i, Quaternion::IDENTITY, 1e-15), "{i:?}");
        // Not unit length: the rotation is the same, and normalizing fixes the length and sign.
        let big = Quaternion::new(-3.0 * a.w, -3.0 * a.x, -3.0 * a.y, -3.0 * a.z);
        assert!(max_diff(big.to_rotation(), a.to_rotation()) < 2e-15);
        let n = big.normalized();
        assert!(close(n, a, 1e-15), "{n:?}");
        assert!(Quaternion::new(0.0, 0.0, 0.0, 0.0).normalized().w.is_nan());
    }
}
