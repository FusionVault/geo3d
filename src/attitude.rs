//! [`Attitude`]: which way a vehicle faces (yaw, pitch and roll), and [`Frd`] offsets in its own
//! body frame (forward, right, down).

use crate::frame::LocalFrame;
use crate::geodesic::{atan2d, sincosd};
use crate::local::{wrap_360, Ned};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;
use crate::vector::{vector_ops, Vec3};

/// A 3×3 rotation, row by row: `out = M · in` for column vectors.
type Rotation = [[f64; 3]; 3];

/// An offset in a vehicle's own **body frame**, metres: forward (out of the nose), right (out of
/// the right wing) and down (out of the belly). The aerospace "FRD" convention. Supports `+`, `-`,
/// `* f64` and unary `-`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Frd {
    /// Forward, metres.
    pub forward: f64,
    /// Right, metres.
    pub right: f64,
    /// Down, metres.
    pub down: f64,
}

impl Frd {
    /// Build from components in metres.
    #[inline]
    pub const fn new(forward: f64, right: f64, down: f64) -> Self {
        Frd { forward, right, down }
    }

    /// Length, metres.
    #[inline]
    pub fn norm(self) -> f64 {
        self.vec().norm()
    }

    /// As a plain vector `(forward, right, down)`.
    #[inline]
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.forward, self.right, self.down)
    }
}

vector_ops!(Frd { forward, right, down });

/// Which way a vehicle faces: **yaw** (its heading, degrees clockwise from north), **pitch**
/// (degrees nose up) and **roll** (degrees right wing down), applied in that order. This is the
/// aerospace convention, measured from North-East-Down.
///
/// It turns offsets in the vehicle's body frame ([`Frd`]) into NED and back. With a
/// [`LocalFrame`] at the vehicle, it also turns them into positions and ECEF.
///
/// ```
/// use geo3d::{Attitude, Frd};
/// // Heading east, nose 10° up.
/// let plane = Attitude::new(90.0, 10.0, 0.0);
/// let ahead = plane.body_to_ned(Frd::new(100.0, 0.0, 0.0));   // 100 m out of the nose
/// assert!((ahead.east - 98.48).abs() < 0.01 && (ahead.down + 17.36).abs() < 0.01); // east, 17 m up
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Attitude {
    /// Heading, degrees clockwise from north.
    pub yaw_deg: f64,
    /// Nose up, degrees.
    pub pitch_deg: f64,
    /// Right wing down, degrees.
    pub roll_deg: f64,
}

impl Attitude {
    /// Build from yaw, pitch and roll in degrees.
    #[inline]
    pub const fn new(yaw_deg: f64, pitch_deg: f64, roll_deg: f64) -> Self {
        Attitude {
            yaw_deg,
            pitch_deg,
            roll_deg,
        }
    }

    /// The rotation from the body frame to NED, row by row: `ned = M · frd`. Its transpose turns
    /// NED back into the body frame. Exact at whole multiples of 90°.
    pub fn to_rotation(self) -> [[f64; 3]; 3] {
        let (sy, cy) = sincosd(self.yaw_deg);
        let (sp, cp) = sincosd(self.pitch_deg);
        let (sr, cr) = sincosd(self.roll_deg);
        [
            [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
            [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
            [-sp, cp * sr, cp * cr],
        ]
    }

    /// The attitude of a body-to-NED rotation matrix (row by row, as from
    /// [`to_rotation`](Self::to_rotation)). Yaw comes out in `[0, 360)`, pitch in `[-90, 90]` and
    /// roll in `(-180, 180]`.
    ///
    /// It stays precise right up to ±90° pitch. At exactly ±90°, yaw and roll turn about the same
    /// axis, so only their combination means anything: roll comes out 0 and yaw takes the whole
    /// turn. Either way, the angles rebuild the same rotation to within rounding. A matrix with a NaN
    /// or infinite entry gives NaN angles.
    ///
    /// ```
    /// use geo3d::Attitude;
    /// let up = Attitude::new(30.0, 89.999_999_9, 10.0);
    /// let back = Attitude::from_rotation(up.to_rotation());
    /// assert!((back.pitch_deg - 89.999_999_9).abs() < 1e-12);
    /// ```
    pub fn from_rotation(m: [[f64; 3]; 3]) -> Attitude {
        // Not every entry is read below, so check them all.
        if !m.iter().flatten().all(|x| x.is_finite()) {
            return Attitude::new(f64::NAN, f64::NAN, f64::NAN);
        }
        // Roll from the bottom row, then yaw from the full-size entries given that roll, so the
        // three angles always rebuild `m`, even where yaw and roll can't be told apart.
        let roll = if m[2][1] == 0.0 && m[2][2] == 0.0 {
            0.0
        } else {
            atan2d(m[2][1], m[2][2])
        };
        let (sr, cr) = sincosd(roll);
        let yaw = atan2d(sr * m[0][2] - cr * m[0][1], cr * m[1][1] - sr * m[1][2]);
        // atan2 rather than asin(−m[2][0]), which loses precision near ±90° pitch.
        let pitch = atan2d(-m[2][0], m[2][1].hypot(m[2][2]));
        Attitude::new(wrap_360(yaw), pitch, roll)
    }

    /// An offset in the body frame as North-East-Down.
    pub fn body_to_ned(self, v: Frd) -> Ned {
        let [n, e, d] = apply(self.to_rotation(), v.into());
        Ned::new(n, e, d)
    }

    /// A North-East-Down offset in the body frame.
    pub fn ned_to_body(self, v: Ned) -> Frd {
        Frd::from(apply(transpose(self.to_rotation()), [v.north, v.east, v.down]))
    }
}

impl LocalFrame {
    /// The rotation from this frame's ENU axes to ECEF, row by row: `ecef = M · enu` for a
    /// direction. Its columns are the east, north and up unit vectors; its transpose turns ECEF
    /// into ENU.
    pub fn enu_to_ecef_rotation(&self) -> [[f64; 3]; 3] {
        let (e, n, u) = self.axes();
        [[e.x, n.x, u.x], [e.y, n.y, u.y], [e.z, n.z, u.z]]
    }

    /// The rotation from this frame's NED axes to ECEF, row by row: `ecef = M · ned` for a
    /// direction. Its transpose turns ECEF into NED.
    pub fn ned_to_ecef_rotation(&self) -> [[f64; 3]; 3] {
        let (e, n, u) = self.axes();
        [[n.x, e.x, -u.x], [n.y, e.y, -u.y], [n.z, e.z, -u.z]]
    }

    /// The rotation from the body frame of a vehicle at this frame's origin, with `attitude`, to
    /// ECEF (row by row). Cesium calls this `headingPitchRollToFixedFrame`.
    pub fn body_to_ecef_rotation(&self, attitude: Attitude) -> [[f64; 3]; 3] {
        mul(self.ned_to_ecef_rotation(), attitude.to_rotation())
    }

    /// The yaw, pitch and roll, seen from this frame, of a body-to-ECEF rotation: the reverse of
    /// [`body_to_ecef_rotation`](Self::body_to_ecef_rotation). Use it to read an attitude given in
    /// ECEF, say from an inertial navigation system, as an ordinary heading, pitch and roll.
    pub fn attitude_from_ecef_rotation(&self, m: [[f64; 3]; 3]) -> Attitude {
        Attitude::from_rotation(mul(transpose(self.ned_to_ecef_rotation()), m))
    }

    /// The position at an offset in the body frame of a vehicle at this frame's origin, with
    /// `attitude`: where a sensor, a lever arm or a target seen from the vehicle is.
    ///
    /// ```
    /// use geo3d::{Attitude, Frd, Geodetic, LocalFrame};
    /// let drone = LocalFrame::new(Geodetic::new(-33.94, 151.18, 120.0));
    /// let heading = Attitude::new(90.0, 0.0, 0.0);            // facing east, level
    /// let camera = drone.geodetic_of_body(heading, Frd::new(0.2, 0.0, 0.1));
    /// assert!(camera.lon_deg > 151.18 && camera.height_m < 120.0);
    /// ```
    pub fn geodetic_of_body(&self, attitude: Attitude, offset: Frd) -> Geodetic {
        self.geodetic_of_ned(attitude.body_to_ned(offset))
    }

    /// Where a position is in the body frame of a vehicle at this frame's origin, with
    /// `attitude`: how far forward, right and down of it.
    pub fn body_of(&self, attitude: Attitude, target: Geodetic) -> Frd {
        attitude.ned_to_body(self.ned_of(target))
    }
}

fn apply(m: Rotation, v: [f64; 3]) -> [f64; 3] {
    m.map(|row| row[0] * v[0] + row[1] * v[1] + row[2] * v[2])
}

fn transpose(m: Rotation) -> Rotation {
    core::array::from_fn(|r| core::array::from_fn(|c| m[c][r]))
}

fn mul(a: Rotation, b: Rotation) -> Rotation {
    core::array::from_fn(|r| core::array::from_fn(|c| a[r][0] * b[0][c] + a[r][1] * b[1][c] + a[r][2] * b[2][c]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::WGS84;

    fn max_diff(a: Rotation, b: Rotation) -> f64 {
        (0..9)
            .map(|i| (a[i / 3][i % 3] - b[i / 3][i % 3]).abs())
            .fold(0.0, f64::max)
    }

    #[test]
    fn the_axes_turn_the_right_way() {
        let nose = Frd::new(1.0, 0.0, 0.0);
        let wing = Frd::new(0.0, 1.0, 0.0);
        // Yaw 90°: the nose points east. Pitch 90°: the nose points up. Roll 90°: the right wing
        // points down. All exact.
        assert_eq!(Attitude::new(90.0, 0.0, 0.0).body_to_ned(nose), Ned::new(0.0, 1.0, 0.0));
        assert_eq!(
            Attitude::new(0.0, 90.0, 0.0).body_to_ned(nose),
            Ned::new(0.0, 0.0, -1.0)
        );
        assert_eq!(Attitude::new(0.0, 0.0, 90.0).body_to_ned(wing), Ned::new(0.0, 0.0, 1.0));
        // Heading south, the right wing points west.
        assert_eq!(
            Attitude::new(180.0, 0.0, 0.0).body_to_ned(wing),
            Ned::new(0.0, -1.0, 0.0)
        );
        // Applied yaw, then pitch, then roll: heading east with the nose up 30°, then rolled 90°
        // right, the right wing takes the place of "down": east and down, at right angles to the
        // nose (east and up).
        let a = Attitude::new(90.0, 30.0, 90.0);
        let (n, w) = (a.body_to_ned(nose), a.body_to_ned(wing));
        assert!((n.vec().dot(w.vec())).abs() < 1e-15);
        let want = Ned::new(0.0, 0.5, 3f64.sqrt() / 2.0);
        assert!((w - want).norm() < 1e-15, "{w:?}");
    }

    #[test]
    fn rotations_are_orthonormal_and_invert() {
        let a = Attitude::new(123.4, -56.7, 89.0);
        let m = a.to_rotation();
        let i = mul(m, transpose(m));
        assert!(max_diff(i, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) < 4e-16);
        let v = Frd::new(3.0, -4.0, 12.0);
        let back = a.ned_to_body(a.body_to_ned(v));
        assert!((back - v).norm() < 1e-14);
    }

    /// Angles come back where they're well defined, and the rotation always does, right up to and
    /// at ±90° pitch.
    #[test]
    fn from_rotation_round_trips_up_to_the_vertical() {
        let mut s: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut u = || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            (s.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 / (1u64 << 53) as f64
        };
        for k in 0..20_000 {
            let (yaw, roll) = (u() * 360.0, u() * 360.0 - 180.0);
            // Spread over every pitch, with half of them within 1e-15° to 1e-1° of ±90°.
            let pitch = if k % 2 == 0 {
                u() * 180.0 - 90.0
            } else {
                let near = 90.0 - 10f64.powf(-u() * 15.0);
                if k % 4 == 1 {
                    near
                } else {
                    -near
                }
            };
            let a = Attitude::new(yaw, pitch, roll);
            let m = a.to_rotation();
            let b = Attitude::from_rotation(m);
            assert!(max_diff(b.to_rotation(), m) < 2e-15, "{a:?} -> {b:?}");
            assert!((b.pitch_deg - pitch).abs() < 1e-13, "{a:?} -> {b:?}");
            if pitch.abs() < 89.9 {
                assert!(
                    (b.yaw_deg - yaw).abs() < 1e-11 && (b.roll_deg - roll).abs() < 1e-11,
                    "{a:?} -> {b:?}"
                );
            }
        }
        // Exactly vertical: roll is 0 and yaw takes the whole turn (yaw − roll nose up, yaw + roll
        // nose down).
        let up = Attitude::from_rotation(Attitude::new(30.0, 90.0, 10.0).to_rotation());
        assert_eq!((up.pitch_deg, up.roll_deg), (90.0, 0.0));
        assert!((up.yaw_deg - 20.0).abs() < 1e-12, "{up:?}");
        let down = Attitude::from_rotation(Attitude::new(30.0, -90.0, 10.0).to_rotation());
        assert_eq!((down.pitch_deg, down.roll_deg), (-90.0, 0.0));
        assert!((down.yaw_deg - 40.0).abs() < 1e-12, "{down:?}");
        // Ranges.
        let b = Attitude::from_rotation(Attitude::new(-10.0, 0.0, 190.0).to_rotation());
        assert!(
            (b.yaw_deg - 350.0).abs() < 1e-12 && (b.roll_deg + 170.0).abs() < 1e-12,
            "{b:?}"
        );
    }

    #[test]
    fn ecef_rotations_agree_with_the_frame() {
        let frame = LocalFrame::new(Geodetic::new(-33.94, 151.18, 21.0));
        let a = Attitude::new(250.0, 12.0, -7.0);
        let v = Frd::new(10.0, -3.0, 2.0);
        // Body → NED → ENU → ECEF by hand, and in one rotation.
        let by_hand = frame.dir_to_ecef(a.body_to_ned(v).to_enu());
        let m = frame.body_to_ecef_rotation(a);
        let r = apply(m, v.into());
        assert!((Vec3::new(r[0], r[1], r[2]) - by_hand).norm() < 1e-13);
        // And the attitude comes back out of the ECEF rotation.
        let b = frame.attitude_from_ecef_rotation(m);
        assert!(
            (b.yaw_deg - 250.0).abs() < 1e-11 && (b.pitch_deg - 12.0).abs() < 1e-11 && (b.roll_deg + 7.0).abs() < 1e-11
        );
        // ENU and NED rotations are the frame's axes.
        let enu = frame.enu_to_ecef_rotation();
        let up = frame.dir_to_ecef(crate::Enu::new(0.0, 0.0, 1.0));
        assert_eq!([enu[0][2], enu[1][2], enu[2][2]], [up.x, up.y, up.z]);
        let ned = frame.ned_to_ecef_rotation();
        assert_eq!([ned[0][2], ned[1][2], ned[2][2]], [-up.x, -up.y, -up.z]);
    }

    #[test]
    fn body_offsets_and_back() {
        let frame = LocalFrame::on(WGS84, Geodetic::new(51.5, -0.12, 100.0));
        let a = Attitude::new(45.0, 5.0, -3.0);
        let target = Geodetic::new(51.51, -0.11, 300.0);
        let frd = frame.body_of(a, target);
        let back = frame.geodetic_of_body(a, frd);
        assert!(crate::distance3d_m(back, target) < 1e-6, "{back:?}");
        // Straight ahead and level, about 1.3 km away to the north-east.
        assert!(frd.forward > 1_000.0 && frd.right.abs() < 300.0, "{frd:?}");
    }
}
