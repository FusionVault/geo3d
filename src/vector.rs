//! [`Vec3`], the plain 3D vector the other types are built on.

#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use core::ops::{Add, Div, Mul, Sub};

use crate::angle::two_sum;

/// Fused multiply-add `a·b + c`, but only where the target has an FMA instruction. Hardware FMA
/// rounds once and is a single instruction; without it `f64::mul_add` becomes a slow software call,
/// several times slower than `a*b + c`. So stock targets use the plain form (the crate's results are
/// unchanged), and `-C target-cpu=native` or `-C target-feature=+fma` builds get the fused path.
#[inline(always)]
pub(crate) fn fma(a: f64, b: f64, c: f64) -> f64 {
    if cfg!(target_feature = "fma") {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

/// The exact product `a·b` as `hi + lo`: one fused multiply–add where the CPU has it, Dekker's
/// split otherwise.
#[inline]
pub(crate) fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    if cfg!(target_feature = "fma") {
        (p, a.mul_add(b, -p))
    } else {
        let split = |x: f64| {
            let c = 134_217_729.0 * x; // 2²⁷ + 1
            let hi = c - (c - x);
            (hi, x - hi)
        };
        let ((ah, al), (bh, bl)) = (split(a), split(b));
        (p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
    }
}

/// A double-double number, `hi + lo`: about 32 significant digits, for sums that would otherwise
/// cancel.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Dd(pub f64, pub f64);

impl Dd {
    /// `a·b`, exactly.
    #[inline]
    pub fn prod(a: f64, b: f64) -> Dd {
        let (hi, lo) = two_prod(a, b);
        Dd(hi, lo)
    }

    /// `x·y`, exactly (three products, summed in double-double).
    #[inline]
    pub fn dot(x: Vec3, y: Vec3) -> Dd {
        Dd::prod(x.x, y.x) + Dd::prod(x.y, y.y) + Dd::prod(x.z, y.z)
    }

    /// Rounded to the nearest `f64`.
    #[inline]
    pub fn f64(self) -> f64 {
        self.0 + self.1
    }
}

impl Add for Dd {
    type Output = Dd;
    #[inline]
    fn add(self, o: Dd) -> Dd {
        let (s, e) = two_sum(self.0, o.0);
        let e = e + self.1 + o.1;
        let hi = s + e;
        Dd(hi, e - (hi - s))
    }
}

impl Sub for Dd {
    type Output = Dd;
    #[inline]
    fn sub(self, o: Dd) -> Dd {
        self + Dd(-o.0, -o.1)
    }
}

impl Mul for Dd {
    type Output = Dd;
    #[inline]
    fn mul(self, o: Dd) -> Dd {
        let (p, e) = two_prod(self.0, o.0);
        let e = e + self.0 * o.1 + self.1 * o.0;
        let hi = p + e;
        Dd(hi, e - (hi - p))
    }
}

/// `+`, `-`, `* f64`, unary `-` and `[f64; 3]` conversions for a three-field vector type.
macro_rules! vector_ops {
    ($t:ident { $a:ident, $b:ident, $c:ident }) => {
        impl core::ops::Add for $t {
            type Output = $t;
            #[inline]
            fn add(self, o: $t) -> $t {
                $t::new(self.$a + o.$a, self.$b + o.$b, self.$c + o.$c)
            }
        }
        impl core::ops::Sub for $t {
            type Output = $t;
            #[inline]
            fn sub(self, o: $t) -> $t {
                $t::new(self.$a - o.$a, self.$b - o.$b, self.$c - o.$c)
            }
        }
        impl core::ops::Mul<f64> for $t {
            type Output = $t;
            #[inline]
            fn mul(self, k: f64) -> $t {
                $t::new(self.$a * k, self.$b * k, self.$c * k)
            }
        }
        impl core::ops::Neg for $t {
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
pub(crate) use vector_ops;

/// A plain 3D vector: a direction, an offset or a velocity. Metres (or m/s), unitless for a
/// direction. Supports `+`, `-`, `* f64`, `/ f64` and unary `-`.
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
        fma(self.x, o.x, fma(self.y, o.y, self.z * o.z))
    }

    /// Cross product (right-handed). On FMA targets each component fuses `a·b − c·d`, which halves
    /// the rounding error — what keeps the cross of two nearly parallel vectors honest.
    #[inline]
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3::new(
            fma(self.y, o.z, -(self.z * o.y)),
            fma(self.z, o.x, -(self.x * o.z)),
            fma(self.x, o.y, -(self.y * o.x)),
        )
    }

    /// Length.
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

    /// Distance between two points.
    #[inline]
    pub fn distance_to(self, o: Vec3) -> f64 {
        (self - o).norm()
    }

    /// The angle between two vectors, **radians** in `[0, π]`; `0` if either is the zero vector.
    /// Full precision at every angle, including tiny ones (it uses `atan2(|a×b|, a·b)`, not `acos`).
    #[inline]
    pub fn angle_to(self, o: Vec3) -> f64 {
        self.cross(o).norm().atan2(self.dot(o))
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

    /// As `[x, y, z]`.
    #[inline]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

vector_ops!(Vec3 { x, y, z });

impl Div<f64> for Vec3 {
    type Output = Vec3;
    #[inline]
    fn div(self, k: f64) -> Vec3 {
        Vec3::new(self.x / k, self.y / k, self.z / k)
    }
}

#[cfg(test)]
mod tests {
    /// `two_prod` is exact: its two parts add up to the true product (checked against a fused
    /// multiply–add), and double-double sums keep what plain `f64` loses.
    #[test]
    fn exact_products_and_double_doubles() {
        let mut s: u64 = 0x2545_f491_4f6c_dd1d;
        for _ in 0..100_000 {
            let mut r = || {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                let m = (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5;
                m * 10f64.powi((s % 61) as i32 - 30)
            };
            let (a, b) = (r(), r());
            let (hi, lo) = two_prod(a, b);
            assert_eq!(hi, a * b);
            assert_eq!(lo, a.mul_add(b, -hi), "{a} × {b}");
        }
        // (1 + 2⁻³⁰)² − 1 − 2⁻²⁹ = 2⁻⁶⁰: lost in f64, kept in double-double.
        let x = 1.0 + 2f64.powi(-30);
        assert_eq!(x * x - 1.0 - 2f64.powi(-29), 0.0);
        let dd = Dd::prod(x, x) - Dd(1.0, 0.0) - Dd(2f64.powi(-29), 0.0);
        assert_eq!(dd.f64(), 2f64.powi(-60));
        assert_eq!((Dd(3.0, 0.0) * Dd(5.0, 0.0)).f64(), 15.0);
        let v = Vec3::new(1e8, 1.0, 1e-8);
        let d = Dd::dot(v, v); // 10¹⁶ + 1 + 10⁻¹⁶: f64 alone keeps only the 10¹⁶
        assert_eq!((d.0, d.1), (1e16, 1.0));
    }

    use super::*;
    use crate::test_util::close;

    #[test]
    fn small_angles_keep_their_precision() {
        // acos(cos θ) loses everything below ~1e-8 rad; atan2(|a×b|, a·b) keeps full precision.
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 1e-10, 0.0);
        assert!(close(a.angle_to(b), 1e-10, 1e-24), "{}", a.angle_to(b));
        assert!(close(a.angle_to(-b), core::f64::consts::PI - 1e-10, 1e-15));
        assert!(close(
            a.angle_to(Vec3::new(0.0, 3.0, 0.0)),
            core::f64::consts::FRAC_PI_2,
            1e-16
        ));
    }

    #[test]
    fn vector_math() {
        let a = Vec3::new(1.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 2.0, 0.0);
        assert!(close(a.angle_to(b), core::f64::consts::FRAC_PI_2, 1e-12));
        assert_eq!(a.angle_to(Vec3::ZERO), 0.0);
        assert_eq!(Vec3::new(3.0, 4.0, 0.0).project_onto(a), Vec3::new(3.0, 0.0, 0.0));
        assert_eq!(Vec3::new(1.0, 1.0, 1.0).project_onto(Vec3::ZERO), Vec3::ZERO);
        assert_eq!(Vec3::new(2.0, 4.0, 6.0) / 2.0, Vec3::new(1.0, 2.0, 3.0));
        assert!(close(Vec3::new(3.0, 4.0, 0.0).distance_to(Vec3::ZERO), 5.0, 1e-12));
        // A unit vector's norm is 1 only to within rounding (here (1,2,2)/3), not bit-exact.
        assert!(close(Vec3::new(1.0, 2.0, 2.0).normalized().unwrap().norm(), 1.0, 1e-15));
        assert!(Vec3::ZERO.normalized().is_none());
        assert_eq!(<[f64; 3]>::from(Vec3::new(1.0, 2.0, 3.0)), [1.0, 2.0, 3.0]);
        assert_eq!(a.cross(b), Vec3::new(0.0, 0.0, 2.0));
        assert_eq!(a + b - b, a);
        assert_eq!(-(a * 2.0), Vec3::from([-2.0, 0.0, 0.0]));
    }
}
