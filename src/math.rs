//! The floating-point maths that the standard library provides as methods on `f64` (`sin`, `sqrt`,
//! ...), supplied by the `libm` crate for builds without it. The methods have the same names, so
//! the rest of the crate reads the same either way: a module brings them into scope with
//! `use crate::math::Float as _;`, and with the standard library present, `f64`'s own methods are
//! used instead.

pub(crate) trait Float {
    fn sqrt(self) -> f64;
    fn cbrt(self) -> f64;
    fn hypot(self, other: f64) -> f64;
    fn sin(self) -> f64;
    fn cos(self) -> f64;
    fn tan(self) -> f64;
    fn sin_cos(self) -> (f64, f64);
    fn atan(self) -> f64;
    fn atan2(self, other: f64) -> f64;
    fn atanh(self) -> f64;
    fn sinh(self) -> f64;
    fn asinh(self) -> f64;
    fn exp(self) -> f64;
    fn exp_m1(self) -> f64;
    fn floor(self) -> f64;
    fn rem_euclid(self, rhs: f64) -> f64;
    fn mul_add(self, a: f64, b: f64) -> f64;
}

impl Float for f64 {
    #[inline]
    fn sqrt(self) -> f64 {
        libm::sqrt(self)
    }
    #[inline]
    fn cbrt(self) -> f64 {
        libm::cbrt(self)
    }
    #[inline]
    fn hypot(self, other: f64) -> f64 {
        libm::hypot(self, other)
    }
    #[inline]
    fn sin(self) -> f64 {
        libm::sin(self)
    }
    #[inline]
    fn cos(self) -> f64 {
        libm::cos(self)
    }
    #[inline]
    fn tan(self) -> f64 {
        libm::tan(self)
    }
    #[inline]
    fn sin_cos(self) -> (f64, f64) {
        libm::sincos(self)
    }
    #[inline]
    fn atan(self) -> f64 {
        libm::atan(self)
    }
    #[inline]
    fn atan2(self, other: f64) -> f64 {
        libm::atan2(self, other)
    }
    #[inline]
    fn atanh(self) -> f64 {
        libm::atanh(self)
    }
    #[inline]
    fn sinh(self) -> f64 {
        libm::sinh(self)
    }
    #[inline]
    fn asinh(self) -> f64 {
        libm::asinh(self)
    }
    #[inline]
    fn exp(self) -> f64 {
        libm::exp(self)
    }
    #[inline]
    fn exp_m1(self) -> f64 {
        libm::expm1(self)
    }
    #[inline]
    fn floor(self) -> f64 {
        libm::floor(self)
    }
    /// As the standard library defines it.
    #[inline]
    fn rem_euclid(self, rhs: f64) -> f64 {
        let r = self % rhs;
        if r < 0.0 {
            r + rhs.abs()
        } else {
            r
        }
    }
    #[inline]
    fn mul_add(self, a: f64, b: f64) -> f64 {
        libm::fma(self, a, b)
    }
}
