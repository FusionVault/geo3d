//! `Display` for the coordinate types, so `println!("{p}")` gives something readable.
//!
//! - **`{}`** prints every number losslessly, the way Rust prints an `f64`: the shortest digits
//!   that read back as exactly the same value.
//! - **`{:.N}`** prints `N` decimals. Where a type mixes degrees and metres ([`Geodetic`], [`Aer`]),
//!   the metres get four fewer, since 1e-4° is about 10 m on the ground.
//! - **`{:#}`** on a [`Geodetic`] writes the hemispheres as letters: `33.87° S, 151.21° E`.
//!
//! Negative zero prints as `0`. Width and alignment flags are ignored.

use core::fmt;

use crate::attitude::{Attitude, Frd};
use crate::local::{Aer, Enu, Ned};
use crate::position::{Ecef, Eci, Geodetic};
use crate::projection::UtmUps;
use crate::quaternion::Quaternion;
use crate::vector::Vec3;

/// Write `x` with `decimals` places, or losslessly if `None`; never `-0`.
fn num(f: &mut fmt::Formatter<'_>, x: f64, decimals: Option<usize>) -> fmt::Result {
    let x = x + 0.0; // -0 + 0 = +0
    match decimals {
        Some(d) => write!(f, "{x:.d$}"),
        None => write!(f, "{x}"),
    }
}

/// The decimals for metres in a type that also has degrees: four fewer than the degrees get.
fn metres_decimals(f: &fmt::Formatter<'_>) -> Option<usize> {
    f.precision().map(|p| p.saturating_sub(4))
}

/// Write `label value unit` for each field, comma separated, all with the formatter's precision.
fn fields(f: &mut fmt::Formatter<'_>, parts: [(&str, f64); 3], unit: &str) -> fmt::Result {
    let decimals = f.precision();
    for (i, (label, value)) in parts.into_iter().enumerate() {
        if i > 0 {
            f.write_str(", ")?;
        }
        write!(f, "{label} ")?;
        num(f, value, decimals)?;
        f.write_str(unit)?;
    }
    Ok(())
}

/// `-33.87°, 151.21°, 20 m`: latitude, longitude, height. With `{:#}`, `33.87° S, 151.21° E, 20 m`.
impl fmt::Display for Geodetic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision();
        if f.alternate() {
            for (value, plus, minus) in [(self.lat_deg, " N", " S"), (self.lon_deg, " E", " W")] {
                num(f, value.abs(), decimals)?;
                f.write_str("°")?;
                if !value.is_nan() {
                    // Zero counts as north and east.
                    f.write_str(if value.is_sign_negative() && value != 0.0 {
                        minus
                    } else {
                        plus
                    })?;
                }
                f.write_str(", ")?;
            }
        } else {
            num(f, self.lat_deg, decimals)?;
            f.write_str("°, ")?;
            num(f, self.lon_deg, decimals)?;
            f.write_str("°, ")?;
        }
        num(f, self.height_m, metres_decimals(f))?;
        f.write_str(" m")
    }
}

/// `az 35.6°, el 7.3°, 19254 m`: azimuth, elevation, range.
impl fmt::Display for Aer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision();
        f.write_str("az ")?;
        num(f, self.azimuth_deg, decimals)?;
        f.write_str("°, el ")?;
        num(f, self.elevation_deg, decimals)?;
        f.write_str("°, ")?;
        num(f, self.range_m, metres_decimals(f))?;
        f.write_str(" m")
    }
}

/// `x -4646032 m, y 2553122 m, z -3534494 m`.
impl fmt::Display for Ecef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fields(f, [("x", self.x), ("y", self.y), ("z", self.z)], " m")
    }
}

/// `x 6524863 m, y 0 m, z 0 m` (in the inertial frame).
impl fmt::Display for Eci {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fields(f, [("x", self.x), ("y", self.y), ("z", self.z)], " m")
    }
}

/// `east 0 m, north 11131 m, up 2960 m`.
impl fmt::Display for Enu {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fields(f, [("east", self.east), ("north", self.north), ("up", self.up)], " m")
    }
}

/// `north 11131 m, east 0 m, down -2960 m`.
impl fmt::Display for Ned {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fields(
            f,
            [("north", self.north), ("east", self.east), ("down", self.down)],
            " m",
        )
    }
}

/// `forward 500 m, right 0 m, down 30 m`.
impl fmt::Display for Frd {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fields(
            f,
            [("forward", self.forward), ("right", self.right), ("down", self.down)],
            " m",
        )
    }
}

/// `yaw 90°, pitch 10°, roll 0°`.
impl fmt::Display for Attitude {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fields(
            f,
            [
                ("yaw", self.yaw_deg),
                ("pitch", self.pitch_deg),
                ("roll", self.roll_deg),
            ],
            "°",
        )
    }
}

/// `(1, 2, 3)`: no unit, since a vector can be metres, m/s or a direction.
impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision();
        f.write_str("(")?;
        num(f, self.x, decimals)?;
        f.write_str(", ")?;
        num(f, self.y, decimals)?;
        f.write_str(", ")?;
        num(f, self.z, decimals)?;
        f.write_str(")")
    }
}

/// `0.7071 + 0i + 0j + 0.7071k`, with minus signs where the parts are negative.
impl fmt::Display for Quaternion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision();
        num(f, self.w, decimals)?;
        for (value, unit) in [(self.x, "i"), (self.y, "j"), (self.z, "k")] {
            let value = value + 0.0;
            f.write_str(if value.is_sign_negative() { " - " } else { " + " })?;
            num(f, value.abs(), decimals)?;
            f.write_str(unit)?;
        }
        Ok(())
    }
}

/// GeographicLib's form: zone and hemisphere, then easting and northing in metres, as in
/// `56s 334435.706 6250816.398`. UPS has no zone number: `n 2000000 2000000`.
impl fmt::Display for UtmUps {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision();
        if self.zone != 0 {
            write!(f, "{}", self.zone)?;
        }
        f.write_str(if self.north { "n " } else { "s " })?;
        num(f, self.easting_m, decimals)?;
        f.write_str(" ")?;
        num(f, self.northing_m, decimals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_is_lossless_and_precision_rounds() {
        let g = Geodetic::new(-33.87, 151.21, 20.0);
        assert_eq!(format!("{g}"), "-33.87°, 151.21°, 20 m");
        assert_eq!(format!("{g:.5}"), "-33.87000°, 151.21000°, 20.0 m");
        assert_eq!(format!("{g:.7}"), "-33.8700000°, 151.2100000°, 20.000 m");
        assert_eq!(format!("{g:.2}"), "-33.87°, 151.21°, 20 m");
        assert_eq!(format!("{g:#}"), "33.87° S, 151.21° E, 20 m");
        assert_eq!(
            format!("{:#.1}", Geodetic::new(51.5, -0.12, 30.0)),
            "51.5° N, 0.1° W, 30 m"
        );
        let a = Aer::new(35.6, 7.3, 19_254.4);
        assert_eq!(format!("{a:.1}"), "az 35.6°, el 7.3°, 19254 m");
        assert_eq!(format!("{a}"), "az 35.6°, el 7.3°, 19254.4 m");
    }

    #[test]
    fn every_type() {
        assert_eq!(format!("{}", Ecef::new(1.5, -2.0, 3.0)), "x 1.5 m, y -2 m, z 3 m");
        assert_eq!(format!("{:.1}", Eci::new(1.0, 2.0, 3.0)), "x 1.0 m, y 2.0 m, z 3.0 m");
        assert_eq!(
            format!("{}", Enu::new(0.0, 11_131.0, 2_960.0)),
            "east 0 m, north 11131 m, up 2960 m"
        );
        assert_eq!(
            format!("{}", Ned::new(1.0, 2.0, -3.0)),
            "north 1 m, east 2 m, down -3 m"
        );
        assert_eq!(
            format!("{}", Frd::new(500.0, 0.0, 30.0)),
            "forward 500 m, right 0 m, down 30 m"
        );
        assert_eq!(
            format!("{}", Attitude::new(90.0, 10.0, 0.0)),
            "yaw 90°, pitch 10°, roll 0°"
        );
        assert_eq!(format!("{:.2}", Vec3::new(1.0, 2.0, 3.0)), "(1.00, 2.00, 3.00)");
        assert_eq!(
            format!("{:.1}", Quaternion::new(0.5, -0.5, 0.5, -0.5)),
            "0.5 - 0.5i + 0.5j - 0.5k"
        );
        assert_eq!(
            format!("{:.3}", UtmUps::new(56, false, 334_435.706_14, 6_250_816.397_8)),
            "56s 334435.706 6250816.398"
        );
        assert_eq!(format!("{}", UtmUps::new(0, true, 2e6, 2e6)), "n 2000000 2000000");
    }

    #[test]
    fn zero_nan_and_infinity() {
        assert_eq!(format!("{}", Geodetic::new(-0.0, -0.0, -0.0)), "0°, 0°, 0 m");
        assert_eq!(format!("{:#}", Geodetic::new(-0.0, 0.0, 0.0)), "0° N, 0° E, 0 m");
        assert_eq!(format!("{:#}", Geodetic::new(f64::NAN, -1.0, 0.0)), "NaN°, 1° W, 0 m");
        assert_eq!(
            format!("{}", Ecef::new(f64::INFINITY, f64::NEG_INFINITY, 0.0)),
            "x inf m, y -inf m, z 0 m"
        );
        assert_eq!(format!("{}", Quaternion::new(1.0, -0.0, 0.0, 0.0)), "1 + 0i + 0j + 0k");
    }

    /// The plain form reads back as exactly the same numbers, for any value.
    #[test]
    fn plain_form_round_trips_exactly() {
        let mut s: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut u = || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            let bits = s.wrapping_mul(0x2545_f491_4f6c_dd1d);
            // Any 53-bit mantissa, times 10^-20 to 10^20, either sign.
            let m = (bits >> 11) as f64 / (1u64 << 53) as f64;
            let e = (bits % 41) as i32 - 20;
            if bits & 1 == 0 {
                m * 10f64.powi(e)
            } else {
                -m * 10f64.powi(e)
            }
        };
        for _ in 0..10_000 {
            let g = Geodetic::new(u(), u(), u());
            let text = format!("{g}");
            let numbers: Vec<f64> = text
                .split(['°', ',', 'm'])
                .map(str::trim)
                .filter(|t| !t.is_empty())
                .map(|t| t.parse().unwrap())
                .collect();
            assert_eq!(numbers, [g.lat_deg + 0.0, g.lon_deg + 0.0, g.height_m + 0.0], "{text}");
        }
    }
}
