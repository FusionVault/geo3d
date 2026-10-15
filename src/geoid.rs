//! Heights above sea level: the [`Geoid`] trait, and [`GeoidGrid`], which reads the standard
//! EGM84, EGM96 and EGM2008 geoid grids in GeographicLib's PGM format.

use core::fmt;

use crate::angle::{ang_normalize, lat_fix};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Geodetic;

/// The geoid: the level surface that mean sea level follows, as a height above the ellipsoid.
///
/// GPS heights ([`Geodetic::height_m`]) are above the WGS84 ellipsoid. Maps, altimeters and
/// aviation use heights above sea level. The geoid's height above the ellipsoid, `N`, links the two:
///
/// - height above sea level `H = h − N` ([`sea_level_height_m`](Geoid::sea_level_height_m));
/// - height above the ellipsoid `h = H + N` ([`geodetic_from_sea_level`](Geoid::geodetic_from_sea_level)).
///
/// `N` ranges from about −106 m (south of India) to +85 m (New Guinea).
///
/// Use [`GeoidGrid`] with one of the standard EGM grids, or any function or closure
/// `Fn(lat_deg, lon_deg) -> metres`, which is a `Geoid` too. That lets a national model, or
/// another crate's geoid, plug straight in:
///
/// ```
/// use geo3d::{Geodetic, Geoid};
///
/// let geoid = |_lat: f64, _lon: f64| 22.46; // N near Sydney, in metres
/// let fix = Geodetic::new(-33.87, 151.21, 42.5); // from GPS
/// assert!((geoid.sea_level_height_m(fix) - 20.04).abs() < 1e-9);
///
/// let roof = geoid.geodetic_from_sea_level(-33.87, 151.21, 20.04);
/// assert!((roof.height_m - 42.5).abs() < 1e-9);
/// ```
///
/// The geoid is the sea level of a calm ocean: real mean sea level differs from it by up to about
/// 2 m (currents and winds). National height datums, like NAVD88 or AHD, also differ by up to a
/// metre or so. For those, implement `Geoid` with the national grid.
pub trait Geoid {
    /// The geoid's height above the ellipsoid, `N`, in metres, at a latitude and longitude in
    /// degrees. NaN for NaN input or a latitude beyond ±90°.
    #[doc(alias = "geoid_height")]
    fn undulation_m(&self, lat_deg: f64, lon_deg: f64) -> f64;

    /// The height of `g` above the geoid (roughly, above mean sea level), metres: `h − N`.
    #[doc(alias = "orthometric_height")]
    #[doc(alias = "msl")]
    #[inline]
    fn sea_level_height_m(&self, g: Geodetic) -> f64 {
        g.height_m - self.undulation_m(g.lat_deg, g.lon_deg)
    }

    /// The position at `lat_deg`, `lon_deg` (degrees) that is `sea_level_height_m` metres above the
    /// geoid, with its height above the ellipsoid: `H + N`. The reverse of
    /// [`sea_level_height_m`](Geoid::sea_level_height_m).
    #[doc(alias = "ellipsoidal_height")]
    #[inline]
    fn geodetic_from_sea_level(&self, lat_deg: f64, lon_deg: f64, sea_level_height_m: f64) -> Geodetic {
        Geodetic::new(
            lat_deg,
            lon_deg,
            sea_level_height_m + self.undulation_m(lat_deg, lon_deg),
        )
    }
}

/// Any `Fn(lat_deg, lon_deg) -> metres` is a geoid.
impl<F: Fn(f64, f64) -> f64> Geoid for F {
    #[inline]
    fn undulation_m(&self, lat_deg: f64, lon_deg: f64) -> f64 {
        self(lat_deg, lon_deg)
    }
}

/// How [`GeoidGrid`] fills in between grid points.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum GeoidInterpolation {
    /// A cubic fitted to 12 nearby points (the default, and GeographicLib's). On egm96-5 it's
    /// within 3 mm of the EGM96 model everywhere.
    #[default]
    Cubic,
    /// Straight lines between the 4 corners of the cell: quicker, and continuous across cell edges,
    /// but less accurate. On egm96-5 it's up to 12 cm out, against the cubic's 3 mm. Use it to match
    /// other software that interpolates bilinearly.
    Bilinear,
}

/// Why [`GeoidGrid::from_pgm`] couldn't read a file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum GeoidGridError {
    /// The bytes don't start with a binary PGM header: `P5`, width, height and maxval.
    NotPgm,
    /// The header has no `# Offset` or `# Scale` comment, or the scale isn't positive.
    NoOffsetOrScale,
    /// The pixels aren't 16-bit (maxval 65535).
    Not16Bit,
    /// The grid isn't a whole-world grid: the width must be even, and the height odd and at least 3.
    NotGlobal,
    /// The data after the header isn't exactly width × height × 2 bytes.
    WrongLength,
}

impl fmt::Display for GeoidGridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            GeoidGridError::NotPgm => "not a binary PGM file (it should start with P5)",
            GeoidGridError::NoOffsetOrScale => "the PGM header has no Offset or Scale, or the scale isn't positive",
            GeoidGridError::Not16Bit => "the PGM pixels aren't 16-bit (maxval 65535)",
            GeoidGridError::NotGlobal => "not a whole-world grid (the width must be even and the height odd)",
            GeoidGridError::WrongLength => "the PGM data isn't width × height × 2 bytes long",
        })
    }
}

impl core::error::Error for GeoidGridError {}

/// A whole-world geoid grid in GeographicLib's PGM format: EGM84, EGM96 or EGM2008 at any of the
/// published resolutions. The answers match GeographicLib's `Geoid` class and its GeoidEval tool.
///
/// The grid reads straight from the file's bytes, without copying or allocating, so it holds any
/// `B: AsRef<[u8]>`: a `Vec<u8>` from `std::fs::read`, a `&'static [u8]` from `include_bytes!` on a
/// microcontroller, or a memory-mapped file. Lookups don't change it, so threads can share one.
///
/// Download the grids from
/// [GeographicLib](https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst). Their
/// errors, against the gravity model they come from, with the default cubic interpolation:
///
/// | Grid | Size | Max error | RMS error |
/// |---|---|---|---|
/// | egm96-15 | 2.1 MB | 0.17 m | 7 mm |
/// | egm96-5 | 19 MB | 3 mm | 0.7 mm |
/// | egm2008-5 | 19 MB | 0.29 m | 4.5 mm |
/// | egm2008-2_5 | 75 MB | 31 mm | 0.8 mm |
/// | egm2008-1 | 470 MB | 2.2 mm | 0.7 mm |
///
/// [`max_error_m`](Self::max_error_m) and [`rms_error_m`](Self::rms_error_m) read these from the
/// file. EGM2008 is the most accurate geoid, and EGM96 the most widely used: most GPS receivers give
/// heights above it. The two differ by about 0.5 m RMS.
///
/// ```no_run
/// use geo3d::{Geodetic, Geoid, GeoidGrid};
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let egm96 = GeoidGrid::from_pgm(std::fs::read("egm96-5.pgm")?)?;
/// let fix = Geodetic::new(-33.87, 151.21, 42.5); // from GPS
/// let n = egm96.undulation_m(fix.lat_deg, fix.lon_deg); // 22.4595 m
/// let above_sea_level = egm96.sea_level_height_m(fix); // 20.04 m
/// # Ok(()) }
/// ```
#[derive(Clone, Debug)]
pub struct GeoidGrid<B> {
    data: B,
    /// Byte offset of the first pixel.
    start: usize,
    width: usize,
    height: usize,
    offset: f64,
    scale: f64,
    /// Columns per degree of longitude, and rows per degree of latitude.
    lon_res: f64,
    lat_res: f64,
    interpolation: GeoidInterpolation,
    /// Byte range of the description in the header.
    description: Option<(usize, usize)>,
    /// Max and RMS error for bilinear, then for cubic, from the header.
    errors: [Option<f64>; 4],
}

impl<B: AsRef<[u8]>> GeoidGrid<B> {
    /// Read a GeographicLib `.pgm` geoid file from its bytes, with cubic interpolation. Only the
    /// header is read now; the rest is read as it's needed.
    ///
    /// Fails with a [`GeoidGridError`] unless the bytes are a whole-world, 16-bit PGM grid whose
    /// header gives its `Offset` and `Scale`, as GeographicLib's files do.
    pub fn from_pgm(data: B) -> Result<Self, GeoidGridError> {
        let bytes = data.as_ref();
        let mut header = Header::default();
        // "P5", then whitespace or a comment.
        if !bytes.starts_with(b"P5") || !bytes.get(2).is_some_and(|&b| b.is_ascii_whitespace() || b == b'#') {
            return Err(GeoidGridError::NotPgm);
        }
        let mut pos = 2;
        let mut numbers = [0usize; 3]; // width, height, maxval
        for n in &mut numbers {
            pos = header.skip_space_and_comments(bytes, pos);
            let digits = bytes[pos..].iter().take_while(|b| b.is_ascii_digit()).count();
            *n = core::str::from_utf8(&bytes[pos..pos + digits])
                .ok()
                .and_then(|s| s.parse().ok())
                .ok_or(GeoidGridError::NotPgm)?;
            pos += digits;
        }
        // Exactly one whitespace byte separates the header from the pixels.
        if !bytes.get(pos).is_some_and(u8::is_ascii_whitespace) {
            return Err(GeoidGridError::NotPgm);
        }
        let start = pos + 1;
        let [width, height, maxval] = numbers;
        let (offset, scale) = match (header.offset, header.scale) {
            (Some(o), Some(s)) if o.is_finite() && s.is_finite() && s > 0.0 => (o, s),
            _ => return Err(GeoidGridError::NoOffsetOrScale),
        };
        if maxval != 0xffff {
            return Err(GeoidGridError::Not16Bit);
        }
        // An even width lets the grid continue over the poles; an odd height puts a row on the equator.
        if width < 2 || width % 2 != 0 || height < 3 || height % 2 == 0 {
            return Err(GeoidGridError::NotGlobal);
        }
        let len = width
            .checked_mul(height)
            .and_then(|n| n.checked_mul(2))
            .and_then(|n| n.checked_add(start));
        if len != Some(bytes.len()) {
            return Err(GeoidGridError::WrongLength);
        }
        Ok(GeoidGrid {
            start,
            width,
            height,
            offset,
            scale,
            lon_res: width as f64 / 360.0,
            lat_res: (height - 1) as f64 / 180.0,
            interpolation: GeoidInterpolation::Cubic,
            description: header.description,
            errors: header.errors,
            data,
        })
    }

    /// The same grid with another [`GeoidInterpolation`].
    #[must_use]
    pub fn with_interpolation(mut self, interpolation: GeoidInterpolation) -> Self {
        self.interpolation = interpolation;
        self
    }

    /// The interpolation in use.
    pub fn interpolation(&self) -> GeoidInterpolation {
        self.interpolation
    }

    /// The description in the file, like `WGS84 EGM96, 5-minute grid`.
    pub fn description(&self) -> Option<&str> {
        let (a, b) = self.description?;
        core::str::from_utf8(self.data.as_ref().get(a..b)?).ok()
    }

    /// The grid's largest error against its gravity model with the interpolation in use, metres, as
    /// the file states it.
    pub fn max_error_m(&self) -> Option<f64> {
        self.errors[self.error_index()]
    }

    /// The grid's RMS error against its gravity model with the interpolation in use, metres, as the
    /// file states it.
    pub fn rms_error_m(&self) -> Option<f64> {
        self.errors[self.error_index() + 1]
    }

    fn error_index(&self) -> usize {
        match self.interpolation {
            GeoidInterpolation::Bilinear => 0,
            GeoidInterpolation::Cubic => 2,
        }
    }

    /// The pixel in column `ix` (wrapping round in longitude) and row `iy` (reflecting over the
    /// poles), as stored. NaN only if the bytes changed size since `from_pgm` checked them.
    fn raw(&self, data: &[u8], ix: isize, iy: isize) -> f64 {
        let (w, h) = (self.width as isize, self.height as isize);
        let mut ix = if ix < 0 {
            ix + w
        } else if ix >= w {
            ix - w
        } else {
            ix
        };
        let mut iy = iy;
        if iy < 0 || iy >= h {
            // Over the pole: the row the same distance back, on the opposite meridian.
            iy = if iy < 0 { -iy } else { 2 * (h - 1) - iy };
            ix += if ix < w / 2 { w / 2 } else { -w / 2 };
        }
        let i = self.start + 2 * (iy as usize * self.width + ix as usize);
        match data.get(i..i + 2) {
            Some(&[hi, lo]) => f64::from(u16::from_be_bytes([hi, lo])),
            _ => f64::NAN,
        }
    }
}

impl<B: AsRef<[u8]>> Geoid for GeoidGrid<B> {
    /// `N` interpolated from the grid, exactly as GeographicLib's `Geoid` does it.
    fn undulation_m(&self, lat_deg: f64, lon_deg: f64) -> f64 {
        let lat = lat_fix(lat_deg);
        let lon = ang_normalize(lon_deg);
        if lat.is_nan() || lon.is_nan() {
            return f64::NAN;
        }
        let data = self.data.as_ref();
        let half = (self.height as isize - 1) / 2;
        let mut fx = lon * self.lon_res;
        let mut fy = -lat * self.lat_res;
        let ix = fx.floor() as isize;
        let iy = (half - 1).min(fy.floor() as isize);
        fx -= ix as f64;
        fy -= iy as f64;
        let iy = iy + half;
        let w = self.width as isize;
        let ix = if ix < 0 {
            ix + w
        } else if ix >= w {
            ix - w
        } else {
            ix
        };
        let v = |dx: isize, dy: isize| self.raw(data, ix + dx, iy + dy);
        let h = match self.interpolation {
            GeoidInterpolation::Bilinear => {
                let a = (1.0 - fx) * v(0, 0) + fx * v(1, 0);
                let b = (1.0 - fx) * v(0, 1) + fx * v(1, 1);
                (1.0 - fy) * a + fy * b
            }
            GeoidInterpolation::Cubic => {
                let (wu, s0) = (self.width, self.start);
                let stencil = if ix >= 1 && ix + 2 < w && iy >= 1 && iy + 2 < self.height as isize {
                    // Away from the seam and the poles, read the four rows straight.
                    let at = |dy: isize| s0 + 2 * ((iy + dy) as usize * wu + ix as usize - 1);
                    let row = |dy: isize| data.get(at(dy)..at(dy) + 8);
                    match (row(-1), row(0), row(1), row(2)) {
                        (Some(a), Some(b), Some(c), Some(d)) => {
                            let p = |r: &[u8], k: usize| f64::from(u16::from_be_bytes([r[2 * k], r[2 * k + 1]]));
                            [
                                p(a, 1),
                                p(a, 2),
                                p(b, 0),
                                p(b, 1),
                                p(b, 2),
                                p(b, 3),
                                p(c, 0),
                                p(c, 1),
                                p(c, 2),
                                p(c, 3),
                                p(d, 1),
                                p(d, 2),
                            ]
                        }
                        _ => [f64::NAN; 12],
                    }
                } else {
                    [
                        v(0, -1),
                        v(1, -1),
                        v(-1, 0),
                        v(0, 0),
                        v(1, 0),
                        v(2, 0),
                        v(-1, 1),
                        v(0, 1),
                        v(1, 1),
                        v(2, 1),
                        v(0, 2),
                        v(1, 2),
                    ]
                };
                // At the poles, fits that don't depend on longitude there.
                let (c3, c0) = if iy == 0 {
                    (&C3N, C0N)
                } else if iy == self.height as isize - 2 {
                    (&C3S, C0S)
                } else {
                    (&C3, C0)
                };
                // The sums are of integers, so exact in any order: this one vectorises.
                let mut t = [0.0; 10];
                for (vj, row) in stencil.iter().zip(c3) {
                    for (ti, c) in t.iter_mut().zip(row) {
                        *ti += vj * c;
                    }
                }
                for ti in &mut t {
                    *ti /= c0;
                }
                t[0] + fx * (t[1] + fx * (t[3] + fx * t[6]))
                    + fy * (t[2] + fx * (t[4] + fx * t[7]) + fy * (t[5] + fx * t[8] + fy * t[9]))
            }
        };
        self.offset + self.scale * h
    }
}

/// What the header comments say.
#[derive(Default)]
struct Header {
    offset: Option<f64>,
    scale: Option<f64>,
    description: Option<(usize, usize)>,
    errors: [Option<f64>; 4],
}

impl Header {
    /// Skip whitespace and `#` comments from `pos`, noting the comments GeographicLib uses.
    fn skip_space_and_comments(&mut self, bytes: &[u8], mut pos: usize) -> usize {
        loop {
            match bytes.get(pos) {
                Some(b) if b.is_ascii_whitespace() => pos += 1,
                Some(b'#') => {
                    let end = bytes[pos..]
                        .iter()
                        .position(|&b| b == b'\n' || b == b'\r')
                        .map_or(bytes.len(), |n| pos + n);
                    self.comment(bytes, pos + 1, end);
                    pos = end;
                }
                _ => return pos,
            }
        }
    }

    fn comment(&mut self, bytes: &[u8], start: usize, end: usize) {
        let Ok(text) = core::str::from_utf8(&bytes[start..end]) else {
            return;
        };
        let mut words = text.split_ascii_whitespace();
        let (Some(key), value) = (words.next(), words.next()) else {
            return;
        };
        let number = value.and_then(|v| v.parse::<f64>().ok());
        match key {
            "Offset" => self.offset = number,
            "Scale" => self.scale = number,
            "MaxBilinearError" => self.errors[0] = number,
            "RMSBilinearError" => self.errors[1] = number,
            "MaxCubicError" => self.errors[2] = number,
            "RMSCubicError" => self.errors[3] = number,
            "Description" => {
                // Everything after the key, trimmed.
                let after_key = text.find(key).map_or(0, |k| k + key.len());
                let rest = &text[after_key..];
                let a = start + after_key + (rest.len() - rest.trim_start().len());
                let len = rest.trim().len();
                if len > 0 {
                    self.description = Some((a, a + len));
                }
            }
            _ => {}
        }
    }
}

// GeographicLib's least-squares cubic fits on the 12-point stencil (Geoid.cpp): row j gives the
// weight of stencil point j in each of the 10 cubic coefficients, over a common denominator. C3N
// and C3S are for the cells touching the north and south poles.
const C0: f64 = 240.0;
#[rustfmt::skip]
const C3: [[f64; 10]; 12] = [
    [   9.0,  -18.0,  -88.0,    0.0,   96.0,   90.0,    0.0,    0.0,  -60.0,  -20.0],
    [  -9.0,   18.0,    8.0,    0.0,  -96.0,   30.0,    0.0,    0.0,   60.0,  -20.0],
    [   9.0,  -88.0,  -18.0,   90.0,   96.0,    0.0,  -20.0,  -60.0,    0.0,    0.0],
    [ 186.0,  -42.0,  -42.0, -150.0,  -96.0, -150.0,   60.0,   60.0,   60.0,   60.0],
    [  54.0,  162.0,  -78.0,   30.0,  -24.0,  -90.0,  -60.0,   60.0,  -60.0,   60.0],
    [  -9.0,  -32.0,   18.0,   30.0,   24.0,    0.0,   20.0,  -60.0,    0.0,    0.0],
    [  -9.0,    8.0,   18.0,   30.0,  -96.0,    0.0,  -20.0,   60.0,    0.0,    0.0],
    [  54.0,  -78.0,  162.0,  -90.0,  -24.0,   30.0,   60.0,  -60.0,   60.0,  -60.0],
    [ -54.0,   78.0,   78.0,   90.0,  144.0,   90.0,  -60.0,  -60.0,  -60.0,  -60.0],
    [   9.0,   -8.0,  -18.0,  -30.0,  -24.0,    0.0,   20.0,   60.0,    0.0,    0.0],
    [  -9.0,   18.0,  -32.0,    0.0,   24.0,   30.0,    0.0,    0.0,  -60.0,   20.0],
    [   9.0,  -18.0,   -8.0,    0.0,  -24.0,  -30.0,    0.0,    0.0,   60.0,   20.0],
];
const C0N: f64 = 372.0;
#[rustfmt::skip]
const C3N: [[f64; 10]; 12] = [
    [   0.0,    0.0, -131.0,    0.0,  138.0,  144.0,    0.0,    0.0, -102.0,  -31.0],
    [   0.0,    0.0,    7.0,    0.0, -138.0,   42.0,    0.0,    0.0,  102.0,  -31.0],
    [  62.0,    0.0,  -31.0,    0.0,    0.0,  -62.0,    0.0,    0.0,    0.0,   31.0],
    [ 124.0,    0.0,  -62.0,    0.0,    0.0, -124.0,    0.0,    0.0,    0.0,   62.0],
    [ 124.0,    0.0,  -62.0,    0.0,    0.0, -124.0,    0.0,    0.0,    0.0,   62.0],
    [  62.0,    0.0,  -31.0,    0.0,    0.0,  -62.0,    0.0,    0.0,    0.0,   31.0],
    [   0.0,    0.0,   45.0,    0.0, -183.0,   -9.0,    0.0,   93.0,   18.0,    0.0],
    [   0.0,    0.0,  216.0,    0.0,   33.0,   87.0,    0.0,  -93.0,   12.0,  -93.0],
    [   0.0,    0.0,  156.0,    0.0,  153.0,   99.0,    0.0,  -93.0,  -12.0,  -93.0],
    [   0.0,    0.0,  -45.0,    0.0,   -3.0,    9.0,    0.0,   93.0,  -18.0,    0.0],
    [   0.0,    0.0,  -55.0,    0.0,   48.0,   42.0,    0.0,    0.0,  -84.0,   31.0],
    [   0.0,    0.0,   -7.0,    0.0,  -48.0,  -42.0,    0.0,    0.0,   84.0,   31.0],
];
const C0S: f64 = 372.0;
#[rustfmt::skip]
const C3S: [[f64; 10]; 12] = [
    [  18.0,  -36.0, -122.0,    0.0,  120.0,  135.0,    0.0,    0.0,  -84.0,  -31.0],
    [ -18.0,   36.0,   -2.0,    0.0, -120.0,   51.0,    0.0,    0.0,   84.0,  -31.0],
    [  36.0, -165.0,  -27.0,   93.0,  147.0,   -9.0,    0.0,  -93.0,   18.0,    0.0],
    [ 210.0,   45.0, -111.0,  -93.0,  -57.0, -192.0,    0.0,   93.0,   12.0,   93.0],
    [ 162.0,  141.0,  -75.0,  -93.0, -129.0, -180.0,    0.0,   93.0,  -12.0,   93.0],
    [ -36.0,  -21.0,   27.0,   93.0,   39.0,    9.0,    0.0,  -93.0,  -18.0,    0.0],
    [   0.0,    0.0,   62.0,    0.0,    0.0,   31.0,    0.0,    0.0,    0.0,  -31.0],
    [   0.0,    0.0,  124.0,    0.0,    0.0,   62.0,    0.0,    0.0,    0.0,  -62.0],
    [   0.0,    0.0,  124.0,    0.0,    0.0,   62.0,    0.0,    0.0,    0.0,  -62.0],
    [   0.0,    0.0,   62.0,    0.0,    0.0,   31.0,    0.0,    0.0,    0.0,  -31.0],
    [ -18.0,   36.0,  -64.0,    0.0,   66.0,   51.0,    0.0,    0.0, -102.0,   31.0],
    [  18.0,  -36.0,    2.0,    0.0,  -66.0,  -51.0,    0.0,    0.0,  102.0,   31.0],
];

#[cfg(test)]
mod tests {
    use super::*;

    const OFFSET: f64 = -108.0;
    const SCALE: f64 = 0.003;

    /// A whole-world PGM grid, `width` × `height`, with pixel values `f(column, row)`.
    fn pgm(width: usize, height: usize, f: impl Fn(usize, usize) -> u16) -> Vec<u8> {
        let mut out = format!(
            "P5\n# Description test grid\n# MaxCubicError 0.003\n# RMSCubicError 0.001\n\
             # Offset {OFFSET}\n# Scale {SCALE}\n{width} {height}\n65535\n"
        )
        .into_bytes();
        for row in 0..height {
            for col in 0..width {
                out.extend_from_slice(&f(col, row).to_be_bytes());
            }
        }
        out
    }

    /// Grid coordinates (column, row) of a latitude and longitude, on a `width` × `height` grid.
    fn xy(width: usize, height: usize, lat: f64, lon: f64) -> (f64, f64) {
        (
            lon.rem_euclid(360.0) * width as f64 / 360.0,
            (90.0 - lat) * (height - 1) as f64 / 180.0,
        )
    }

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    /// The cubic fit is exact for any cubic, and bilinear for any bilinear function, away from the
    /// poles and the seam at 0°.
    #[test]
    fn interpolation_is_exact_for_polynomials() {
        let (w, h) = (16, 13);
        let cubic = |x: f64, y: f64| {
            30_000.0 + 2.0 * x.powi(3) - 3.0 * x * x * y + x * y * y + 4.0 * y.powi(3) + 7.0 * x * x - 5.0 * x * y
                + 9.0 * y * y
                + 11.0 * x
                - 13.0 * y
        };
        let bilinear = |x: f64, y: f64| 30_000.0 + 101.0 * x - 37.0 * y + 13.0 * x * y;
        let mut seed = 1;
        for (f, interpolation) in [
            (&cubic as &dyn Fn(f64, f64) -> f64, GeoidInterpolation::Cubic),
            (&bilinear, GeoidInterpolation::Bilinear),
        ] {
            let bytes = pgm(w, h, |c, r| f(c as f64, r as f64) as u16);
            let grid = GeoidGrid::from_pgm(&bytes[..])
                .unwrap()
                .with_interpolation(interpolation);
            for _ in 0..10_000 {
                // Cells whose stencils stay inside the grid.
                let x = 1.0 + lcg(&mut seed) * (w - 4) as f64;
                let y = 1.0 + lcg(&mut seed) * (h - 4) as f64;
                let lon = x * 360.0 / w as f64;
                let lat = 90.0 - y * 180.0 / (h - 1) as f64;
                let (gx, gy) = xy(w, h, lat, lon);
                let want = OFFSET + SCALE * f(gx, gy);
                let got = grid.undulation_m(lat, lon);
                assert!(
                    (got - want).abs() < 1e-9,
                    "{interpolation:?} at {lat}, {lon}: {got} vs {want}"
                );
            }
        }
    }

    #[test]
    fn grid_points_wrap_and_poles() {
        let (w, h) = (36, 19); // 10° grid
        let bytes = pgm(w, h, |c, r| {
            if r == 0 {
                40_000 // the north pole is one point
            } else if r == h - 1 {
                20_000
            } else {
                (30_000.0 + 3_000.0 * ((c as f64) * 10f64.to_radians()).sin() + 500.0 * r as f64) as u16
            }
        });
        for interpolation in [GeoidInterpolation::Cubic, GeoidInterpolation::Bilinear] {
            let g = GeoidGrid::from_pgm(bytes.as_slice())
                .unwrap()
                .with_interpolation(interpolation);
            let raw = |c: usize, r: usize| {
                let i = bytes.len() - 2 * w * h + 2 * (r * w + c);
                OFFSET + SCALE * f64::from(u16::from_be_bytes([bytes[i], bytes[i + 1]]))
            };
            // Bilinear goes through every grid point.
            if interpolation == GeoidInterpolation::Bilinear {
                assert_eq!(g.undulation_m(0.0, 0.0), raw(0, 9));
                assert_eq!(g.undulation_m(-20.0, 250.0), raw(25, 11));
            }
            // Longitude wraps: the same point three ways.
            for lat in [-73.0, -5.5, 0.0, 41.0, 88.0] {
                let n = g.undulation_m(lat, 123.25);
                assert_eq!(n, g.undulation_m(lat, 123.25 - 360.0));
                assert_eq!(n, g.undulation_m(lat, 123.25 + 720.0));
            }
            // At the poles every longitude gives the same height.
            for (lat, want) in [(90.0, raw(0, 0)), (-90.0, raw(0, h - 1))] {
                for lon in [-180.0, -97.0, 0.0, 33.3, 180.0] {
                    let n = g.undulation_m(lat, lon);
                    assert!((n - want).abs() < 1e-9, "{interpolation:?} {lat}, {lon}: {n} vs {want}");
                }
            }
            assert!(g.undulation_m(f64::NAN, 0.0).is_nan());
            assert!(g.undulation_m(0.0, f64::NAN).is_nan());
            assert!(g.undulation_m(0.0, f64::INFINITY).is_nan());
            assert!(g.undulation_m(90.000001, 0.0).is_nan());
        }
    }

    /// Turning the grid by 7° turns the answers by 7°, exactly: across the seam at 0°, and in the
    /// cells next to the poles, whose stencils reach over to the other side.
    #[test]
    fn turning_the_grid_turns_the_answers() {
        let (w, h) = (360, 181); // 1°, so grid columns are whole degrees
        let f = |c: usize, r: usize| (20_000 + (c * 7919 + r * 104_729) % 30_000) as u16;
        let k = 7;
        let a = pgm(w, h, f);
        let b = pgm(w, h, |c, r| f((c + w - k) % w, r));
        for interpolation in [GeoidInterpolation::Cubic, GeoidInterpolation::Bilinear] {
            let ga = GeoidGrid::from_pgm(&a).unwrap().with_interpolation(interpolation);
            let gb = GeoidGrid::from_pgm(&b).unwrap().with_interpolation(interpolation);
            let mut seed = 7;
            for i in 0..20_000 {
                // Quarter-degree points, so adding 7° is exact; every 4th one is near a pole.
                let lat = if i % 4 == 0 {
                    89.0
                } else {
                    -89.0 + 178.0 * lcg(&mut seed)
                };
                let lat = (lat * 4.0).round() / 4.0 * if i % 8 == 0 { -1.0 } else { 1.0 };
                let lon = (-180.0 + 360.0 * lcg(&mut seed)) * 4.0;
                let lon = lon.round() / 4.0;
                let lon = if i % 3 == 0 { -0.25 - (i % 12) as f64 } else { lon }; // the seam
                assert_eq!(
                    ga.undulation_m(lat, lon),
                    gb.undulation_m(lat, lon + k as f64),
                    "{interpolation:?} {lat}, {lon}"
                );
            }
        }
    }

    #[test]
    fn header() {
        let bytes = pgm(4, 3, |_, _| 36_000);
        let g = GeoidGrid::from_pgm(bytes.clone()).unwrap();
        assert_eq!(g.description(), Some("test grid"));
        assert_eq!((g.max_error_m(), g.rms_error_m()), (Some(0.003), Some(0.001)));
        assert_eq!(g.interpolation(), GeoidInterpolation::Cubic);
        let b = g.with_interpolation(GeoidInterpolation::Bilinear);
        assert_eq!((b.max_error_m(), b.rms_error_m()), (None, None));
        assert!((b.undulation_m(12.0, 34.0) - (OFFSET + SCALE * 36_000.0)).abs() < 1e-12);

        // GeographicLib's real egm96-5 header, with Windows line endings for good measure.
        let real = "P5\r\n# Geoid file in PGM format for the GeographicLib::Geoid class\r\n\
            # Description WGS84 EGM96, 5-minute grid\r\n\
            # URL https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84#tab_egm96\r\n\
            # DateTime 2009-08-29 18:45:03\r\n# MaxBilinearError 0.140\r\n# RMSBilinearError 0.005\r\n\
            # MaxCubicError 0.003\r\n# RMSCubicError 0.001\r\n# Offset -108\r\n# Scale 0.003\r\n\
            # Origin 90N 0E\r\n# AREA_OR_POINT Point\r\n# Vertical_Datum WGS84\r\n4320    2161\r\n65535\n";
        let mut bytes = real.as_bytes().to_vec();
        bytes.resize(real.len() + 4320 * 2161 * 2, 0);
        let g = GeoidGrid::from_pgm(&bytes).unwrap();
        assert_eq!(g.description(), Some("WGS84 EGM96, 5-minute grid"));
        assert_eq!((g.max_error_m(), g.rms_error_m()), (Some(0.003), Some(0.001)));
        assert_eq!(g.undulation_m(10.0, 20.0), -108.0);
    }

    #[test]
    fn bad_files_are_refused() {
        let good = pgm(4, 3, |_, _| 0);
        assert!(GeoidGrid::from_pgm(&good).is_ok());
        let with = |from: &str, to: &str| {
            let text = String::from_utf8_lossy(&good[..good.len() - 24]).replace(from, to);
            let mut bytes = text.into_bytes();
            bytes.extend_from_slice(&good[good.len() - 24..]);
            GeoidGrid::from_pgm(bytes).err()
        };
        use GeoidGridError::*;
        assert_eq!(with("P5", "P2"), Some(NotPgm));
        assert_eq!(with("# Scale 0.003", ""), Some(NoOffsetOrScale));
        assert_eq!(with("# Scale 0.003", "# Scale -1"), Some(NoOffsetOrScale));
        assert_eq!(with("# Offset -108", "# Offset x"), Some(NoOffsetOrScale));
        assert_eq!(with("65535", "255"), Some(Not16Bit));
        assert_eq!(with("4 3", "3 4"), Some(NotGlobal));
        assert_eq!(with("4 3", "2 6"), Some(NotGlobal));
        assert_eq!(with("4 3", "4 1"), Some(NotGlobal));
        assert_eq!(with("4 3", "6 3"), Some(WrongLength));
        assert_eq!(with("4 3", "99999999999999999999 3"), Some(NotPgm));
        // A width that parses, but whose grid can't fit in memory, on 32- and 64-bit alike.
        let huge = format!("{} 3", usize::MAX / 2 + 1);
        assert_eq!(with("4 3", &huge), Some(WrongLength));
        assert_eq!(with("65535\n", "65535"), Some(NotPgm));
        assert_eq!(GeoidGrid::from_pgm(&good[..good.len() - 1]).err(), Some(WrongLength));
        assert_eq!(GeoidGrid::from_pgm(b"").err(), Some(NotPgm));
        assert_eq!(GeoidGrid::from_pgm(b"P5").err(), Some(NotPgm));
        assert_eq!(GeoidGrid::from_pgm(b"P5 # Offset 1").err(), Some(NotPgm));
        assert!(format!("{}", WrongLength).contains("width × height × 2"));
    }

    #[test]
    fn closures_and_conversions() {
        let geoid = |lat: f64, _lon: f64| 10.0 + lat;
        let g = Geodetic::new(20.0, 5.0, 100.0);
        assert_eq!(geoid.sea_level_height_m(g), 70.0);
        let back = geoid.geodetic_from_sea_level(20.0, 5.0, 70.0);
        assert_eq!(back, g);
        // A trait object works too.
        let dynamic: &dyn Geoid = &geoid;
        assert_eq!(dynamic.undulation_m(1.0, 0.0), 11.0);
    }
}
