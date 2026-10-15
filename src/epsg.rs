//! Every ellipsoid in the EPSG registry, by its code: [`Ellipsoid::from_epsg`].

use crate::ellipsoid::*;

/// EPSG's older units, in metres (EPSG units 9005, 9080, 9031 and 9002).
const CLARKE_FOOT: f64 = 0.304_797_265_4;
const INDIAN_FOOT: f64 = 12.0 / 39.370_142;
const GERMAN_LEGAL_METRE: f64 = 1.000_013_596_5;
const FOOT: f64 = 0.304_8;

impl Ellipsoid {
    /// The ellipsoid with this code in the [EPSG registry](https://epsg.org), the standard
    /// catalogue of datums and coordinate systems; `None` for any other code. All 50 are here,
    /// each with the numbers EPSG defines it by (as of October 2026), converted to metres where
    /// EPSG gives feet.
    ///
    /// ```
    /// use geo3d::{Ellipsoid, KRASSOVSKY_1940};
    /// assert_eq!(Ellipsoid::from_epsg(7024), Some(KRASSOVSKY_1940));
    /// let everest = Ellipsoid::from_epsg(7044).unwrap(); // Everest 1830 (1962 Definition), Pakistan
    /// assert_eq!(everest.a(), 6_377_301.243);
    /// assert_eq!(Ellipsoid::from_epsg(4326), None); // that's a coordinate system, not an ellipsoid
    /// ```
    ///
    /// | Code | Ellipsoid | Defined by | Named |
    /// |---|---|---|---|
    /// | 1024 | CGCS2000 | `a`, `1/f` | [`CGCS2000`](crate::CGCS2000) |
    /// | 1025 | GSK-2011 | `a`, `1/f` | [`GSK_2011`](crate::GSK_2011) |
    /// | 1026 | Zach 1812 | `a`, `1/f` |  |
    /// | 7001 | Airy 1830 | `a`, `1/f` | [`AIRY_1830`](crate::AIRY_1830) |
    /// | 7002 | Airy Modified 1849 | `a`, `1/f` |  |
    /// | 7003 | Australian National Spheroid | `a`, `1/f` | [`AUSTRALIAN_NATIONAL`](crate::AUSTRALIAN_NATIONAL) |
    /// | 7004 | Bessel 1841 | `a`, `1/f` | [`BESSEL_1841`](crate::BESSEL_1841) |
    /// | 7005 | Bessel Modified | `a`, `1/f` |  |
    /// | 7007 | Clarke 1858 | axes, in Clarke's feet |  |
    /// | 7008 | Clarke 1866 | axes | [`CLARKE_1866`](crate::CLARKE_1866) |
    /// | 7010 | Clarke 1880 (Benoit) | axes |  |
    /// | 7011 | Clarke 1880 (IGN) | axes |  |
    /// | 7012 | Clarke 1880 (RGS) | `a`, `1/f` |  |
    /// | 7013 | Clarke 1880 (Arc) | `a`, `1/f` |  |
    /// | 7014 | Clarke 1880 (SGA 1922) | `a`, `1/f` |  |
    /// | 7015 | Everest 1830 (1937 Adjustment) | `a`, `1/f` | [`EVEREST_1830`](crate::EVEREST_1830) |
    /// | 7016 | Everest 1830 (1967 Definition) | `a`, `1/f` |  |
    /// | 7018 | Everest 1830 Modified | `a`, `1/f` | [`EVEREST_1830_MODIFIED`](crate::EVEREST_1830_MODIFIED) |
    /// | 7019 | GRS 1980 | `a`, `1/f` | [`GRS80`](crate::GRS80) |
    /// | 7020 | Helmert 1906 | `a`, `1/f` |  |
    /// | 7021 | Indonesian National Spheroid | `a`, `1/f` |  |
    /// | 7022 | International 1924 | `a`, `1/f` | [`INTERNATIONAL_1924`](crate::INTERNATIONAL_1924) |
    /// | 7024 | Krassowsky 1940 | `a`, `1/f` | [`KRASSOVSKY_1940`](crate::KRASSOVSKY_1940) |
    /// | 7025 | NWL 9D | `a`, `1/f` | [`WGS66`](crate::WGS66) |
    /// | 7027 | Plessis 1817 | `a`, `1/f` |  |
    /// | 7028 | Struve 1860 | `a`, `1/f` |  |
    /// | 7029 | War Office | `a`, `1/f` |  |
    /// | 7030 | WGS 84 | `a`, `1/f` | [`WGS84`](crate::WGS84) |
    /// | 7031 | GEM 10C | `a`, `1/f` |  |
    /// | 7032 | OSU86F | `a`, `1/f` |  |
    /// | 7033 | OSU91A | `a`, `1/f` |  |
    /// | 7034 | Clarke 1880 | axes, in Clarke's feet |  |
    /// | 7036 | GRS 1967 | `a`, `1/f` | [`GRS67`](crate::GRS67) |
    /// | 7041 | Average Terrestrial System 1977 | `a`, `1/f` |  |
    /// | 7042 | Everest (1830 Definition) | axes, in Indian feet |  |
    /// | 7043 | WGS 72 | `a`, `1/f` | [`WGS72`](crate::WGS72) |
    /// | 7044 | Everest 1830 (1962 Definition) | `a`, `1/f` |  |
    /// | 7045 | Everest 1830 (1975 Definition) | `a`, `1/f` |  |
    /// | 7046 | Bessel Namibia (GLM) | `a`, `1/f`, in German legal metres |  |
    /// | 7048 | GRS 1980 Authalic Sphere | axes |  |
    /// | 7049 | IAG 1975 | `a`, `1/f` | [`IAG_1975`](crate::IAG_1975) |
    /// | 7050 | GRS 1967 Modified | `a`, `1/f` |  |
    /// | 7051 | Danish 1876 | `a`, `1/f` |  |
    /// | 7052 | Clarke 1866 Authalic Sphere | axes |  |
    /// | 7053 | Hough 1960 | `a`, `1/f` |  |
    /// | 7054 | PZ-90 | `a`, `1/f` | [`PZ90`](crate::PZ90) |
    /// | 7055 | Clarke 1880 (international foot) | axes, in international feet |  |
    /// | 7056 | Everest 1830 (RSO 1969) | `a`, `1/f` |  |
    /// | 7057 | International 1924 Authalic Sphere | axes |  |
    /// | 7058 | Hughes 1980 | axes |  |
    pub const fn from_epsg(code: u32) -> Option<Ellipsoid> {
        Some(match code {
            1024 => CGCS2000,                                                                      // CGCS2000
            1025 => GSK_2011,                                                                      // GSK-2011
            1026 => Ellipsoid::new(6_376_045.0, 1.0 / 310.0),                                      // Zach 1812
            7001 => AIRY_1830,                                                                     // Airy 1830
            7002 => Ellipsoid::new(6_377_340.189, 1.0 / 299.324_964_6),                            // Airy Modified 1849
            7003 => AUSTRALIAN_NATIONAL, // Australian National Spheroid
            7004 => BESSEL_1841,         // Bessel 1841
            7005 => Ellipsoid::new(6_377_492.018, 1.0 / 299.152_812_8), // Bessel Modified
            7007 => Ellipsoid::from_axes(20_926_348.0 * CLARKE_FOOT, 20_855_233.0 * CLARKE_FOOT), // Clarke 1858
            7008 => CLARKE_1866,         // Clarke 1866
            7010 => Ellipsoid::from_axes(6_378_300.789, 6_356_566.435), // Clarke 1880 (Benoit)
            7011 => Ellipsoid::from_axes(6_378_249.2, 6_356_515.0), // Clarke 1880 (IGN)
            7012 => Ellipsoid::new(6_378_249.145, 1.0 / 293.465), // Clarke 1880 (RGS)
            7013 => Ellipsoid::new(6_378_249.145, 1.0 / 293.466_307_7), // Clarke 1880 (Arc)
            7014 => Ellipsoid::new(6_378_249.2, 1.0 / 293.465_98), // Clarke 1880 (SGA 1922)
            7015 => EVEREST_1830,        // Everest 1830 (1937 Adjustment)
            7016 => Ellipsoid::new(6_377_298.556, 1.0 / 300.801_7), // Everest 1830 (1967 Definition)
            7018 => EVEREST_1830_MODIFIED, // Everest 1830 Modified
            7019 => GRS80,               // GRS 1980
            7020 => Ellipsoid::new(6_378_200.0, 1.0 / 298.3), // Helmert 1906
            7021 => Ellipsoid::new(6_378_160.0, 1.0 / 298.247), // Indonesian National Spheroid
            7022 => INTERNATIONAL_1924,  // International 1924
            7024 => KRASSOVSKY_1940,     // Krassowsky 1940
            7025 => WGS66,               // NWL 9D
            7027 => Ellipsoid::new(6_376_523.0, 1.0 / 308.64), // Plessis 1817
            7028 => Ellipsoid::new(6_378_298.3, 1.0 / 294.73), // Struve 1860
            7029 => Ellipsoid::new(6_378_300.0, 1.0 / 296.0), // War Office
            7030 => WGS84,               // WGS 84
            7031 => Ellipsoid::new(6_378_137.0, 1.0 / 298.257_223_563), // GEM 10C
            7032 => Ellipsoid::new(6_378_136.2, 1.0 / 298.257_223_563), // OSU86F
            7033 => Ellipsoid::new(6_378_136.3, 1.0 / 298.257_223_563), // OSU91A
            7034 => Ellipsoid::from_axes(20_926_202.0 * CLARKE_FOOT, 20_854_895.0 * CLARKE_FOOT), // Clarke 1880
            7036 => GRS67,               // GRS 1967
            7041 => Ellipsoid::new(6_378_135.0, 1.0 / 298.257), // Average Terrestrial System 1977
            7042 => Ellipsoid::from_axes(20_922_931.8 * INDIAN_FOOT, 20_853_374.58 * INDIAN_FOOT), // Everest (1830 Definition)
            7043 => WGS72,                                                                         // WGS 72
            7044 => Ellipsoid::new(6_377_301.243, 1.0 / 300.801_725_5), // Everest 1830 (1962 Definition)
            7045 => Ellipsoid::new(6_377_299.151, 1.0 / 300.801_725_5), // Everest 1830 (1975 Definition)
            7046 => Ellipsoid::new(6_377_397.155 * GERMAN_LEGAL_METRE, 1.0 / 299.152_812_8), // Bessel Namibia (GLM)
            7048 => Ellipsoid::from_axes(6_371_007.0, 6_371_007.0),     // GRS 1980 Authalic Sphere
            7049 => IAG_1975,                                           // IAG 1975
            7050 => Ellipsoid::new(6_378_160.0, 1.0 / 298.25),          // GRS 1967 Modified
            7051 => Ellipsoid::new(6_377_019.27, 1.0 / 300.0),          // Danish 1876
            7052 => Ellipsoid::from_axes(6_370_997.0, 6_370_997.0),     // Clarke 1866 Authalic Sphere
            7053 => Ellipsoid::new(6_378_270.0, 1.0 / 297.0),           // Hough 1960
            7054 => PZ90,                                               // PZ-90
            7055 => Ellipsoid::from_axes(20_926_202.0 * FOOT, 20_854_895.0 * FOOT), // Clarke 1880 (international foot)
            7056 => Ellipsoid::new(6_377_295.664, 1.0 / 300.801_7),     // Everest 1830 (RSO 1969)
            7057 => Ellipsoid::from_axes(6_371_228.0, 6_371_228.0),     // International 1924 Authalic Sphere
            7058 => Ellipsoid::from_axes(6_378_273.0, 6_356_889.449),   // Hughes 1980
            _ => return None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The EPSG registry's own numbers (apps.epsg.org, October 2026), straight from its API:
    /// code, semi-major axis, semi-minor axis or NaN, inverse flattening or NaN, and the unit in
    /// metres.
    #[rustfmt::skip]
    #[allow(clippy::excessive_precision)]
    const REGISTRY: [(u32, f64, f64, f64, f64); 50] = [
    (1024, 6378137.0, f64::NAN, 298.257222101, 1.0),
    (1025, 6378136.5, f64::NAN, 298.2564151, 1.0),
    (1026, 6376045.0, f64::NAN, 310.0, 1.0),
    (7001, 6377563.396, f64::NAN, 299.3249646, 1.0),
    (7002, 6377340.189, f64::NAN, 299.3249646, 1.0),
    (7003, 6378160.0, f64::NAN, 298.25, 1.0),
    (7004, 6377397.155, f64::NAN, 299.1528128, 1.0),
    (7005, 6377492.018, f64::NAN, 299.1528128, 1.0),
    (7007, 20926348.0, 20855233.0, f64::NAN, 0.304_797_265_4),
    (7008, 6378206.4, 6356583.8, f64::NAN, 1.0),
    (7010, 6378300.789, 6356566.435, f64::NAN, 1.0),
    (7011, 6378249.2, 6356515.0, f64::NAN, 1.0),
    (7012, 6378249.145, f64::NAN, 293.465, 1.0),
    (7013, 6378249.145, f64::NAN, 293.4663077, 1.0),
    (7014, 6378249.2, f64::NAN, 293.46598, 1.0),
    (7015, 6377276.345, f64::NAN, 300.8017, 1.0),
    (7016, 6377298.556, f64::NAN, 300.8017, 1.0),
    (7018, 6377304.063, f64::NAN, 300.8017, 1.0),
    (7019, 6378137.0, f64::NAN, 298.257222101, 1.0),
    (7020, 6378200.0, f64::NAN, 298.3, 1.0),
    (7021, 6378160.0, f64::NAN, 298.247, 1.0),
    (7022, 6378388.0, f64::NAN, 297.0, 1.0),
    (7024, 6378245.0, f64::NAN, 298.3, 1.0),
    (7025, 6378145.0, f64::NAN, 298.25, 1.0),
    (7027, 6376523.0, f64::NAN, 308.64, 1.0),
    (7028, 6378298.3, f64::NAN, 294.73, 1.0),
    (7029, 6378300.0, f64::NAN, 296.0, 1.0),
    (7030, 6378137.0, f64::NAN, 298.257223563, 1.0),
    (7031, 6378137.0, f64::NAN, 298.257223563, 1.0),
    (7032, 6378136.2, f64::NAN, 298.257223563, 1.0),
    (7033, 6378136.3, f64::NAN, 298.257223563, 1.0),
    (7034, 20926202.0, 20854895.0, f64::NAN, 0.304_797_265_4),
    (7036, 6378160.0, f64::NAN, 298.247167427, 1.0),
    (7041, 6378135.0, f64::NAN, 298.257, 1.0),
    (7042, 20922931.8, 20853374.58, f64::NAN, 12.0 / 39.370_142),
    (7043, 6378135.0, f64::NAN, 298.26, 1.0),
    (7044, 6377301.243, f64::NAN, 300.8017255, 1.0),
    (7045, 6377299.151, f64::NAN, 300.8017255, 1.0),
    (7046, 6377397.155, f64::NAN, 299.1528128, 1.000_013_596_5),
    (7048, 6371007.0, 6371007.0, f64::NAN, 1.0),
    (7049, 6378140.0, f64::NAN, 298.257, 1.0),
    (7050, 6378160.0, f64::NAN, 298.25, 1.0),
    (7051, 6377019.27, f64::NAN, 300.0, 1.0),
    (7052, 6370997.0, 6370997.0, f64::NAN, 1.0),
    (7053, 6378270.0, f64::NAN, 297.0, 1.0),
    (7054, 6378136.0, f64::NAN, 298.257839303, 1.0),
    (7055, 20926202.0, 20854895.0, f64::NAN, 0.304_8),
    (7056, 6377295.664, f64::NAN, 300.8017, 1.0),
    (7057, 6371228.0, 6371228.0, f64::NAN, 1.0),
    (7058, 6378273.0, 6356889.449, f64::NAN, 1.0),
    ];

    #[test]
    fn every_registry_ellipsoid_matches_its_definition() {
        for (code, a, b, rf, unit) in REGISTRY {
            let e = Ellipsoid::from_epsg(code).unwrap_or_else(|| panic!("{code} missing"));
            assert!((e.a() - a * unit).abs() <= 1e-15 * e.a(), "{code}: a {}", e.a());
            if rf.is_nan() {
                // Defined by its axes: the polar radius comes back to within a nanometre.
                assert!((e.b() - b * unit).abs() < 1e-8, "{code}: b {} vs {}", e.b(), b * unit);
            } else {
                assert!((1.0 / e.f() - rf).abs() <= 1e-13 * rf, "{code}: 1/f {}", 1.0 / e.f());
            }
        }
    }

    #[test]
    fn nothing_else_is_an_ellipsoid_code() {
        let found = (0..100_000).filter(|&c| Ellipsoid::from_epsg(c).is_some()).count();
        assert_eq!(found, REGISTRY.len());
    }

    #[test]
    fn named_ones_and_derived_values() {
        let named = [
            (7030, WGS84),
            (7019, GRS80),
            (1024, CGCS2000),
            (7054, PZ90),
            (1025, GSK_2011),
            (7024, KRASSOVSKY_1940),
            (7049, IAG_1975),
            (7043, WGS72),
            (7025, WGS66),
            (7036, GRS67),
            (7003, AUSTRALIAN_NATIONAL),
            (7050, AUSTRALIAN_NATIONAL),
            (7022, INTERNATIONAL_1924),
            (7001, AIRY_1830),
            (7004, BESSEL_1841),
            (7008, CLARKE_1866),
            (7015, EVEREST_1830),
            (7018, EVEREST_1830_MODIFIED),
        ];
        for (code, e) in named {
            assert_eq!(Ellipsoid::from_epsg(code), Some(e), "{code}");
        }
        // Values EPSG and PROJ quote for the ones it defines in other units.
        let everest = Ellipsoid::from_epsg(7042).unwrap(); // Indian feet
        assert!((everest.a() - 6_377_299.365_595_4).abs() < 1e-6, "{}", everest.a());
        assert!(
            (1.0 / everest.f() - 300.801_725_5).abs() < 1e-6,
            "{}",
            1.0 / everest.f()
        );
        let namibia = Ellipsoid::from_epsg(7046).unwrap(); // German legal metres
        assert!((namibia.a() - 6_377_483.865).abs() < 1e-3, "{}", namibia.a());
        let clarke_1880 = Ellipsoid::from_epsg(7034).unwrap(); // Clarke's feet
        assert!((1.0 / clarke_1880.f() - 293.466_307_656).abs() < 1e-6);
        // Spheres have no flattening.
        assert_eq!(Ellipsoid::from_epsg(7048).unwrap().f(), 0.0);
        // WGS 60 isn't in EPSG; its figure is the US Department of Defense's.
        assert_eq!((WGS60.a(), 1.0 / WGS60.f()), (6_378_165.0, 298.3));
    }
}
