//! Optional conversions between the coordinate types and the common math / GIS crates.
//!
//! Each block is behind its own default-off feature, so the crate stays dependency-free unless you
//! ask for a bridge:
//!
//! - **`mint`** — [`From`]/[`Into`] with `mint::Vector3<f64>` (all vector types) and `mint::Point3<f64>`
//!   (the position types). `mint` is the Math-Interoperability-Types lingua franca, so this one bridge
//!   also reaches `glam`, `nalgebra` and `cgmath`, which all speak `mint`.
//! - **`glam`** — with `glam::DVec3`.
//! - **`nalgebra`** — with `nalgebra::Vector3<f64>` and `nalgebra::Point3<f64>`.
//! - **`geo-types`** — [`Geodetic`](crate::Geodetic) ↔ `geo_types::Point<f64>` / `Coord<f64>` for 2D
//!   GIS pipelines. **Longitude is `x`, latitude is `y`** (the `geo` convention) and the ellipsoidal
//!   height is dropped on the way out and taken as `0` on the way in.
//!
//! ```
//! # #[cfg(feature = "glam")] {
//! use geo3d::Ecef;
//! let e = Ecef::new(1.0, 2.0, 3.0);
//! let v: glam::DVec3 = e.into();
//! assert_eq!(Ecef::from(v), e);
//! # }
//! ```

// The macros below only expand under a feature, so they are otherwise-unused without one.
#![allow(unused_macros)]

/// `From`/`Into` both ways between a geo3d 3-component type and a foreign 3-vector type. The forward
/// direction reads the named fields; the reverse always builds through `Type::new(x, y, z)`.
macro_rules! vec3_bridge {
    ($geo:ident { $a:ident, $b:ident, $c:ident } <=> $foreign:ty : $ctor:expr, |$v:ident| $rev:expr) => {
        impl From<crate::$geo> for $foreign {
            #[inline]
            fn from(g: crate::$geo) -> Self {
                let ctor = $ctor;
                ctor(g.$a, g.$b, g.$c)
            }
        }
        impl From<$foreign> for crate::$geo {
            #[inline]
            fn from($v: $foreign) -> Self {
                let (x, y, z) = $rev;
                crate::$geo::new(x, y, z)
            }
        }
    };
}

#[cfg(feature = "mint")]
mod mint_impls {
    vec3_bridge!(Vec3  { x, y, z }            <=> mint::Vector3<f64> : |x, y, z| mint::Vector3 { x, y, z }, |v| (v.x, v.y, v.z));
    vec3_bridge!(Enu   { east, north, up }    <=> mint::Vector3<f64> : |x, y, z| mint::Vector3 { x, y, z }, |v| (v.x, v.y, v.z));
    vec3_bridge!(Ned   { north, east, down }  <=> mint::Vector3<f64> : |x, y, z| mint::Vector3 { x, y, z }, |v| (v.x, v.y, v.z));
    vec3_bridge!(Ecef  { x, y, z }            <=> mint::Vector3<f64> : |x, y, z| mint::Vector3 { x, y, z }, |v| (v.x, v.y, v.z));
    vec3_bridge!(Eci   { x, y, z }            <=> mint::Vector3<f64> : |x, y, z| mint::Vector3 { x, y, z }, |v| (v.x, v.y, v.z));
    // The position types also convert to/from `Point3`.
    vec3_bridge!(Ecef  { x, y, z }            <=> mint::Point3<f64>  : |x, y, z| mint::Point3 { x, y, z },  |v| (v.x, v.y, v.z));
    vec3_bridge!(Eci   { x, y, z }            <=> mint::Point3<f64>  : |x, y, z| mint::Point3 { x, y, z },  |v| (v.x, v.y, v.z));
}

#[cfg(feature = "glam")]
mod glam_impls {
    vec3_bridge!(Vec3  { x, y, z }            <=> glam::DVec3 : |x, y, z| glam::DVec3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Enu   { east, north, up }    <=> glam::DVec3 : |x, y, z| glam::DVec3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Ned   { north, east, down }  <=> glam::DVec3 : |x, y, z| glam::DVec3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Ecef  { x, y, z }            <=> glam::DVec3 : |x, y, z| glam::DVec3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Eci   { x, y, z }            <=> glam::DVec3 : |x, y, z| glam::DVec3::new(x, y, z), |v| (v.x, v.y, v.z));
}

#[cfg(feature = "nalgebra")]
mod nalgebra_impls {
    vec3_bridge!(Vec3  { x, y, z }            <=> nalgebra::Vector3<f64> : |x, y, z| nalgebra::Vector3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Enu   { east, north, up }    <=> nalgebra::Vector3<f64> : |x, y, z| nalgebra::Vector3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Ned   { north, east, down }  <=> nalgebra::Vector3<f64> : |x, y, z| nalgebra::Vector3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Ecef  { x, y, z }            <=> nalgebra::Vector3<f64> : |x, y, z| nalgebra::Vector3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Eci   { x, y, z }            <=> nalgebra::Vector3<f64> : |x, y, z| nalgebra::Vector3::new(x, y, z), |v| (v.x, v.y, v.z));
    vec3_bridge!(Ecef  { x, y, z }            <=> nalgebra::Point3<f64>  : |x, y, z| nalgebra::Point3::new(x, y, z),  |v| (v.coords.x, v.coords.y, v.coords.z));
    vec3_bridge!(Eci   { x, y, z }            <=> nalgebra::Point3<f64>  : |x, y, z| nalgebra::Point3::new(x, y, z),  |v| (v.coords.x, v.coords.y, v.coords.z));
}

#[cfg(feature = "geo-types")]
mod geo_types_impls {
    use crate::Geodetic;

    // The `geo` convention is (x = longitude, y = latitude). Height is dropped going out, 0 coming in.
    impl From<Geodetic> for geo_types::Point<f64> {
        #[inline]
        fn from(g: Geodetic) -> Self {
            geo_types::Point::new(g.lon_deg, g.lat_deg)
        }
    }
    impl From<Geodetic> for geo_types::Coord<f64> {
        #[inline]
        fn from(g: Geodetic) -> Self {
            geo_types::Coord {
                x: g.lon_deg,
                y: g.lat_deg,
            }
        }
    }
    impl From<geo_types::Point<f64>> for Geodetic {
        #[inline]
        fn from(p: geo_types::Point<f64>) -> Self {
            Geodetic::new(p.y(), p.x(), 0.0)
        }
    }
    impl From<geo_types::Coord<f64>> for Geodetic {
        #[inline]
        fn from(c: geo_types::Coord<f64>) -> Self {
            Geodetic::new(c.y, c.x, 0.0)
        }
    }
}

#[cfg(all(
    test,
    any(feature = "mint", feature = "glam", feature = "nalgebra", feature = "geo-types")
))]
mod tests {
    #[allow(unused_imports)]
    use crate::{Ecef, Eci, Enu, Geodetic, Ned, Vec3};

    #[cfg(feature = "mint")]
    #[test]
    fn mint_round_trips_and_orders_fields() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let m: mint::Vector3<f64> = v.into();
        assert_eq!((m.x, m.y, m.z), (1.0, 2.0, 3.0));
        assert_eq!(Vec3::from(m), v);
        // ENU field order maps east→x, north→y, up→z, and round-trips.
        let enu = Enu::new(4.0, 5.0, 6.0);
        let me: mint::Vector3<f64> = enu.into();
        assert_eq!((me.x, me.y, me.z), (4.0, 5.0, 6.0));
        assert_eq!(Enu::from(me), enu);
        // Ned maps north→x, east→y, down→z.
        let ned = Ned::new(7.0, 8.0, 9.0);
        let mn: mint::Vector3<f64> = ned.into();
        assert_eq!((mn.x, mn.y, mn.z), (7.0, 8.0, 9.0));
        assert_eq!(Ned::from(mn), ned);
        // Positions also bridge to Point3.
        let e = Ecef::new(10.0, 20.0, 30.0);
        let p: mint::Point3<f64> = e.into();
        assert_eq!(Ecef::from(p), e);
        let i = Eci::new(-1.0, -2.0, -3.0);
        let pi: mint::Point3<f64> = i.into();
        assert_eq!(Eci::from(pi), i);
    }

    #[cfg(feature = "glam")]
    #[test]
    fn glam_round_trips() {
        let e = Ecef::new(1.5, -2.5, 3.5);
        let g: glam::DVec3 = e.into();
        assert_eq!((g.x, g.y, g.z), (1.5, -2.5, 3.5));
        assert_eq!(Ecef::from(g), e);
        let enu = Enu::new(1.0, 2.0, 3.0);
        assert_eq!(Enu::from(glam::DVec3::from(enu)), enu);
    }

    #[cfg(feature = "nalgebra")]
    #[test]
    fn nalgebra_round_trips_vector_and_point() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let nv: nalgebra::Vector3<f64> = v.into();
        assert_eq!((nv.x, nv.y, nv.z), (1.0, 2.0, 3.0));
        assert_eq!(Vec3::from(nv), v);
        let e = Ecef::new(10.0, 20.0, 30.0);
        let np: nalgebra::Point3<f64> = e.into();
        assert_eq!(Ecef::from(np), e);
    }

    #[cfg(feature = "geo-types")]
    #[test]
    fn geo_types_uses_lon_x_lat_y_and_drops_height() {
        let g = Geodetic::new(-33.87, 151.21, 20.0); // lat, lon, height
        let p: geo_types::Point<f64> = g.into();
        assert_eq!((p.x(), p.y()), (151.21, -33.87)); // x = lon, y = lat
        let back = Geodetic::from(p);
        assert_eq!((back.lat_deg, back.lon_deg, back.height_m), (-33.87, 151.21, 0.0));
        let c: geo_types::Coord<f64> = g.into();
        assert_eq!((c.x, c.y), (151.21, -33.87));
        assert_eq!(Geodetic::from(c).lat_deg, -33.87);
    }
}
