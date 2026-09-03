//! 3D geospatial core (→ the CesiumJS 3D app).
//!
//! Our 2D stack (the `geo` crate, geodatafusion, OpenSearch `geo_point`) is 2D and stays that way. True 3D
//! geometry cannot be done on geodetic `(lat, lon, alt)` directly — it is NOT a Euclidean space — so this
//! crate provides the **metric-frame foundation**: WGS84 geodetic → **ECEF** (Earth-Centered Earth-Fixed,
//! metres — the frame CesiumJS renders in) and a local **ENU** (East-North-Up) tangent frame, then
//! closed-form volumetric primitives on top (slant range, extruded prism, sphere, cone). All pure `f64`,
//! zero deps. Higher tiers (line-of-sight over terrain, arbitrary meshes) can later be backed by a 3D lib
//! without changing these signatures.
//!
//! DOMAIN-NEUTRAL: these are generic geometry — an altitude band is aviation airspace OR a warehouse
//! mezzanine range; a cone is a camera field-of-view OR an antenna beam OR a spotlight; a sphere is a
//! delivery radius OR a proximity alert; closest-approach is aircraft separation OR two ships' CPA. No
//! primitive here is specific to any one domain.

/// WGS84 geodetic point: latitude/longitude in **degrees**, altitude in **metres** above the ellipsoid.
pub type Geo = (f64, f64, f64);
/// A Cartesian vector in metres (ECEF, or ENU relative to a reference).
pub type Vec3 = [f64; 3];

// ── WGS84 ellipsoid constants ──────────────────────────────────────────────────────────────────────
const WGS84_A: f64 = 6_378_137.0; // semi-major axis (equatorial radius), metres
const WGS84_F: f64 = 1.0 / 298.257_223_563; // flattening
const WGS84_E2: f64 = WGS84_F * (2.0 - WGS84_F); // first eccentricity squared

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: Vec3) -> f64 {
    dot(a, a).sqrt()
}

/// WGS84 geodetic → ECEF Cartesian (metres). The exact, standard closed form. This is CesiumJS's native
/// frame, so a track's `(lat, lon, altitudeM)` maps straight to what the globe draws.
pub fn geodetic_to_ecef((lat_deg, lon_deg, alt_m): Geo) -> Vec3 {
    let (lat, lon) = (lat_deg.to_radians(), lon_deg.to_radians());
    let (slat, clat) = lat.sin_cos();
    let (slon, clon) = lon.sin_cos();
    let n = WGS84_A / (1.0 - WGS84_E2 * slat * slat).sqrt(); // prime-vertical radius of curvature
    [
        (n + alt_m) * clat * clon,
        (n + alt_m) * clat * slon,
        (n * (1.0 - WGS84_E2) + alt_m) * slat,
    ]
}

/// Rotate an ECEF point into the local **ENU** (East, North, Up) frame anchored at `reference`. Distances
/// are metres; the frame is right-handed with Up = the ellipsoid normal. Used by directional primitives
/// (cones, look angles) where "az/el from a reference point" is the natural description.
pub fn ecef_to_enu(reference: Geo, target_ecef: Vec3) -> Vec3 {
    let (lat, lon) = (reference.0.to_radians(), reference.1.to_radians());
    let (slat, clat) = lat.sin_cos();
    let (slon, clon) = lon.sin_cos();
    let d = sub(target_ecef, geodetic_to_ecef(reference));
    [
        dot([-slon, clon, 0.0], d),               // East
        dot([-slat * clon, -slat * slon, clat], d), // North
        dot([clat * clon, clat * slon, slat], d),   // Up
    ]
}

/// The local ENU offset (metres) of `point` as seen from `reference`.
pub fn enu_offset(reference: Geo, point: Geo) -> Vec3 {
    ecef_to_enu(reference, geodetic_to_ecef(point))
}

/// True 3D **slant range** between two geodetic points, in metres (straight-line ECEF chord — the real
/// through-space distance between two things at altitude, unlike the 2D great-circle `haversineKm`).
pub fn distance3d_m(a: Geo, b: Geo) -> f64 {
    norm(sub(geodetic_to_ecef(a), geodetic_to_ecef(b)))
}

/// Altitude-band test (metres). `None` bound = unbounded on that side, so a floor-only or ceiling-only
/// volume works. The vertical half of an extruded-prism volume.
pub fn altitude_between(alt_m: f64, floor_m: Option<f64>, ceiling_m: Option<f64>) -> bool {
    alt_m >= floor_m.unwrap_or(f64::NEG_INFINITY) && alt_m <= ceiling_m.unwrap_or(f64::INFINITY)
}

/// Whether `point` lies inside a **sphere** of `radius_m` about `center` — a spherical range/coverage test
/// (a service radius, a proximity alert, an exclusion bubble, …).
pub fn within_sphere(center: Geo, radius_m: f64, point: Geo) -> bool {
    distance3d_m(center, point) <= radius_m
}

/// A right-circular **cone** — a directional 3D volume: a coverage sector, a field-of-view, a beam or
/// spotlight, a broadcast/emission cone, a scan arc, etc. Apex at a reference point, pointing along a
/// compass `azimuth_deg` (0 = North, clockwise) and `elevation_deg` (0 = horizontal, +up), with a
/// `half_angle_deg` spread and a slant `range_m`.
pub struct Cone {
    pub apex: Geo,
    pub azimuth_deg: f64,
    pub elevation_deg: f64,
    pub half_angle_deg: f64,
    pub range_m: f64,
}

impl Cone {
    /// Unit axis direction in the apex's ENU frame.
    fn axis_enu(&self) -> Vec3 {
        let (az, el) = (self.azimuth_deg.to_radians(), self.elevation_deg.to_radians());
        let (saz, caz) = az.sin_cos();
        let cel = el.cos();
        [saz * cel, caz * cel, el.sin()] // East, North, Up
    }

    /// Whether `point` is inside the cone: within `range_m` of the apex AND within `half_angle` of the axis.
    pub fn contains(&self, point: Geo) -> bool {
        let p = enu_offset(self.apex, point);
        let r = norm(p);
        if r > self.range_m {
            return false;
        }
        if r == 0.0 {
            return true; // the apex itself
        }
        let cos_angle = dot(p, self.axis_enu()) / r; // axis is unit length
        cos_angle >= self.half_angle_deg.to_radians().cos()
    }
}

/// An **extruded prism** volume: a 2D footprint (owned/tested by the caller via the 2D `geo` stack) plus a
/// vertical band. This function is the vertical half — compose it with a 2D containment test for the full
/// "inside this airspace volume". Kept here so the band semantics live with the rest of the 3D core.
pub fn within_prism(footprint_contains_2d: bool, alt_m: f64, floor_m: Option<f64>, ceiling_m: Option<f64>) -> bool {
    footprint_contains_2d && altitude_between(alt_m, floor_m, ceiling_m)
}

// ── ECEF ↔ geodetic (the inverse) + ray/kinematic/datum utilities ────────────────────────────────────

/// ECEF (metres) → WGS84 geodetic `(lat°, lon°, alt_m)` — Bowring's closed-form method (sub-mm). The
/// inverse of [`geodetic_to_ecef`]; required to turn any computed ECEF point (e.g. where a pointing ray
/// meets the surface) back into a lat/lon a map/store can use.
pub fn ecef_to_geodetic([x, y, z]: Vec3) -> Geo {
    let b = WGS84_A * (1.0 - WGS84_F); // polar radius
    let ep2 = (WGS84_A * WGS84_A - b * b) / (b * b); // second eccentricity²
    let p = (x * x + y * y).sqrt();
    let lon = y.atan2(x);
    if p < 1e-9 {
        let lat = if z >= 0.0 { std::f64::consts::FRAC_PI_2 } else { -std::f64::consts::FRAC_PI_2 };
        return (lat.to_degrees(), lon.to_degrees(), z.abs() - b);
    }
    let theta = (z * WGS84_A).atan2(p * b);
    let (st, ct) = theta.sin_cos();
    let lat = (z + ep2 * b * st * st * st).atan2(p - WGS84_E2 * WGS84_A * ct * ct * ct);
    let (slat, clat) = lat.sin_cos();
    let n = WGS84_A / (1.0 - WGS84_E2 * slat * slat).sqrt();
    (lat.to_degrees(), lon.to_degrees(), p / clat - n)
}

/// Rotate a local ENU vector (East, North, Up) at `reference` into an ECEF vector — the inverse rotation
/// of [`ecef_to_enu`] (translation dropped, so it suits ray DIRECTIONS). e/n/u are orthonormal.
pub fn enu_to_ecef_dir(reference: Geo, enu: Vec3) -> Vec3 {
    let (lat, lon) = (reference.0.to_radians(), reference.1.to_radians());
    let (slat, clat) = lat.sin_cos();
    let (slon, clon) = lon.sin_cos();
    let e = [-slon, clon, 0.0];
    let n = [-slat * clon, -slat * slon, clat];
    let u = [clat * clon, clat * slon, slat];
    [
        enu[0] * e[0] + enu[1] * n[0] + enu[2] * u[0],
        enu[0] * e[1] + enu[1] * n[1] + enu[2] * u[1],
        enu[0] * e[2] + enu[1] * n[2] + enu[2] * u[2],
    ]
}

/// Nearest ECEF intersection of the ray `origin + t·dir` (t ≥ 0) with the WGS84 ellipsoid, or `None` if
/// it misses (points into space). Scales space to a unit sphere (÷a, ÷a, ÷b) and solves the ray-sphere
/// quadratic. `dir` need not be unit-length. The oblate-spheroid version of "where does this ray hit the
/// Earth" — line-of-sight ground point, a camera/antenna footprint, a projected bearing, etc.
pub fn ray_ellipsoid_intersect(origin: Vec3, dir: Vec3) -> Option<Vec3> {
    let b = WGS84_A * (1.0 - WGS84_F);
    let s = [1.0 / WGS84_A, 1.0 / WGS84_A, 1.0 / b];
    let o = [origin[0] * s[0], origin[1] * s[1], origin[2] * s[2]];
    let d = [dir[0] * s[0], dir[1] * s[1], dir[2] * s[2]];
    let (qa, qb, qc) = (dot(d, d), 2.0 * dot(o, d), dot(o, o) - 1.0);
    let disc = qb * qb - 4.0 * qa * qc;
    if disc < 0.0 || qa == 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let (t1, t2) = ((-qb - sq) / (2.0 * qa), (-qb + sq) / (2.0 * qa));
    let t = if t1 >= 0.0 {
        t1
    } else if t2 >= 0.0 {
        t2
    } else {
        return None; // ellipsoid is behind the ray
    };
    Some([origin[0] + t * dir[0], origin[1] + t * dir[1], origin[2] + t * dir[2]])
}

/// Where a ray from `observer`, pointing along compass `azimuth_deg` (0 = N, CW) and `elevation_deg`
/// (0 = horizon, +up), strikes the WGS84 surface — the geodetic ground point, or `None` if it points above
/// the horizon into space. Composes ENU→ECEF pointing with the ray-ellipsoid solve and the inverse.
pub fn ground_strike(observer: Geo, azimuth_deg: f64, elevation_deg: f64) -> Option<Geo> {
    let (az, el) = (azimuth_deg.to_radians(), elevation_deg.to_radians());
    let cel = el.cos();
    let dir = enu_to_ecef_dir(observer, [az.sin() * cel, az.cos() * cel, el.sin()]);
    ray_ellipsoid_intersect(geodetic_to_ecef(observer), dir).map(ecef_to_geodetic)
}

/// Azimuth (° from N, CW), elevation (° above horizon), and range (metres) from `observer` to `target` —
/// the look angles to point toward a target (the inverse of the cone axis).
pub fn look_angles(observer: Geo, target: Geo) -> (f64, f64, f64) {
    let [e, n, u] = enu_offset(observer, target);
    let range = (e * e + n * n + u * u).sqrt();
    let az = e.atan2(n).to_degrees().rem_euclid(360.0);
    let el = if range == 0.0 { 0.0 } else { (u / range).asin().to_degrees() };
    (az, el, range)
}

/// Closest point of approach of two entities in LINEAR motion — ECEF positions (m) + velocities (m/s).
/// Returns `(t_seconds, min_distance_m)`. Assumes constant velocity (a first-order conflict/collision
/// estimate); `t` may be negative when the closest approach is already in the past. Feed `min_distance_m`
/// to a threshold for a proximity/separation/conflict check (aircraft separation, vessel CPA, rendezvous).
pub fn closest_approach(p1: Vec3, v1: Vec3, p2: Vec3, v2: Vec3) -> (f64, f64) {
    let prel = sub(p1, p2);
    let vrel = sub(v1, v2);
    let vv = dot(vrel, vrel);
    let t = if vv == 0.0 { 0.0 } else { -dot(prel, vrel) / vv };
    let at = [prel[0] + vrel[0] * t, prel[1] + vrel[1] * t, prel[2] + vrel[2] * t];
    (t, norm(at))
}

/// 7-parameter Helmert (Bursa-Wolf) datum-shift parameters: translations (metres), rotations (RADIANS),
/// and scale (unitless, e.g. ppm × 1e-6). Values are per-datum — look them up in EPSG (e.g. a local grid
/// → WGS84 set); deliberately NOT hardcoded here (a wrong constant is worse than none).
#[derive(Clone, Copy, Debug, Default)]
pub struct Helmert7 {
    pub tx: f64,
    pub ty: f64,
    pub tz: f64,
    pub rx: f64,
    pub ry: f64,
    pub rz: f64,
    pub scale: f64,
}

/// Apply a Helmert transform to an ECEF point (source datum → WGS84). Pure matrix mult over `Vec3`; the
/// approximate inverse is the same call with every parameter negated (valid for the usual small angles).
pub fn helmert(p: Vec3, h: &Helmert7) -> Vec3 {
    let m = 1.0 + h.scale;
    [
        h.tx + m * (p[0] - h.rz * p[1] + h.ry * p[2]),
        h.ty + m * (h.rz * p[0] + p[1] - h.rx * p[2]),
        h.tz + m * (-h.ry * p[0] + h.rx * p[1] + p[2]),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn ecef_hits_the_known_axis_points() {
        // Equator/prime-meridian at the surface → (a, 0, 0).
        let e = geodetic_to_ecef((0.0, 0.0, 0.0));
        assert!(close(e[0], WGS84_A, 1e-3) && close(e[1], 0.0, 1e-6) && close(e[2], 0.0, 1e-6));
        // North pole → (0, 0, polar radius b = a(1-f)).
        let n = geodetic_to_ecef((90.0, 0.0, 0.0));
        assert!(close(n[0], 0.0, 1e-3) && close(n[1], 0.0, 1e-3) && close(n[2], WGS84_A * (1.0 - WGS84_F), 1.0));
        // Altitude adds along the outward normal at the equator.
        let up = geodetic_to_ecef((0.0, 0.0, 1000.0));
        assert!(close(up[0], WGS84_A + 1000.0, 1e-3));
    }

    #[test]
    fn distance3d_is_vertical_then_horizontal() {
        // Pure vertical separation ≈ the altitude difference.
        assert!(close(distance3d_m((0.0, 0.0, 0.0), (0.0, 0.0, 1000.0)), 1000.0, 1e-3));
        // 1° of longitude at the equator is a ~111.3 km chord (slightly under the arc length).
        let d = distance3d_m((0.0, 0.0, 0.0), (0.0, 1.0, 0.0));
        assert!(d > 111_000.0 && d < 111_400.0, "got {d}");
    }

    #[test]
    fn enu_orients_east_north_up() {
        // A point directly above the reference → straight Up.
        let u = enu_offset((10.0, 20.0, 0.0), (10.0, 20.0, 500.0));
        assert!(close(u[0], 0.0, 1e-3) && close(u[1], 0.0, 1e-3) && close(u[2], 500.0, 1e-2));
        // A point slightly north → positive North, ~zero East.
        let n = enu_offset((10.0, 20.0, 0.0), (10.01, 20.0, 0.0));
        assert!(n[1] > 0.0 && close(n[0], 0.0, 1.0));
    }

    #[test]
    fn altitude_band_and_sphere() {
        assert!(altitude_between(5000.0, Some(3000.0), Some(10000.0)));
        assert!(!altitude_between(2000.0, Some(3000.0), Some(10000.0)));
        assert!(altitude_between(999_999.0, None, None)); // unbounded
        assert!(altitude_between(1000.0, None, Some(2000.0))); // ceiling only
        // A 10 km sphere about a site: a point 1 km up is inside, 20 km away is outside.
        let site = (34.0, -118.0, 0.0);
        assert!(within_sphere(site, 10_000.0, (34.0, -118.0, 1000.0)));
        assert!(!within_sphere(site, 10_000.0, (34.2, -118.0, 0.0)));
    }

    #[test]
    fn cone_directional_volume() {
        // A cone pointing straight up (elevation 90°), 30° half-angle, 20 km range, from a ground site.
        let c = Cone { apex: (0.0, 0.0, 0.0), azimuth_deg: 0.0, elevation_deg: 90.0, half_angle_deg: 30.0, range_m: 20_000.0 };
        assert!(c.contains((0.0, 0.0, 5000.0))); // straight up, in range
        assert!(!c.contains((0.0, 0.0, 25_000.0))); // straight up but beyond range
        assert!(!c.contains((0.2, 0.0, 1000.0))); // low + off to the side → outside the 30° spread
        // A cone aimed due east at the horizon includes a point to the east, not one to the north.
        let east = Cone { apex: (0.0, 0.0, 1000.0), azimuth_deg: 90.0, elevation_deg: 0.0, half_angle_deg: 20.0, range_m: 100_000.0 };
        assert!(east.contains((0.0, 0.3, 1000.0)));
        assert!(!east.contains((0.3, 0.0, 1000.0)));
    }

    #[test]
    fn prism_composes_footprint_and_band() {
        assert!(within_prism(true, 8000.0, Some(5000.0), Some(12000.0)));
        assert!(!within_prism(false, 8000.0, Some(5000.0), Some(12000.0))); // outside footprint
        assert!(!within_prism(true, 3000.0, Some(5000.0), Some(12000.0))); // below the floor
    }

    #[test]
    fn ecef_to_geodetic_round_trips() {
        for g in [(0.0, 0.0, 0.0), (34.05, -118.24, 100.0), (-33.87, 151.21, 2000.0), (89.9, 10.0, 500.0)] {
            let (lat, lon, alt) = ecef_to_geodetic(geodetic_to_ecef(g));
            assert!(close(lat, g.0, 1e-8) && close(lon, g.1, 1e-8) && close(alt, g.2, 1e-4), "{g:?} -> {:?}", (lat, lon, alt));
        }
    }

    #[test]
    fn ray_ellipsoid_and_ground_strike() {
        // Straight down from 10 km over the equator/prime meridian hits ~(0,0,0).
        let g = ground_strike((0.0, 0.0, 10_000.0), 0.0, -90.0).unwrap();
        assert!(close(g.0, 0.0, 1e-6) && close(g.1, 0.0, 1e-6) && close(g.2, 0.0, 1e-3));
        // Pointing UP into space → no intersection.
        assert!(ground_strike((0.0, 0.0, 10_000.0), 0.0, 90.0).is_none());
        // A shallow look east from altitude strikes the surface somewhere east + below.
        let east = ground_strike((0.0, 0.0, 10_000.0), 90.0, -5.0).unwrap();
        assert!(east.1 > 0.0 && east.2.abs() < 1.0);
    }

    #[test]
    fn look_angles_point_toward_a_target() {
        // Target due north → azimuth ~0.
        let (az, _el, r) = look_angles((0.0, 0.0, 0.0), (0.1, 0.0, 0.0));
        assert!((az < 1.0 || az > 359.0) && r > 0.0);
        // Target due east → azimuth ~90.
        let (az, ..) = look_angles((0.0, 0.0, 0.0), (0.0, 0.1, 0.0));
        assert!(close(az, 90.0, 1.0));
        // Target straight up → elevation ~90.
        let (_az, el, r) = look_angles((10.0, 20.0, 0.0), (10.0, 20.0, 5000.0));
        assert!(close(el, 90.0, 1e-2) && close(r, 5000.0, 1e-2));
    }

    #[test]
    fn closest_approach_time_and_distance() {
        // Two objects converging on the x-axis, offset 5 m on y; CPA at t=5 s, min distance = 5 m.
        let (t, d) = closest_approach([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [10.0, 5.0, 0.0], [-1.0, 0.0, 0.0]);
        assert!(close(t, 5.0, 1e-9) && close(d, 5.0, 1e-9));
        // Parallel, non-closing → t=0, distance = current separation.
        let (t, d) = closest_approach([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(close(t, 0.0, 1e-9) && close(d, 3.0, 1e-9));
    }

    #[test]
    fn helmert_identity_and_translation() {
        let p = [1000.0, 2000.0, 3000.0];
        assert_eq!(helmert(p, &Helmert7::default()), p); // no params = identity
        let shifted = helmert(p, &Helmert7 { tx: 10.0, ty: -5.0, tz: 2.0, ..Default::default() });
        assert_eq!(shifted, [1010.0, 1995.0, 3002.0]); // pure translation
    }
}
