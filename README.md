# geo3d

[![crates.io](https://img.shields.io/crates/v/geo3d.svg)](https://crates.io/crates/geo3d)
[![docs.rs](https://docs.rs/geo3d/badge.svg)](https://docs.rs/geo3d)
[![ci](https://github.com/FusionVault/geo3d/actions/workflows/ci.yml/badge.svg)](https://github.com/FusionVault/geo3d/actions/workflows/ci.yml)

3D geodesy in pure `f64`, zero dependencies.

Latitude, longitude and height are not a Euclidean space, so you cannot subtract two of them and
get a distance, an angle or a direction. `geo3d` gives you the metric frames that make 3D geometry
work — Earth-centred (ECEF) and local East-North-Up — and the everyday calculations on top of them:
slant range, look angles, where a line of sight meets the ground, whether a point is inside a
volume, when two moving objects get closest, geodesic distance along the ellipsoid.

Everything is typed (`Geodetic`, `Ecef`, `Enu`, `Aer`, …), `Copy`, allocation-free and
`#![forbid(unsafe_code)]`. Angles in the API are degrees, lengths are metres. ECEF is the frame
CesiumJS renders in, so a `(lat, lon, height)` track maps straight onto a globe.

## Install

```bash
cargo add geo3d
```

Optional: `cargo add geo3d --features serde` derives `Serialize`/`Deserialize` on every coordinate type.

## Quick start

```rust
use geo3d::{distance3d_m, Geodetic, LocalFrame, Sphere, Volume};

let site = Geodetic::new(34.05, -118.24, 100.0);      // latitude°, longitude°, height m
let aircraft = Geodetic::new(34.30, -118.10, 3_000.0);

let slant = distance3d_m(site, aircraft);              // straight-line metres through space
let aer = LocalFrame::new(site).look_angles(aircraft); // azimuth°, elevation°, range m
let near = Sphere::new(site, 40_000.0).contains(aircraft);

assert!(slant > 30_000.0 && slant < 32_000.0);
assert!(aer.elevation_deg > 4.0 && near);
```

Every example in this README is compiled and run as a doctest by CI, so they are safe to copy.

## Guide

### Coordinates and ellipsoids

`Geodetic` is latitude/longitude in degrees and ellipsoidal height in metres. `Ecef` is metres from
the Earth's centre. An `Ellipsoid` is a value — `WGS84` is the default and has shorthand free
functions; `GRS80`, `WGS72`, `AIRY_1830`, `BESSEL_1841`, `CLARKE_1866`, `INTERNATIONAL_1924` and
`PZ90` are there too, or build your own with `Ellipsoid::new(a, f)`.

```rust
use geo3d::{ecef_to_geodetic, geodetic_to_ecef, Ecef, Geodetic, GRS80, WGS84};

let sydney = Geodetic::new(-33.87, 151.21, 20.0);

let e: Ecef = geodetic_to_ecef(sydney);       // WGS84 shorthand
let e2 = WGS84.to_ecef(sydney);               // the same, spelled out
let back = ecef_to_geodetic(e);
assert_eq!(e, e2);
assert!((back.lat_deg - sydney.lat_deg).abs() < 1e-11);

let on_grs80 = GRS80.to_ecef(sydney);         // any ellipsoid works the same way
assert!((on_grs80.x - e.x).abs() < 1e-3);
```

### Local frames: ENU and NED

A `LocalFrame` is a tangent plane anchored at a reference point. It precomputes its rotation once,
so converting many targets against one observer is cheap. Every conversion exists in both the
East-North-Up and North-East-Down conventions.

```rust
use geo3d::{Enu, Geodetic, LocalFrame};

let radar = LocalFrame::new(Geodetic::new(51.5, -0.12, 30.0));
let target = Geodetic::new(51.6, -0.12, 3_000.0);   // ~11 km due north, 3 km up

let enu = radar.enu_of(target);                     // east, north, up in metres
let ned = radar.ned_of(target);                     // north, east, down
assert!(enu.east.abs() < 1.0 && enu.north > 11_000.0 && ned.down < 0.0);

// ...and back: an offset in the frame as a geodetic position.
let there = radar.geodetic_of(Enu::new(0.0, 11_000.0, 2_970.0));
assert!(there.lat_deg > 51.59);
```

Use `to_enu` / `to_ecef` when you already have ECEF positions, and `dir_to_enu` / `dir_to_ecef`
for directions (rotation only, no translation).

### Look angles, pointing and ground strike

`Aer` is a pointing from an observer: azimuth (degrees clockwise from north), elevation (degrees
above the horizon) and range. `look_angles` gives it, `destination` inverts it, and `ground_strike`
finds where a pointing meets the ellipsoid surface.

```rust
use geo3d::{Aer, Geodetic, LocalFrame};

let station = LocalFrame::new(Geodetic::new(-33.94, 151.18, 21.0));
let aircraft = Geodetic::new(-33.80, 151.30, 2_500.0);

let aer = station.look_angles(aircraft);          // ≈ az 35.6°, el 7.3°, range 19 254 m
assert!(aer.azimuth_deg > 35.0 && aer.azimuth_deg < 36.0 && aer.is_above_horizon());

let again = station.destination(aer);             // the aircraft, from the pointing
assert!((again.lat_deg - aircraft.lat_deg).abs() < 1e-9);

// A camera looking 2° below the horizon along that bearing: where does it hit the ground?
let ground = station.ground_strike(aer.azimuth_deg, -2.0);
assert!(ground.is_some());
assert!(station.ground_strike(0.0, 10.0).is_none()); // above the horizon: never meets the surface

// Pointings convert to local vectors and back.
let enu = Aer::new(90.0, 0.0, 1_000.0).to_enu();  // 1 km due east
assert!((enu.east - 1_000.0).abs() < 1e-9);
```

### Distances: straight line, great circle, geodesic

Three different questions, three functions. Pick by what you mean:

- `distance3d_m` — the straight line through space (the ECEF chord). Right for anything at height:
  aircraft separation, sensor range, satellite slant range.
- `great_circle_distance_m` — along the surface on a mean-radius sphere. Fast, within about 0.3 %
  of the true value; fine for sorting, thresholds and display.
- `geodesic_inverse` — along the surface of the ellipsoid (Vincenty), accurate to about 0.5 mm.
  Returns `None` only for nearly antipodal points.

```rust
use geo3d::{distance3d_m, Geodetic, WGS84};

let lax = Geodetic::new(33.9425, -118.4081, 0.0);
let jfk = Geodetic::new(40.6413, -73.7781, 0.0);

let chord = distance3d_m(lax, jfk);
let spherical = WGS84.great_circle_distance_m(lax, jfk);
let geodesic = WGS84.geodesic_inverse(lax, jfk).unwrap();

assert!(chord < spherical);
assert!((geodesic.distance_m - spherical).abs() / spherical < 0.005);

// The direct problem: go that far on that bearing and you arrive at JFK.
let (end, _final_bearing) = WGS84.geodesic_direct(lax, geodesic.initial_bearing_deg, geodesic.distance_m);
assert!((end.lat_deg - jfk.lat_deg).abs() < 1e-8);

// Spherical helpers for bearing, destination and midpoint are there too.
let bearing = WGS84.initial_bearing_deg(lax, jfk);
let halfway = WGS84.great_circle_midpoint(lax, jfk);
assert!(bearing > 60.0 && bearing < 70.0 && halfway.lon_deg > lax.lon_deg);
```

### Volumes

Everything that can answer "is this point inside?" implements `Volume`. Volumes precompute their
frame at construction, so testing many points against one volume is cheap.

- `AltitudeBand` — a height range with optional floor and ceiling. Compose it with any 2D footprint
  test (from `geo`, geodatafusion, PostGIS, …) via `within_prism` for an extruded airspace or geofence.
- `Sphere` — a range or proximity bubble.
- `Cone` — a directional volume: a camera field of view, an antenna beam, a scan sector.
- `Cylinder` — a circular footprint through a height band.

```rust
use geo3d::{AltitudeBand, Cone, Cylinder, Geodetic, Sphere, Volume};

let site = Geodetic::new(0.0, 0.0, 0.0);

let band = AltitudeBand::new(Some(1_000.0), Some(3_000.0));
let inside_footprint = true; // from your 2D geometry library
assert!(band.within_prism(inside_footprint, 2_000.0));
assert!(!band.contains_height(500.0));

let bubble = Sphere::new(site, 10_000.0);
let camera = Cone::new(site, 90.0, 0.0, 20.0, 5_000.0);        // due east, 20° half-angle, 5 km
let zone = Cylinder::new(site, 8_000.0, AltitudeBand::below(2_000.0));

let p = Geodetic::new(0.0, 0.02, 500.0);                       // ~2.2 km east, 500 m up
assert!(bubble.contains(p) && camera.contains(p) && zone.contains(p));
assert!(!camera.contains(Geodetic::new(0.02, 0.0, 500.0)));    // north: outside the beam
```

### Tracks, courses and conflicts

`Track` is a position and velocity in ECEF. `Course` is what receivers report — course over ground,
ground speed, climb rate — and a `LocalFrame` converts between the two at the object's position.

```rust
use geo3d::{Course, Geodetic, LocalFrame, Track, WGS84};

let a_pos = Geodetic::new(-33.80, 151.30, 2_500.0);
let b_pos = Geodetic::new(-33.90, 151.25, 2_400.0);

// Course over ground → ECEF velocity, each at its own position.
let a = Track::new(WGS84.to_ecef(a_pos), LocalFrame::new(a_pos).velocity_of_course(Course::new(225.0, 120.0, 0.0)));
let b = Track::new(WGS84.to_ecef(b_pos), LocalFrame::new(b_pos).velocity_of_course(Course::new(0.0, 100.0, 0.0)));

let cpa = a.closest_approach(b);            // cpa.time_s (seconds from now), cpa.distance_m
let alert = a.conflict(b, 1_000.0);         // Some(cpa) only if it is in the future and within 1 km
let in_a_minute = WGS84.to_geodetic(a.at(60.0)); // dead reckoning

assert!(cpa.time_s > 0.0 && alert.is_some());
assert!(in_a_minute.lat_deg < a_pos.lat_deg);    // heading south-west

// Back from a velocity to a course, e.g. from two fixes ten seconds apart.
let fixes = Track::from_fixes(WGS84.to_ecef(a_pos), a.at(10.0), 10.0);
let course = LocalFrame::new(in_a_minute).course_of_velocity(fixes.velocity);
assert!((course.course_deg - 225.0).abs() < 0.1);
```

### Line of sight and horizon

Geometric visibility over a smooth ellipsoid — no terrain, no refraction.

```rust
use geo3d::{Geodetic, WGS84};

let observer = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 10_000.0));   // 10 km up
let near = WGS84.to_ecef(Geodetic::new(0.0, 1.0, 10_000.0));       // ~111 km away
let far = WGS84.to_ecef(Geodetic::new(0.0, 10.0, 10_000.0));       // ~1 100 km away

assert!(WGS84.line_of_sight(observer, near));
assert!(!WGS84.line_of_sight(observer, far));                      // the Earth is in the way

let horizon = WGS84.horizon_distance_m(10_000.0);                  // ≈ 357 km
assert!(horizon > 350_000.0 && horizon < 360_000.0);
```

### Inertial frame

The `eci` module rotates between Earth-fixed and Earth-centred inertial coordinates by Greenwich
mean sidereal time (IAU 1982). It is the plain rotation satellite-tracking code uses; precession,
nutation and polar motion are out of scope.

```rust
use geo3d::{eci, Ecef};

let gmst = eci::gmst_from_unix(1_800_000_000.0);   // Unix seconds → sidereal angle in radians
let fixed = Ecef::new(7_000_000.0, 0.0, 0.0);

let inertial = fixed.to_eci(gmst);
let back = inertial.to_ecef(gmst);
assert!((back.x - fixed.x).abs() < 1e-6);
```

### Datum shifts

`Helmert7` applies a 7-parameter transformation between ECEF frames. Build it from parameters as
EPSG publishes them, in either rotation convention.

```rust
use geo3d::{Ecef, Helmert7};

// Position-vector convention: tx ty tz in metres, rx ry rz in arc-seconds, scale in ppm.
let h = Helmert7::position_vector(-446.448, 125.157, -542.060, -0.1502, -0.2470, -0.8421, 20.4894);
// Coordinate-frame convention (flipped rotation signs): Helmert7::coordinate_frame(...)

let p = Ecef::new(3_874_938.0, 116_218.0, 5_047_168.0);
let shifted = h.apply(p);
let back = h.apply_inverse(shifted);           // exact; `h.inverse()` is the first-order approximation
assert!(back.distance_to(p) < 1e-6);
```

### Serde

With the `serde` feature every coordinate type serialises as plain named fields:

```rust,ignore
// cargo add geo3d --features serde
let json = serde_json::to_string(&geo3d::Geodetic::new(34.05, -118.24, 100.0))?;
// {"lat_deg":34.05,"lon_deg":-118.24,"height_m":100.0}
```

### Interoperability

The coordinate types are deliberately small and `Copy`, but you will want to move them into the math
and GIS crates you already use. Optional, **default-off** features add `From`/`Into` conversions, so
the default build stays dependency-free:

- `mint` — `mint::Vector3<f64>` (every vector type) and `mint::Point3<f64>` (the position types).
  `mint` is the common interop layer, so this one bridge also reaches `glam`, `nalgebra` and `cgmath`.
- `glam` — `glam::DVec3`.
- `nalgebra` — `nalgebra::Vector3<f64>` / `Point3<f64>`.
- `geo-types` — `Geodetic` ↔ `geo_types::Point<f64>` / `Coord<f64>`, for piping into the 2D
  [`geo`](https://crates.io/crates/geo) ecosystem (polygon tests, PostGIS, shapefiles). **Longitude
  is `x`, latitude is `y`** (the `geo` convention); the ellipsoidal height is dropped.

```rust
# #[cfg(all(feature = "glam", feature = "geo-types"))] {
use geo3d::{Ecef, Geodetic};

// An ECEF position straight into a glam DVec3 and back.
let e = Ecef::new(4_510_731.0, 4_510_731.0, 0.0);
let v: glam::DVec3 = e.into();
assert_eq!(Ecef::from(v), e);

// A geodetic point into geo-types for 2D GIS (x = lon, y = lat).
let p: geo_types::Point<f64> = Geodetic::new(-33.87, 151.21, 20.0).into();
assert!((p.x() - 151.21).abs() < 1e-9 && (p.y() + 33.87).abs() < 1e-9);
# }
```

```bash
cargo add geo3d --features glam,geo-types
```

## Precision

Each function states its accuracy class, and the tests hold it:

- **Exact closed forms** — geodetic → ECEF, the local frames, look angles, ray–ellipsoid, volumes,
  closest approach, Helmert `apply` / `apply_inverse`. Round trips are held to 1e-11° and 1e-6 m
  over a global grid from 11 km below the surface to geostationary height, and the conversions are
  cross-validated against an independent implementation (`map_3d`, a dev-only dependency).
- **ECEF → geodetic** — Bowring's method iterated to machine precision, so the inverse holds at
  orbital heights as well as on the ground.
- **Geodesics (Vincenty)** — about 0.5 mm. Pinned to GeographicLib reference values: the classic
  Flinders Peak–Buninyong line on two ellipsoids, direct problems from 1 m to 9 000 km, and the exact
  quarter-equator and quarter-meridian lengths. Nearly antipodal pairs return `None` rather than a
  wrong number.
- **Spherical approximations** — `great_circle_*` and `horizon_distance_m` use the mean radius and
  say so; within 0.5 % of the geodesic.
- **Linearised** — `Helmert7::inverse` (use `apply_inverse` for exact), `interpolate` (an ECEF chord),
  `Track` (constant velocity) and `eci` (GMST only) are documented as such on each item.

## Performance

The geometry is closed-form and precomputes per reference, so per-call cost is dominated by the
hardware `sin`/`cos`/`atan2`/`sqrt`. Two things make that measurable and a little faster:

- **Benchmarks.** `cargo bench` runs a [criterion](https://crates.io/crates/criterion) suite over the
  hot paths (conversions, look angles, geodesics, Helmert, ray intersection, vector ops), so a
  regression shows up as a number.
- **FMA, where it pays.** The dot/cross products, the geodetic ↔ ECEF conversions, the Vincenty
  series and the Helmert transform are written as fused multiply-adds — but the fused form is taken
  **only when the target has an FMA instruction**, because `f64::mul_add` lowers to a slow software
  libcall otherwise. On a stock target the code is the plain multiply-add it always was (identical
  speed and results); build with FMA enabled to get the fused path — one instruction, one rounding
  (marginally more accurate), and a few percent faster on the geodesy operations:

  ```bash
  RUSTFLAGS="-C target-cpu=native" cargo build --release   # or -C target-feature=+fma
  ```

## Design

- Zero required dependencies, pure `f64`, every type `Copy`, nothing allocates,
  `#![forbid(unsafe_code)]`. Optional features add `serde` and the interop bridges, nothing else.
- Precompute per reference, not per point: `LocalFrame`, `Sphere`, `Cone` and `Cylinder` do their
  trigonometry once at construction.
- Closed forms wherever one exists; iteration only where it must be, and then to convergence.
- Ellipsoid-generic: everything takes an `Ellipsoid`, WGS84 is the default.
- Degrees and metres at the API; radians never leak out.

Not in scope, by design: planar 2D geometry (use the [`geo`](https://crates.io/crates/geo)
ecosystem), map projections (`proj`, `geodesy`), geodesics to the millimetre over long lines
(`geographiclib`), terrain and meshes (`parry3d`), and precession, nutation or polar motion.

## Develop

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features       # unit, cross-check and doc tests, including this README
cargo run --example tracking    # a ground station tracking an aircraft
cargo publish --dry-run
```

CI runs all of the above plus a build on the declared minimum Rust version (1.85).

## Versioning and license

Plain semver; tags `vX.Y.Z` at the published commit; see [CHANGELOG.md](CHANGELOG.md).
Apache-2.0, see [LICENSE](LICENSE).
