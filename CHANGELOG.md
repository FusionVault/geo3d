# Changelog

All notable changes to geo3d. We use [semantic versioning](https://semver.org), and each
release is tagged `vX.Y.Z`.

## Unreleased

## 0.3.0 — 2026-10-14

An accuracy release. Your code keeps compiling: no public type, method or signature has changed.
Some results are now more accurate, as listed below.

### 🎯 More accurate

- **Distance along the ground now works for every pair of points**
  ([#3](https://github.com/FusionVault/geo3d/issues/3)). `geodesic_inverse` and `geodesic_direct`
  now use Karney's algorithm, the one GeographicLib uses, instead of Vincenty's.
  - It is accurate to about 15 nm, up from about 0.5 mm.
  - `geodesic_inverse` used to return `None` for points on nearly opposite sides of the Earth. Now
    it only returns `None` for NaN or infinite input, or a latitude beyond ±90°.
  - It matches geographiclib-rs to within 4 nm on 25 000 test pairs.
  - The code is ported from GeographicLib and geographiclib-rs, which are both MIT-licensed. See
    NOTICE.
- **Four precision fixes** ([#5](https://github.com/FusionVault/geo3d/issues/5)). Each has a test
  that failed before the fix:
  - **Elevation near straight up.** It used to round to exactly 90°, up to 0.6 m out at
    geostationary range. This affects every look angle.
  - **Great-circle distance.** It was up to 0.9 m short near the far side of the Earth. It is now
    accurate to about a nanometre everywhere.
  - **`Vec3::angle_to`.** Tiny angles used to come out as 0.
  - **Ray–ground intersections.** These lost precision when starting close to the surface. This
    affects `ray_intersect`, `ground_strike` and `line_of_sight`.

### 📖 Clearer docs

- The README and the docs.rs front page were rewritten in plain language. They now have recipes, a
  "find what you need" table, an accuracy table, a glossary, and a lookup table for MATLAB, pymap3d
  and map_3d users.
- There is a new [comparison](comparisons/README.md) against five other Rust crates, with charts.
  It covers accuracy, tricky cases, features and speed, including the places where geo3d is
  slower.
- There is a new [contributing guide](CONTRIBUTING.md).
- A correction: the `eci` frame is TEME, the one SGP4 satellite models use. It is not J2000, which
  differs by about 0.4°. Treating UTC as UT1 can add up to 400 m. The old docs promised "tens of
  metres", which wasn't true.
- docs.rs now shows the optional features. Searching for familiar names like `geodetic2aer`,
  `vdist` or `haversine` now finds the matching geo3d function.

### ⚡ Speed

- `geodesic_direct` is a little faster than before, and about 13 % faster than geographiclib-rs.
- `geodesic_inverse` takes about twice as long as before (about 0.86 µs per call). That extra work
  is what makes it exact everywhere. It is within 4 % of geographiclib-rs.
- Look angles and great-circle distance are 15–20 ns slower per call, because of the precision
  fixes.

### 🧹 Housekeeping

- The code is split into one module per topic.
- The crates.io description is simpler. The `no-std::no-alloc` category was removed, because geo3d
  needs `std` today ([#12](https://github.com/FusionVault/geo3d/issues/12) will change that).
- CI now checks that every code change comes with a CHANGELOG line.

## 0.2.3 — 2026-10-06

Fix: `ground_strike` (`LocalFrame::ground_strike` and the WGS84 shorthand) returns the nearest point of the surface **ahead** of the observer.

- An observer exactly on the surface strikes its own point, whichever way it points. Before, rounding put the near root a hair below zero and the far root was taken: a line 10° down from the ground landed about 2,200 km away, on the far side of the Earth.
- An observer below the surface gives `None`. Before, a line from inside the ellipsoid hit it from within.
- A non-finite height or angle gives `None`.
- Above the surface the result is unchanged. `Ellipsoid::ray_intersect` is unchanged.

## 0.2.2 — 2026-09-14

Additive: optional ecosystem interop and an FMA/benchmark performance pass. No public API removed or changed; the default build stays dependency-free and byte-for-byte identical in its results.

Interop (all default-off features):

- `From`/`Into` between the coordinate types and `mint` (`Vector3`/`Point3`), `glam` (`DVec3`) and `nalgebra` (`Vector3`/`Point3`) — features `mint`, `glam`, `nalgebra`.
- `Geodetic` ↔ `geo_types::Point`/`Coord` (longitude `x`, latitude `y`; height dropped) — feature `geo-types`.
- The interop features track their upstream crates' minimum Rust; the zero-dependency core's MSRV is unchanged (1.85).

Performance:

- Dot/cross products, geodetic ↔ ECEF, the Vincenty series and the Helmert transform are written as fused multiply-adds, gated on `target_feature = "fma"`: stock targets keep the plain multiply-add (identical speed and results — `f64::mul_add` is a slow libcall without FMA hardware), while FMA-enabled builds (`-C target-cpu=native`) get one instruction, one rounding, and a few percent on the geodesy operations.
- A `cargo bench` criterion suite (`benches/geo3d.rs`) over the hot paths.

## 0.2.1 — 2026-09-13

Documentation only: the rewritten README (task-by-task guide, every example doctested) ships in the crate. No code changes.

## 0.2.0 — 2026-09-13

Additive. Closes the gaps against `map_3d` and adds what it does not have.

- NED side everywhere: `LocalFrame::{ned_of, to_ned, ecef_of_ned, geodetic_of_ned}`, `Aer::{from_ned, to_ned}`, `Ned::{vec, to_aer}`; `+ - * neg` operators and array conversions on `Enu`/`Ned`.
- ECEF targets: `LocalFrame::{look_angles_ecef, ecef_of_aer}`.
- `eci` module: Julian date, IAU 1982 GMST, `ecef_to_eci` / `eci_to_ecef` and `Eci` type.
- `Ellipsoid`: `meridional_radius`, `mean_radius`, `geocentric_radius`, `geocentric_latitude`, `surface_normal`, `horizon_distance_m`, `line_of_sight`, `interpolate`, and spherical great-circle `distance / initial_bearing / destination / midpoint`.
- Motion: `Course` (course over ground, ground speed, climb) with `LocalFrame::{course_of_velocity, velocity_of_course}`; `Track` with `from_fixes`, `at`, `closest_approach`, `conflict`.
- Vincenty `geodesic_inverse` / `geodesic_direct` on the ellipsoid (~0.5 mm) with a `Geodesic` result; `asin` inputs clamped; round-trip tests tightened to 1e-11° / 1e-6 m over a global grid and cross-validated against `map_3d` (dev-only).
- `Cylinder` volume; `Aer::{is_above_horizon, normalized}`; `Geodetic::{normalized, antipode}`; `Vec3::{distance_to, angle_to, project_onto}` and `/`.

## 0.1.0 — 2026-09-13

Initial release.

- Typed coordinates: `Geodetic`, `Ecef`, `Enu`, `Ned`, `Aer`, `Vec3`.
- `Ellipsoid` with eight named constants; geodetic ↔ ECEF closed forms; ray–ellipsoid intersection.
- `LocalFrame`: precomputed ENU frame; ENU ↔ ECEF ↔ geodetic, look angles and destination, ground strike.
- Volumes: `AltitudeBand` (+ prism composition), `Sphere`, `Cone`, behind the `Volume` trait.
- `closest_approach` kinematics; `Helmert7` datum shifts in both EPSG rotation conventions.
- Optional `serde` feature.
