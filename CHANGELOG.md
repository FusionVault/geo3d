# Changelog

All notable changes to geo3d. We use [semantic versioning](https://semver.org), and each
release is tagged `vX.Y.Z`.

## Unreleased

- 🔌 **Runs without the standard library** ([#12](https://github.com/FusionVault/geo3d/issues/12)).
  On microcontrollers and flight controllers, `default-features = false, features = ["libm"]` takes
  the maths from `libm`. Everything else is the same, every optional feature works too, and CI
  builds it for an Arm Cortex-M4F. Nothing changes for anyone else: the new `std` feature is on by
  default, and the default build still has no dependencies.
  - ⚠️ If you already set `default-features = false`, add `features = ["std"]` (or `libm`), or
    the build will stop with a message saying so.
  - Keep `std` where you have it: its sine, cosine and atan2 round the last digit more carefully
    than libm's (worst 0.5 against 0.7–1.2 units in the last place). libm is 10–20 % faster.
- 🧭 **Route tools** ([#7](https://github.com/FusionVault/geo3d/issues/7)), built on n-vectors so
  the poles and the antimeridian need no special care:
  - `cross_track_distance_m` and `along_track_distance_m`: how far off a route, and how far along;
  - `great_circle_intersection`: where two routes cross;
  - `mean_position`: the centre of a group of positions;
  - `great_circle_interpolate`: a point part-way along a route, at constant speed.

  They stay exact on short paths. The textbook n-vector formula can be kilometres out on a path a
  few micrometres long; geo3d builds the path from the latitude and longitude differences instead.
- 🐛 **`great_circle_destination` near the poles.** It computed the latitude with `asin`, which
  loses precision near ±90°: up to 11 cm out arriving close to a pole. It now uses n-vectors and
  `atan2`, and is exact everywhere.
- 📚 **The NavLab examples** (`geo3d::navlab`): nine of FFI's ten classic n-vector problems, solved
  with geo3d and checked against the published answers. (Example 2 needs attitude, issue #13.)
- 🎯 **Uncertainty between frames** ([#6](https://github.com/FusionVault/geo3d/issues/6)). `Cov3`
  is a 3×3 covariance. It moves between frames with the position:
  - ENU ↔ ECEF (`LocalFrame::enu_cov_to_ecef`, `ecef_cov_to_enu`) and ENU ↔ NED, exactly;
  - AER ↔ ENU (`Aer::cov_to_enu`, `Enu::cov_to_aer`) and geodetic ↔ ECEF
    (`Ellipsoid::geodetic_cov_to_ecef`, `ecef_cov_to_geodetic`), to first order.

  The Jacobians are checked against finite differences, and the results against a Monte Carlo
  run. `Cov3` converts to `nalgebra::Matrix3`, `glam::DMat3` and the `mint` matrices.
- 🗺️ **Polygon area** ([#4](https://github.com/FusionVault/geo3d/issues/4)).
  `Ellipsoid::polygon_area` gives the area and perimeter of a polygon whose edges are geodesics. It
  is a port of GeographicLib's method, and matches it to within 0.1 m² on any polygon, including
  ones around a pole or across the antimeridian. Counter-clockwise polygons have a positive area.
- 🛰️ **Satellite time** ([#8](https://github.com/FusionVault/geo3d/issues/8)):
  - `eci::gmst_from_unix_with_dut1` corrects UTC to UT1, which removes up to 400 m of error.
  - `eci::era_rad` and `eci::era_from_unix_with_dut1` give the Earth Rotation Angle (IAU 2000).
  - Sidereal time from Unix time now resolves 0.2 µs instead of 40 µs.
  - Both formulas match SOFA's reference values to 1e-12 rad.
  - A new example, `examples/satellite.rs`, tracks the ISS with the `sgp4` crate.
- 📡 **Refraction for radar and radio** ([#9](https://github.com/FusionVault/geo3d/issues/9)).
  `horizon_distance_with_refraction_m` and `line_of_sight_with_refraction` take an
  effective-Earth-radius factor `k`. `STANDARD_REFRACTION_K` is 4/3, the ITU-R standard for radio,
  which reaches about 15 % further than a straight line.
- 🐛 **`line_of_sight` between points exactly on the surface.** It used to say two points on the
  ground could see each other at any distance, even when the line between them went hundreds of
  metres underground. Now it checks the deepest point of the line, so they only see each other if
  the line stays above the surface.
- 🛡️ **An edge-case test suite** ([#10](https://github.com/FusionVault/geo3d/issues/10)). It runs
  every public conversion at the exact poles, the antimeridian, the Earth's centre, 100 000 km up,
  deep underground, exact antipodes, and NaN or infinite input. It found and fixed two bugs:
  - Bearings and azimuths could come out as exactly 360° (from a tiny negative angle). They are
    now always in `[0, 360)`. This affects look angles, `Course`, `initial_bearing_deg` and
    geodesic bearings.
  - `ground_strike` returned `Some(NaN position)` for a NaN latitude or longitude at height 0.
    It now returns `None`.
- 🧹 CI now runs [cargo-semver-checks](https://github.com/obi1kenobi/cargo-semver-checks) on every
  pull request, so accidental breaking changes are caught before release
  ([#11](https://github.com/FusionVault/geo3d/issues/11)).
- 🎯 **Faster, more accurate x, y, z → latitude/longitude/height**
  ([#28](https://github.com/FusionVault/geo3d/issues/28)). `to_geodetic` now uses an exact formula
  (Vermeille's, refined by Karney for GeographicLib) and one correction step. It no longer
  iterates.
  - It is about 40 % faster.
  - It is at least as accurate at every height: within 3 nm near the Earth and 20 nm at satellite
    heights. This was measured at 50-digit precision on 5 000 points.
  - Deep inside the Earth, the old method could be up to 53 m out. The new one stays within 3 nm.
  - A point exactly on the polar axis now always reports longitude 0.

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
