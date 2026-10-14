# Changelog

All notable changes to `geo3d`. Plain semver; tags `vX.Y.Z` at the published commit.

## Unreleased

- **Geodesics now use Karney's algorithm**, as in GeographicLib, instead of Vincenty's ([#3](https://github.com/FusionVault/geo3d/issues/3)).
  `geodesic_inverse` now works for every pair of points. It used to return `None` for points on
  nearly opposite sides of the Earth. Both `geodesic_inverse` and `geodesic_direct` are accurate to
  about 15 nm, up from 0.5 mm. The signatures are unchanged. `geodesic_inverse` returns `None` only
  for non-finite input or a latitude beyond ±90°. Ported from GeographicLib via geographiclib-rs
  (both MIT; see NOTICE), and cross-checked against geographiclib-rs on 25 000 pairs on four
  ellipsoids, where they agree to within 4 nm.
- Speed: `geodesic_direct` is slightly faster than before (about 0.22 µs per call) and about 10 %
  faster than geographiclib-rs. `geodesic_inverse` costs about 0.74 µs on a typical pair, roughly
  twice Vincenty's (it does the extra work that makes it exact everywhere) and within a few
  percent of geographiclib-rs. WGS84's series coefficients are built at compile time.

## 0.2.4 — unreleased

Documentation and internal structure. No public API change: every type, method and trait impl is
unchanged, and results are identical.

- Rewritten README and crate docs: a short pitch, a "find what you need" table, task-by-task
  recipes, an accuracy table, conventions (units, ranges, height reference), a glossary, and a
  lookup table for people coming from MATLAB, pymap3d or map_3d.
- Corrected the `eci` docs: `Eci` is the TEME frame that SGP4 / TLE propagators output, not
  J2000 / GCRF (about 0.4° apart), and `gmst_from_unix` treats UTC as UT1 (up to about 400 m at
  the equator). The old "tens of metres" claim did not hold against J2000 / GCRF data.
- docs.rs now builds with all features, so the `serde` and interop impls show up; the `tracking`
  example is scraped into the docs; search aliases for the MATLAB / pymap3d names
  (`geodetic2ecef`, `ecef2enu`, `geodetic2aer`, `vdist`, `vreckon`, `haversine`, …).
- crates.io metadata: a plainer description; `gps` and `coordinates` keywords in place of `enu`
  and `geospatial`; dropped the `no-std::no-alloc` category (the crate needs `std` for its float
  maths); added `aerospace` and `science::robotics`.
- Source split into one module per topic (vectors, positions, local frames, ellipsoid, geodesics,
  great circles, visibility, motion, volumes, datum, ECI); `Course` moved next to `Track`.

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
