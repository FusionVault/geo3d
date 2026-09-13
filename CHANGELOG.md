# Changelog

All notable changes to `geo3d`. Plain semver; tags `vX.Y.Z` at the published commit.

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
