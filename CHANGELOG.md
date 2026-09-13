# Changelog

All notable changes to `geo3d`. Plain semver; tags `vX.Y.Z` at the published commit.

## 0.1.0 — 2026-09-13

Initial release.

- Typed coordinates: `Geodetic`, `Ecef`, `Enu`, `Ned`, `Aer`, `Vec3`.
- `Ellipsoid` with eight named constants; geodetic ↔ ECEF closed forms; ray–ellipsoid intersection.
- `LocalFrame`: precomputed ENU frame; ENU ↔ ECEF ↔ geodetic, look angles and destination, ground strike.
- Volumes: `AltitudeBand` (+ prism composition), `Sphere`, `Cone`, behind the `Volume` trait.
- `closest_approach` kinematics; `Helmert7` datum shifts in both EPSG rotation conventions.
- Optional `serde` feature.
