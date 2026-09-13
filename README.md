# geo3d

3D geodesy in pure `f64`, zero dependencies.

Real 3D geometry cannot be done on geodetic `(lat, lon, height)` directly — it is not a Euclidean
space. This crate provides the metric frames and the closed-form primitives on top of them, in a small,
typed, allocation-free API:

| | |
|---|---|
| `Ellipsoid` (`WGS84`, `GRS80`, …) | geodetic ↔ ECEF (closed form both ways), prime-vertical radius, ray–ellipsoid intersection |
| `Geodetic`, `Ecef`, `Enu`, `Ned`, `Aer` | typed coordinates — lat/lon or east/north cannot be swapped by accident |
| `LocalFrame` | a topocentric East-North-Up frame that precomputes its rotation once; ENU ↔ ECEF ↔ geodetic, look angles (az/el/range) and their inverse, ground strike |
| `AltitudeBand`, `Sphere`, `Cone`, `Cylinder` (`Volume`) | closed-form volumes; compose with any 2D footprint test for extruded prisms |
| `Track`, `Course`, `closest_approach` | dead reckoning, course-over-ground ↔ ECEF velocity, first-order conflict detection |
| `Ellipsoid` extras | radii of curvature, geocentric latitude, horizon distance, geometric line of sight, ECEF interpolation, spherical great-circle distance / bearing / destination / midpoint |
| `eci` | Julian date, GMST (IAU 1982), inertial ↔ Earth-fixed |
| `Helmert7` | 7-parameter datum shifts, both EPSG rotation conventions |

Angles in the public API are degrees, lengths metres. ECEF is the frame CesiumJS renders in, so a
`(lat, lon, height)` track maps straight to what a globe draws.

```rust
use geo3d::{distance3d_m, Geodetic, LocalFrame, Sphere, Volume};

let site = Geodetic::new(34.05, -118.24, 100.0);
let aircraft = Geodetic::new(34.30, -118.10, 3_000.0);

let slant = distance3d_m(site, aircraft);          // through-space distance, metres
let frame = LocalFrame::new(site);
let aer = frame.look_angles(aircraft);              // azimuth / elevation / range
let back = frame.destination(aer);                  // ...and back to the aircraft
let ground = frame.ground_strike(aer.azimuth_deg, -3.0); // where a 3° down-look meets the surface
assert!(Sphere::new(site, 40_000.0).contains(aircraft));
```

`cargo run --example tracking` walks through a ground station tracking an aircraft.

## Design

- **Zero dependencies, pure `f64`.** Every type is `Copy`; nothing allocates; `#![forbid(unsafe_code)]`.
- **Precompute per reference, not per point.** `LocalFrame`, `Sphere` and `Cone` do their
  trigonometry once at construction, so "many targets against one observer" is a subtraction and a
  few dot products per point.
- **Closed forms.** Geodetic ↔ ECEF uses the exact forward form and Bowring's inverse
  (sub-millimetre at terrestrial and orbital heights, no iteration).
- **Ellipsoid-generic.** Everything takes an `Ellipsoid`; WGS84 is the default and has shorthand
  free functions.
- **Optional `serde`** feature derives `Serialize`/`Deserialize` on the coordinate types.

Out of scope, by design: planar 2D geometry (use the [`geo`](https://crates.io/crates/geo)
ecosystem), map projections (`proj`, `geodesy`), geodesics to millimetre accuracy over long lines
(`geographiclib`), terrain and meshes (`parry3d`), and precession/nutation/polar motion (the `eci`
module is the plain GMST rotation).

## Develop

```bash
cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all-features
cargo publish --dry-run
```

## Versioning

Plain semver, tags `vX.Y.Z` at the published commit. See [CHANGELOG.md](CHANGELOG.md).
