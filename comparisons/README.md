# geo3d compared

How geo3d compares with other Rust crates on accuracy, robustness, features and speed, with the code
to reproduce every number. This crate is not published and is a workspace of its own, so none of
these crates ever reach geo3d's dependencies.

The numbers are for geo3d with Karney geodesics and the precision fixes of pull requests #35 and
#36. The suite calls geo3d only through its public API.

## Summary

**Accuracy.**

- Round trips (geodetic → ECEF → geodetic) are exact to machine precision at every height from
  11 km underground to geostationary orbit: 2e-14° and 0.02 µm. nav-types and jord match that;
  sguaba is within 3 µm; map_3d is off by up to 1.1e-4° in latitude at geostationary height.
- Geodesics agree with geographiclib-rs (GeographicLib's algorithm) to 3.7 nm on 110 000 pairs,
  10 000 of them nearly antipodal, with no failures.
- ENU positions agree with every other crate to 0.04 µm, and look angles with map_3d to 0.08 µm.
  Straight up and on the horizon, geo3d and map_3d are the only crates that stay exact: jord and
  sguaba are off by up to 5e-7° near straight up (0.13 m at geostationary range), and sguaba snaps
  targets 1e-7° off due east or west onto it.

**Robustness.** geo3d gives a correct, finite answer in every edge case tested: the Earth's centre,
the exact poles, 100 000 km altitude and exact antipodes. NaN input gives NaN output or `None`, and
nothing panics. nav-types returns NaN at the exact poles, sguaba panics on a NaN ECEF input, map_3d
and jord return NaN at the Earth's centre, and map_3d is wrong at 100 000 km.

**Features.** geo3d is the only one of these crates with a ray–ellipsoid ground strike, 3D volumes,
3D closest approach and Helmert datum shifts, with no dependencies and `#![forbid(unsafe_code)]`.
It has no n-vectors yet (nav-types and jord do) and no attitude or body frames (sguaba and jord do);
both are on the roadmap.

**Speed.** geo3d is the fastest at geodetic → ECEF, at ENU from a reused `LocalFrame`, at a
single look-angle call and at the geodesic direct problem, and ties jord for look angles from a
reused frame. Its geodesic inverse is about 5 % slower than geographiclib-rs's. Two places are
slower by choice or for now:

- **ECEF → geodetic** takes about 207 ns against 71–117 ns for the others. geo3d iterates Bowring's
  method to full convergence. A closed form that is both faster and more accurate (nanometres) is
  planned ([#28](https://github.com/FusionVault/geo3d/issues/28)).
- **Great-circle distance** takes 60 ns against 38–68 ns. geo3d uses the angle between unit vectors,
  which is accurate to about a nanometre at every distance; the faster haversine is up to 0.9 m off
  near the antipode.

## What is compared

| Crate | Version | What it is |
|---|---|---|
| geo3d | this repository | 3D geodesy, degrees and metres, no dependencies |
| [map_3d](https://crates.io/crates/map_3d) | 0.1.7 | free functions in radians, in the style of pymap3d and MATLAB |
| [nav-types](https://crates.io/crates/nav-types) | 0.5.2 | typed WGS84 / ECEF / ENU / NED positions on nalgebra |
| [sguaba](https://crates.io/crates/sguaba) | 0.10.7 | typed coordinate frames and rigid-body transforms, with uom units |
| [jord](https://crates.io/crates/jord) | 0.20.0 | n-vector positions, local frames, spherical navigation |
| [geographiclib-rs](https://crates.io/crates/geographiclib-rs) | 0.2.7 | Karney's geodesics; the reference for geodesics |
| [geo](https://crates.io/crates/geo) | 0.33.1 | 2D geometry; only its `Haversine` and `Geodesic` distances are used |

Operations: geodetic → ECEF, ECEF → geodetic, geodetic → ENU, geodetic → look angles (azimuth,
elevation, range), geodesic inverse and direct, and great-circle distance. "—" means the crate
does not offer that operation.

All crates get the same inputs, in degrees and metres. Each one is called through a small adapter
in [`src/lib.rs`](src/lib.rs) that does what its users would write: degrees to radians for
map_3d, building the typed values for nav-types, sguaba and jord, and reading the results back out.
That conversion work is part of every timing. Inputs and results go through `std::hint::black_box`.

- **One observer, many targets**: each crate keeps what its API lets you keep, built once outside
  the timed loop. geo3d keeps a `LocalFrame`, jord an `EnuFrame`, sguaba a `RigidBodyTransform`.
  nav-types keeps the observer's ECEF position, but every subtraction converts it back to latitude
  and longitude to build the rotation. map_3d has no frame object, so it keeps the observer in
  radians.
- **Single call**: the frame is built on every call.
- **Checked constructors**: nav-types is called through `try_from_degrees_and_meters` and sguaba
  through its builder, which both return `None` for a bad latitude. (nav-types'
  `from_degrees_and_meters` panics instead.)
- **Release builds** for both speed and accuracy. In debug builds sguaba also checks every
  conversion with a round trip and panics outside its supported altitude range.
- **No `target-cpu=native`**. geo3d's fused multiply-add paths are not used.

## Reproduce

```bash
cd comparisons
cargo bench                          # speed; each reported time is for 1 000 calls
cargo run --release --bin accuracy   # accuracy and robustness, as Markdown tables
cargo test                           # every adapter gives the same answer on a simple case
```

The inputs come from a fixed-seed xorshift generator, so every run uses the same points.

## Accuracy

From `cargo run --release --bin accuracy`.

### Round trip: geodetic → ECEF → geodetic

66 960 points per height: a 1° grid plus latitudes ±89.9°, ±89.99° and ±89.999°. Each crate uses
its own forward and inverse conversion. No crate failed on any of these points.

Max latitude error (degrees):

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| −11 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 3.6e-14 | 2.1e-14 |
| 0 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 10 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 400 km | 2.1e-14 | 2.5e-9 | 1.4e-14 | 4.3e-14 | 2.1e-14 |
| 20 200 km | 1.4e-14 | 7.0e-5 | 2.1e-14 | 5.0e-14 | 2.1e-14 |
| 35 786 km | 1.8e-14 | 1.1e-4 | 1.4e-14 | 4.3e-14 | 2.1e-14 |

Max longitude error (degrees): 1.4e-14 to 5.7e-14 for every crate at every height.

Max height error (metres):

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| −11 km | 1.9e-9 | 2.0e-9 | 3.0e-9 | 1.9e-7 | 2.7e-9 |
| 0 km | 1.9e-9 | 2.1e-9 | 1.4e-9 | 1.2e-9 | 2.1e-9 |
| 10 km | 1.9e-9 | 2.2e-9 | 2.5e-9 | 2.6e-9 | 3.6e-9 |
| 400 km | 2.8e-9 | 2.2e-9 | 2.7e-9 | 3.2e-7 | 3.2e-9 |
| 20 200 km | 1.1e-8 | 6.3e-6 | 1.1e-8 | 1.2e-6 | 1.1e-8 |
| 35 786 km | 1.5e-8 | 1.4e-5 | 1.5e-8 | 2.2e-6 | 1.5e-8 |

map_3d's inverse takes one fixed correction step. It is exact near the surface, but at GPS-orbit
height it is off by 7e-5° in latitude (about 30 m along the orbit), and by 1.1e-4° at
geostationary height.

### Geodesics, against geographiclib-rs

100 000 random point pairs (uniform over the sphere) and 10 000 nearly antipodal pairs (the second
point within 0.5° of the first one's antipode). A failure is `None`, NaN or a panic. geo's
`Geodesic` calls geographiclib-rs, so it matches it exactly.

Inverse: distance and azimuths between two points.

| Crate | Pairs | Count | Failures | Max distance error (m) | Mean distance error (m) | Max azimuth error (°) |
|---|---|---|---|---|---|---|
| geo3d | random | 100 000 | 0 | 3.7e-9 | 4.7e-13 | 5.7e-14 |
| geo (`Geodesic`) | random | 100 000 | 0 | 0 | 0 | — |
| geo3d | nearly antipodal | 10 000 | 0 | 3.7e-9 | 3.7e-13 | 4.5e-13 |
| geo (`Geodesic`) | nearly antipodal | 10 000 | 0 | 0 | 0 | — |

Direct: travel geographiclib-rs's distance on its initial azimuth, and compare the arrival with its
own direct solution.

| Crate | Lines | Count | Failures | Max position error (m) | Mean position error (m) | Max final azimuth error (°) |
|---|---|---|---|---|---|---|
| geo3d | random | 100 000 | 0 | 0 | 0 | 0 |
| geo (`Geodesic`) | random | 100 000 | 0 | 0 | 0 | — |
| geo3d | nearly antipodal | 10 000 | 0 | 0 | 0 | 0 |
| geo (`Geodesic`) | nearly antipodal | 10 000 | 0 | 0 | 0 | — |

Before pull request #35, geo3d used Vincenty's method: within 0.1 mm of geographiclib-rs when it
succeeded, but it returned `None` for 7 257 of the 10 000 nearly antipodal pairs and for 1 of the
100 000 random ones.

Great-circle (spherical) distance against the ellipsoidal geodesic, all 110 000 pairs:

| Crate | Sphere radius | Max error | Mean error |
|---|---|---|---|
| geo3d | 6 371 008.771 m (mean radius) | 0.56 % | 0.13 % |
| map_3d | 6 371 000 m | 0.56 % | 0.13 % |
| nav-types | 6 378 137 m (equatorial) | 0.67 % | 0.16 % |
| sguaba | 6 378 137 m (equatorial) | 0.67 % | 0.16 % |
| jord | 6 371 000.8 m | 0.56 % | 0.13 % |
| geo | 6 371 008.8 m | 0.56 % | 0.13 % |

### Local frames and look angles

Every crate is compared with geo3d and with map_3d. For targets built from a known pointing, that
pointing is the reference too, and the last column gives the angle between the two pointings.
Disagreements are in micrometres. For look angles, each answer (azimuth, elevation, range) is
turned back into a point, and the table gives the distance between the two points, so a pointing
error counts by how far it moves the target.

15 000 observer–target pairs: nearby aircraft, and targets anywhere up to geostationary height.

| Crate | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs geo3d (µm) | Look angles vs map_3d (µm) |
|---|---|---|---|---|
| geo3d | 0 | 0.017 | 0 | 0.077 |
| map_3d | 0.017 | 0 | 0.077 | 0 |
| nav-types | 0.018 | 0.017 | — | — |
| sguaba | 0.037 | 0.038 | 8.2 | 8.2 |
| jord | 0.022 | 0.020 | 0.081 | 0.060 |

2 400 pairs on the horizon: due north, east, south and west, and 1e-7° and 1e-5° off them, at 10 km
and 1 000 km.

| Crate | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing error (°) |
|---|---|---|---|---|
| geo3d | 0 | 0.002 | 0.006 | 3.3e-11 |
| map_3d | 0.002 | 0 | 0.006 | 3.7e-11 |
| nav-types | 0.001 | 0.002 | — | — |
| sguaba | 0.001 | 0.002 | 1 745 | 1.0e-7 |
| jord | 0.004 | 0.004 | 0.006 | 3.4e-11 |

2 000 pairs straight up and 1e-7° to 1e-3° off it, from 1 km to geostationary range.

| Crate | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing error (°) |
|---|---|---|---|---|
| geo3d | 0 | 0.015 | 0.030 | 2.1e-10 |
| map_3d | 0.015 | 0 | 0.030 | 2.1e-10 |
| nav-types | 0.015 | 0.017 | — | — |
| sguaba | 0.027 | 0.030 | 129 550 | 4.8e-7 |
| jord | 0.023 | 0.023 | 129 550 | 4.8e-7 |

ENU agrees everywhere to a few hundredths of a micrometre. The look angles differ in how the angles
are computed:

- geo3d and map_3d use `atan2` for both angles and stay exact everywhere. (Before pull request #36,
  geo3d computed the elevation as `asin(up / range)` and was off by up to 1e-6° near straight up.)
- jord computes the elevation with `asin`, which loses precision within about 1e-6° of straight up:
  it is off by up to 5e-7° there (0.13 m at geostationary range).
- sguaba uses `acos` for both angles. A target 1e-7° off due east or west comes out exactly due
  east or west (1.7 mm at 1 000 km), and near straight up it is off by up to 5e-7°. Its 8 µm on the
  general targets also comes from its angle formulas: its ENU agrees to 0.04 µm.

## Robustness

From the same run. "ok" means within 1 mm and 1e-8°. For NaN input, NaN out or a rejection
(`None`) are both fine; a panic is not.

| Case | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo |
|---|---|---|---|---|---|---|---|
| ECEF → geodetic at the Earth's centre (0, 0, 0); want ±90°, h = −b | ok | NaN | NaN latitude, h −6 292 742 m | NaN longitude (lat 90°, h −b) | NaN | — | — |
| ECEF → geodetic exactly above the North Pole (0, 0, b + 100 m) | ok | ok | NaN latitude and height | ok | ok | — | — |
| Round trip exactly at the South Pole (−90°, 0°, 100 m) | ok | ok | NaN latitude and height | ok | ok | — | — |
| ENU across the antimeridian, 179.9999° → −179.9999° | ok | ok | ok | ok | ok | — | — |
| Round trip at 100 000 km altitude (45°, 45°) | ok | wrong: latitude 45.000126° | ok | ok | ok | — | — |
| Geodetic → ECEF with a NaN latitude | NaN | NaN | rejected (`None`) | rejected (`None`) | NaN | — | — |
| ECEF → geodetic with a NaN x | NaN | NaN | NaN | **panic** | NaN | — | — |
| Geodesic between exact antipodes (0°, 0°) → (0°, 180°) | ok (20 003 931.459 m) | — | — | — | — | ok | ok |
| Geodesic with a NaN latitude | rejected (`None`) | — | — | — | — | NaN | NaN |

sguaba's panic is an `expect` on its own latitude check, in release builds too.

## Features and footprint

Checked in each crate's source at the versions above. Dependencies are counted with default
features, from `cargo tree -e normal`.

| | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo |
|---|---|---|---|---|---|---|---|
| Direct dependencies | 0 | 0 | 1 | 5 | 0 | 2 | 13 |
| All crates pulled in | 0 | 0 | 19 | 21 | 0 | 5 | 41 |
| `unsafe` | `#![forbid(unsafe_code)]` | none, not forbidden | none, not forbidden | marks frame assumptions: users write `unsafe` to build ECEF ↔ ENU/NED transforms; no memory-unsafe code | `#![deny(unsafe_code)]` | none, not forbidden | 2 blocks |
| MSRV | 1.85 | not declared | not declared | 1.89 | 1.85 | 1.70 | 1.88 |
| Ellipsoids | 8 named, or any | 10 named | WGS84 | WGS84 | 4 named, or any; spheres | any | WGS84 (`Geodesic`); spheres |
| Look angles | yes | yes | no (ENU/NED only) | yes | yes | — | — |
| Ray–ellipsoid ground strike | yes | no | no | no | no | — | — |
| 3D volumes | sphere, cone, cylinder, altitude band | no | no | no | no (caps and polygons on a sphere) | — | no (2D shapes) |
| Closest approach | 3D, constant velocity | no | no | no | on a sphere (surface tracks) | — | no |
| Datum shift (Helmert) | yes | no | no | no | no | — | no |
| Inertial frame | TEME (GMST rotation) | GMST rotation (`ecef2eci`) | no | no | no | — | no |
| Geodesics | Karney | haversine only | haversine only | haversine only | great circles on a sphere | Karney | Karney (geographiclib-rs) |
| n-vector | no | no | yes | no | yes | no | no |
| Attitude (yaw, pitch, roll) | no | no | no | yes | yes | no | no |

## Speed

Measured on:

- AMD Ryzen 9 5950X (16 cores, 32 threads), Ubuntu 24.04.5 LTS, Linux 7.0.0
- rustc 1.98.0 (88d9e12ae 2026-08-18), default `x86_64-unknown-linux-gnu` target
- criterion 0.5.1: 300 ms warm-up, 1 s measurement, 50 samples

Median time per call in nanoseconds (criterion's median for 1 000 calls, divided by 1 000). Runs
agree within about 7 %. Fastest in bold.

| Operation | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo |
|---|---|---|---|---|---|---|---|
| Geodetic → ECEF | **22** | 29 | 37 | 53 | 29 | — | — |
| ECEF → geodetic, −500 m to 12 km | 206 | 117 | 87 | **71** | 105 | — | — |
| ECEF → geodetic, 400 km to 35 786 km | 208 | 111 | **89** | 101 | 107 | — | — |
| ENU, one observer (frame reused) | **28** | 49 | 147 | 59 | 34 | — | — |
| ENU, single call (frame built each time) | 53 | **49** | 187 | 132 | 90 | — | — |
| Look angles, one observer (frame reused) | **67** | 94 | — | 98 | **67** | — | — |
| Look angles, single call | **92** | 93 | — | 171 | 127 | — | — |
| Geodesic inverse | 873 | — | — | — | — | 830 | **767** (distance only) |
| Geodesic direct | **275** | — | — | — | — | 316 | 287 (position only) |
| Great-circle distance | 60 | 47 | 38 | 53 | 68 | — | **38** |

Inputs: 1 000 points anywhere on Earth from −500 m to 12 km high (and from 400 km to 35 786 km for
the second ECEF → geodetic row); a radar near Sydney with 1 000 aircraft within a degree of it, up
to 12 km high; 1 000 random point pairs; 1 000 lines on random bearings up to 20 000 km long.
geo3d and geographiclib-rs return the distance and both azimuths (inverse) or the position and
final azimuth (direct). geo's `Geodesic` calls geographiclib-rs and returns less. The great-circle
distances use different sphere radii; see [Accuracy](#accuracy).
