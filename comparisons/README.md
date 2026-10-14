# 📊 How geo3d compares

We tested geo3d against five other Rust crates that work with 3D positions on Earth. Every crate got
exactly the same inputs. This page shows what we found, and how you can check it yourself.

The crates are [map_3d](https://crates.io/crates/map_3d), [nav-types](https://crates.io/crates/nav-types),
[sguaba](https://crates.io/crates/sguaba), [jord](https://crates.io/crates/jord) and
[geographiclib-rs](https://crates.io/crates/geographiclib-rs). We also used
[geo](https://crates.io/crates/geo) for its distance functions.

## 🏁 The short version

- 🎯 **Accuracy.** geo3d was exact everywhere we looked, including far out in space, straight up,
  and on opposite sides of the Earth.
- 🛡️ **Tricky cases.** geo3d passed all nine. Every other crate that converts positions failed at
  least one.
- 🧰 **Features.** Only geo3d has ground strikes, 3D zones, collision checks and datum shifts. Others
  are ahead on n-vectors and attitude, which are on our roadmap.
- ⚡ **Speed.** geo3d was fastest, or level with the fastest, in five of the ten tests. It's slower
  in the other five, and we explain why below.

## 🎯 Converting x, y, z back to latitude and longitude

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/roundtrip-dark.svg">
  <img alt="Worst round-trip position error up to geostationary height: geo3d 15 nm, map_3d 81 m" src="charts/roundtrip-light.svg">
</picture>

We converted 66 960 positions per height to x, y, z and back, from 11 km underground to
geostationary orbit. geo3d, nav-types and jord came back exact. map_3d stops after one correction
step, so it drifts as you go higher: about 30 m at GPS satellite height and 81 m at geostationary
height.

<details>
<summary>📋 Show the numbers</summary>

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

</details>

## 📏 Distance along the ground

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/surface-distance-dark.svg">
  <img alt="Worst surface-distance error: geo3d 3.7 nm, great-circle crates 38 km to 51 km" src="charts/surface-distance-light.svg">
</picture>

geo3d works out distances on the Earth's true, slightly squashed shape. It agrees with GeographicLib
to 4 nanometres on 110 000 pairs of points, and never fails. The other position crates only measure
on a perfect sphere, which can be up to 51 km out on a long trip.

Before version 0.3, geo3d used an older method called Vincenty's. That method gave no answer for
7 257 of the 10 000 pairs on opposite sides of the Earth.

<details>
<summary>📋 Show the numbers</summary>

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

| Crate | Sphere radius | Max error | Mean error | Max error (km) |
|---|---|---|---|---|
| geo3d | 6 371 008.771 m (mean radius) | 0.56 % | 0.13 % | 37.8 |
| map_3d | 6 371 000 m | 0.56 % | 0.13 % | 37.8 |
| nav-types | 6 378 137 m (equatorial) | 0.67 % | 0.16 % | 50.7 |
| sguaba | 6 378 137 m (equatorial) | 0.67 % | 0.16 % | 50.7 |
| jord | 6 371 000.8 m | 0.56 % | 0.13 % | 37.8 |
| geo | 6 371 008.8 m | 0.56 % | 0.13 % | 37.8 |

</details>

## 📡 Look angles

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/look-angles-dark.svg">
  <img alt="Worst look-angle error on hard cases: geo3d 30 nm, jord 13 cm, sguaba 13 cm" src="charts/look-angles-light.svg">
</picture>

Every crate agrees on ordinary targets, to well under a micrometre. The differences appear when you
point straight up or exactly along a compass direction. geo3d and map_3d use a formula that stays
precise there. jord and sguaba lose precision, and end up 13 cm off at geostationary range.

<details>
<summary>📋 Show the numbers</summary>

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

</details>

## 🛡️ Tricky cases

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/robustness-dark.svg">
  <img alt="Edge-case scorecard: geo3d passes 9 of 9 cases" src="charts/robustness-light.svg">
</picture>

We tried nine places where maths often breaks: the Earth's centre, the exact poles, the date line,
100 000 km up, opposite sides of the Earth, and inputs that are NaN. geo3d gave a correct answer
every time, and never panicked.

<details>
<summary>📋 Show the numbers</summary>

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

</details>

## 🧰 Features

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/features-dark.svg">
  <img alt="Feature matrix: geo3d is the only crate with ground strike, 3D volumes, closest approach and datum shifts" src="charts/features-light.svg">
</picture>

geo3d covers everything from a GPS fix to a decision: where to point, where a line of sight lands,
whether something is inside a zone, when two objects come closest, and how to shift datums. It does
this with no dependencies and no `unsafe` code.

<details>
<summary>📋 Show the numbers</summary>

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

</details>

## ⚡ Speed

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/speed-dark.svg">
  <img alt="Time per call for 10 operations; geo3d is fastest in 3 and tied in 2" src="charts/speed-light.svg">
</picture>

geo3d is the fastest, or level with the fastest, in five of the ten tests. It's slower in the other
five, and we know why:

- **x, y, z back to latitude and longitude** (two tests): geo3d refines its answer until it's exact.
  A method that is just as exact and faster is next on the list
  ([#28](https://github.com/FusionVault/geo3d/issues/28)).
- **One-off ENU conversions:** a `LocalFrame` costs a little to build, and pays off from the second
  target onwards.
- **Quick great-circle distance:** geo3d uses a formula that stays exact on the far side of the
  Earth. The faster formula most crates use can be almost a metre out there.
- **Distance along the ground:** about 4 % behind geographiclib-rs. geo is faster again because it
  skips the compass bearings.

<details>
<summary>📋 Show the numbers</summary>

Measured on:

- AMD Ryzen 9 5950X (16 cores, 32 threads), Ubuntu 24.04.5 LTS, Linux 7.0.0
- rustc 1.98.0 (88d9e12ae 2026-08-18), default `x86_64-unknown-linux-gnu` target
- criterion 0.5.1: 300 ms warm-up, 1 s measurement, 50 samples

Median nanoseconds per call. In bold: the fastest, and anything within 3 % of it.

| Operation | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo |
|---|---|---|---|---|---|---|---|
| Geodetic → ECEF | **21** | 30 | 36 | 53 | 28 | — | — |
| ECEF → geodetic, −500 m to 12 km | 200 | 111 | 86 | **68** | 104 | — | — |
| ECEF → geodetic, 400 km to 35 786 km | 208 | 111 | **89** | 100 | 107 | — | — |
| ENU, one observer (frame reused) | **28** | 49 | 146 | 59 | 34 | — | — |
| ENU, single call (frame built each time) | 53 | **49** | 185 | 132 | 90 | — | — |
| Look angles, one observer (frame reused) | **67** | 94 | — | 98 | **67** | — | — |
| Look angles, single call | **91** | **93** | — | 170 | 125 | — | — |
| Geodesic inverse | 861 | — | — | — | — | 830 | **768 (distance only)** |
| Geodesic direct | **275** | — | — | — | — | 316 | 286 (position only) |
| Great-circle distance | 60 | 46 | **38** | 53 | 67 | — | **38** |

Inputs: 1 000 points anywhere on Earth from −500 m to 12 km high (and from 400 km to 35 786 km for
the second ECEF → geodetic row); a radar near Sydney with 1 000 aircraft within a degree of it, up
to 12 km high; 1 000 random point pairs; 1 000 lines on random bearings up to 20 000 km long.
geo3d and geographiclib-rs return the distance and both azimuths (inverse) or the position and
final azimuth (direct). geo's `Geodesic` calls geographiclib-rs and returns less. The great-circle
distances use different sphere radii; see [Distance along the ground](#-distance-along-the-ground).

</details>

## 🔬 How we tested

<details>
<summary>The crates, versions and rules</summary>

| Crate | Version | What it is |
|---|---|---|
| geo3d | this repository | 3D geodesy, degrees and metres, no dependencies |
| [map_3d](https://crates.io/crates/map_3d) | 0.1.7 | free functions in radians, in the style of pymap3d and MATLAB |
| [nav-types](https://crates.io/crates/nav-types) | 0.5.2 | typed WGS84 / ECEF / ENU / NED positions on nalgebra |
| [sguaba](https://crates.io/crates/sguaba) | 0.10.7 | typed coordinate frames and rigid-body transforms, with uom units |
| [jord](https://crates.io/crates/jord) | 0.20.0 | n-vector positions, local frames, spherical navigation |
| [geographiclib-rs](https://crates.io/crates/geographiclib-rs) | 0.2.7 | Karney's geodesics; the reference for geodesics |
| [geo](https://crates.io/crates/geo) | 0.33.1 | 2D geometry; only its `Haversine` and `Geodesic` distances are used |

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

</details>

## 🔁 Run it yourself

```bash
cd comparisons
cargo run --release --bin accuracy   # accuracy and tricky cases, as tables
cargo bench                          # speed
cargo test                           # checks every crate's adapter gives the same simple answer
```

To redraw the charts after a change:

```bash
cargo run --release --bin accuracy > results/accuracy.md
cargo bench
python3 charts.py
```

The inputs come from a fixed random seed, so every run tests the same points. The charts are drawn
by [`charts.py`](charts.py) from [`results/accuracy.md`](results/accuracy.md) and
[`results/speed.csv`](results/speed.csv). This folder is its own Rust project and isn't published,
so none of these crates become dependencies of geo3d.
