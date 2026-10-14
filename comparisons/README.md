# 📊 How geo3d compares

We tested geo3d against six other Rust crates that work with positions on Earth. Every crate got
exactly the same inputs. This page shows what we found, and how you can check it yourself.

The crates are [map_3d](https://crates.io/crates/map_3d), [nav-types](https://crates.io/crates/nav-types),
[sguaba](https://crates.io/crates/sguaba), [jord](https://crates.io/crates/jord) and
[geographiclib-rs](https://crates.io/crates/geographiclib-rs), and [utm](https://crates.io/crates/utm) for
grid coordinates. We also used
[geo](https://crates.io/crates/geo) for its distance, area and route functions.

## 🏁 The short version

- 🎯 **Accuracy.** geo3d was exact everywhere we looked, including far out in space, straight up,
  and on opposite sides of the Earth.
- 🗺️ **Area.** geo3d's polygon areas match GeographicLib everywhere, even around the poles. Methods
  that treat the Earth as a ball are 0.7 % to 1.4 % out.
- 🛡️ **Tricky cases.** geo3d passed all nine. Every other crate that converts positions failed at
  least one.
- ✈️ **Attitude.** geo3d keeps track of yaw, pitch and roll even pointing straight up. sguaba and
  jord can lose it there.
- 🧰 **Features.** Only geo3d has ground strikes, 3D zones, collision checks, uncertainty and datum
  shifts.
- ⚡ **Speed.** geo3d was fastest, or level with the fastest, in six of the thirteen tests. It's slower
  in the other five, and we explain why below.

## 🎯 Converting x, y, z back to latitude and longitude

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/roundtrip-dark.svg">
  <img alt="Worst round-trip position error up to geostationary height: geo3d 21 nm, map_3d 81 m" src="charts/roundtrip-light.svg">
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
| −11 km | 1.4e-14 | 2.1e-14 | 2.1e-14 | 3.6e-14 | 2.1e-14 |
| 0 km | 1.4e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 10 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 400 km | 1.4e-14 | 2.5e-9 | 1.4e-14 | 4.3e-14 | 2.1e-14 |
| 20 200 km | 1.4e-14 | 7.0e-5 | 2.1e-14 | 5.0e-14 | 2.1e-14 |
| 35 786 km | 2.1e-14 | 1.1e-4 | 1.4e-14 | 4.3e-14 | 2.1e-14 |

Max longitude error (degrees): 1.4e-14 to 5.7e-14 for every crate at every height.

Max height error (metres):

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| −11 km | 2.1e-9 | 2.0e-9 | 3.0e-9 | 1.9e-7 | 2.7e-9 |
| 0 km | 2.1e-9 | 2.1e-9 | 1.4e-9 | 1.2e-9 | 2.1e-9 |
| 10 km | 2.1e-9 | 2.2e-9 | 2.5e-9 | 2.6e-9 | 3.6e-9 |
| 400 km | 2.0e-9 | 2.2e-9 | 2.7e-9 | 3.2e-7 | 3.2e-9 |
| 20 200 km | 7.5e-9 | 6.3e-6 | 1.1e-8 | 1.2e-6 | 1.1e-8 |
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

## 🗺️ Polygon area

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/polygon-area-dark.svg">
  <img alt="Polygon area scorecard: geo3d matches GeographicLib on all four kinds of polygon; spherical methods are 0.7 % to 1.4 % off, and geo's fails across the antimeridian and around the poles" src="charts/polygon-area-light.svg">
</picture>

We measured 4 000 polygons: small fields, regions up to 100 km across, shapes that cross the
antimeridian, and rings around both poles. geo3d matched GeographicLib on all of them, to 0.12 m²
at worst. That's on rings around a pole, where the area is what's left when edge terms of up to
20 million km² cancel, so a rounding step in one sine shows. geo's geodesic area calls
geographiclib-rs, so it matches too.

The spherical methods in geo and jord treat the Earth as a ball, so they are about 0.7 % out, and
jord up to 1.4 % on small fields. geo's spherical method also treats longitude as a plain number
line. Across ±180° that makes some answers thousands of times too big, and around a pole it
returns zero.

<details>
<summary>📋 Show the numbers</summary>

1 000 polygons of each kind, with 3 to 32 vertices, going counter-clockwise. Longitudes are in
[−180°, 180°), as GeoJSON has them. The rings around a pole have their vertices on one parallel. Up
to 100 km across, the path of an edge barely matters, so a geodesic, a great circle and geo's
straight line on an equal-area map all enclose almost the same area. The errors are in the size of
the area, whatever its sign.

| Crate | Polygons | Failures | Max error (m²) | Max error (share of the area) |
|---|---|---|---|---|
| geo3d | fields, 10 m to 1 km across | 0 | 0 | 0 % |
| geo3d | regions, 1 km to 100 km across | 0 | 0 | 0 % |
| geo3d | across the antimeridian, 10 km to 100 km | 0 | 0 | 0 % |
| geo3d | around a pole, 60° to 89.9° | 0 | 3.1e-2 | 6.2e-9 % |
| geo (`GeodesicArea`) | fields, 10 m to 1 km across | 0 | 0 | 0 % |
| geo (`GeodesicArea`) | regions, 1 km to 100 km across | 0 | 0 | 0 % |
| geo (`GeodesicArea`) | across the antimeridian, 10 km to 100 km | 0 | 0 | 0 % |
| geo (`GeodesicArea`) | around a pole, 60° to 89.9° | 0 | 0 | 0 % |
| geo (`ChamberlainDuquetteArea`, sphere) | fields, 10 m to 1 km across | 0 | 3.0e3 | 0.67 % |
| geo (`ChamberlainDuquetteArea`, sphere) | regions, 1 km to 100 km across | 0 | 1.8e7 | 0.71 % |
| geo (`ChamberlainDuquetteArea`, sphere) | across the antimeridian, 10 km to 100 km | 0 | 3.0e12 | 9.4e5 % |
| geo (`ChamberlainDuquetteArea`, sphere) | around a pole, 60° to 89.9° | 0 | 3.4e13 | 100.00 % |
| jord (spherical excess, sphere) | fields, 10 m to 1 km across | 0 | 2.7e3 | 1.35 % |
| jord (spherical excess, sphere) | regions, 1 km to 100 km across | 0 | 2.1e7 | 0.74 % |
| jord (spherical excess, sphere) | across the antimeridian, 10 km to 100 km | 0 | 2.1e7 | 0.73 % |
| jord (spherical excess, sphere) | around a pole, 60° to 89.9° | 0 | 2.5e11 | 0.89 % |

</details>

## 🧭 How far off the path?

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/cross-track-dark.svg">
  <img alt="Worst cross-track error on short paths: geo3d 1.6 nm, geo 3.4 nm, jord 73 µm, textbook n-vector formula 2.0 km" src="charts/cross-track-light.svg">
</picture>

We asked each crate how far a point is from a path, on 20 000 paths from a micrometre to 1 000 km
long. The paths run due north along a meridian, so the exact answer is known.

geo3d stays within 2 nanometres everywhere. It builds each path from the latitude and longitude
differences, so short paths lose nothing. geo is within 4 nm. jord is within about 5 nm on long
paths, but drifts to 73 µm on short ones. The textbook n-vector formula, used in the NavLab examples and
by geo3d itself before #7, can be kilometres out when the path is only micrometres long.

<details>
<summary>📋 Show the numbers</summary>

10 000 paths of each kind, at random places between 70° S and 70° N, with a point 1 m to 1 km to
one side of the middle. Each crate is checked on its own sphere radius. "Textbook" is the formula
written out in [`src/bin/accuracy.rs`](src/bin/accuracy.rs): the path's normal is the unit vector
of n1 × n2.

| Crate | Paths | Failures | Max error (m) |
|---|---|---|---|
| geo3d | short paths, 1 µm to 1 km | 0 | 1.6e-9 |
| jord | short paths, 1 µm to 1 km | 0 | 7.3e-5 |
| geo | short paths, 1 µm to 1 km | 0 | 3.4e-9 |
| textbook n1 × n2 (NavLab) | short paths, 1 µm to 1 km | 0 | 2.0e3 |
| geo3d | long paths, 1 km to 1 000 km | 0 | 1.6e-9 |
| jord | long paths, 1 km to 1 000 km | 0 | 5.1e-9 |
| geo | long paths, 1 km to 1 000 km | 0 | 3.8e-9 |
| textbook n1 × n2 (NavLab) | long paths, 1 km to 1 000 km | 0 | 1.8e-6 |

</details>

## 🗺️ UTM grid coordinates

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/utm-dark.svg">
  <img alt="UTM: geo3d within 1.9 nm of GeographicLib and 3.5 nm grid and back; utm crate 750 µm and 77 cm" src="charts/utm-light.svg">
</picture>

We checked UTM against GeographicLib's own GeoConvert at nanometre precision, then compared geo3d
and the utm crate on 30 000 random points: near zones' central meridians, at their edges, and in
Svalbard's double-width zones.

- geo3d uses Krüger's series to the eighth order (GeographicLib uses the sixth). It's within 2 nm of
  GeographicLib, and a round trip comes back within 4 nm.
- The utm crate uses the older USGS formulas. To the grid it's under a millimetre out, but its way
  back is cruder: a round trip can land 77 cm away.
- At the edges of the rules, the utm crate:
  - puts a point exactly on the equator 10 000 km south;
  - calls the 180° meridian zone 61, which doesn't exist;
  - carries on in UTM beyond 84° N instead of switching to UPS.

<details>
<summary>📋 Show the numbers</summary>

| Library | Points | Count | Wrong zone | Max error (m) |
|---|---|---|---|---|
| geo3d | GeographicLib's reference points | 18 | 0 | 1.9e-9 |
| utm | GeographicLib's reference points | 18 | 0 | 7.5e-4 |
| utm | within 1° of the central meridian, against geo3d | 10000 | 0 | 9.2e-4 |
| geo3d | within 1° of the central meridian, grid and back | 10000 | 0 | 2.5e-9 |
| utm | within 1° of the central meridian, grid and back | 10000 | 0 | 3.3e-1 |
| utm | 2° to 3° from the central meridian, against geo3d | 10000 | 0 | 8.4e-4 |
| geo3d | 2° to 3° from the central meridian, grid and back | 10000 | 0 | 3.2e-9 |
| utm | 2° to 3° from the central meridian, grid and back | 10000 | 0 | 7.7e-1 |
| utm | Svalbard's wide zones, 72° to 84° N, against geo3d | 10000 | 0 | 9.6e-4 |
| geo3d | Svalbard's wide zones, 72° to 84° N, grid and back | 10000 | 0 | 3.5e-9 |
| utm | Svalbard's wide zones, 72° to 84° N, grid and back | 10000 | 0 | 2.0e-1 |

At the edges of its rules:

| Case | geo3d | utm |
|---|---|---|
| Exactly on the equator (0°, 3°) | zone 31, 500000 E, 0 N | zone 31, 500000 E, 10000000 N |
| On the 180° meridian (10°, 180°) | zone 1, 171071 E, 1106909 N | zone 61, 171071 E, 1106909 N |
| Beyond 84° N, where UTM stops (85°, 10°) | UPS (not UTM) | zone 32, 509730 E, 9439817 N |

</details>

## ✈️ Yaw, pitch and roll

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/attitude-dark.svg">
  <img alt="Worst attitude error near ±90° pitch, at 1 000 km: geo3d 0.89 nm, sguaba 1900 km, jord via ECEF 2000 km" src="charts/attitude-light.svg">
</picture>

We turned 30 000 attitudes into rotations and back. Away from the vertical, every crate gets them
back. Near ±90° pitch, yaw and roll turn about the same axis, so only the orientation they rebuild
can be checked.

- geo3d rebuilds the same orientation every time, to within a nanometre at 1 000 km. It reads roll
  first, then yaw from the full-size matrix entries given that roll, and pitch with `atan2`.
- sguaba (through nalgebra) reads pitch with `asin`, which loses precision near the vertical, and
  splits yaw and roll from entries that are almost zero. Its orientation can come back pointing
  1 900 km off at 1 000 km.
- jord's own round trip is fine. But once the rotation has been through ECEF, as one from an
  inertial navigation system has, its tiny entries carry ordinary rounding, and its yaw and roll
  can come back pointing 2 000 km off at 1 000 km.

<details>
<summary>📋 Show the numbers</summary>

10 000 attitudes in each set: yaw and roll anywhere, pitch anywhere or within 1e-15° to 0.1° of
±90°. "Via ECEF" turns the attitude into ECEF at a random place and back with plain matrix
products, then reads the angles out (sguaba has no matrix-to-angles function to call directly). The
turn error is how far the rebuilt orientation points from the original, at 1 000 km.

| Crate | Attitudes | Failures | Max pitch error (°) | Max turn error at 1 000 km (m) |
|---|---|---|---|---|
| geo3d | any pitch | 0 | 1.4e-14 | 6.1e-10 |
| sguaba | any pitch | 0 | 1.4e-12 | 2.4e-8 |
| jord | any pitch | 0 | 2.1e-14 | 1.3e-9 |
| geo3d | within 1e-15° to 0.1° of ±90° pitch | 0 | 0 | 5.6e-10 |
| sguaba | within 1e-15° to 0.1° of ±90° pitch | 0 | 1.9e-6 | 1.9e6 |
| jord | within 1e-15° to 0.1° of ±90° pitch | 0 | 1.4e-14 | 1.4e-9 |
| geo3d | via ECEF, within 1e-15° to 0.1° of ±90° pitch | 0 | 1.4e-14 | 8.9e-10 |
| jord | via ECEF, within 1e-15° to 0.1° of ±90° pitch | 0 | 1.4e-14 | 2.0e6 |

</details>

## 🖨️ Printing a position

Only geo3d and sguaba can print a position. geo3d prints the number you gave it: −33.87 comes out
as `-33.87°`, and `{:.3}` gives three decimals. sguaba keeps angles in radians, so the same latitude
comes back as `33.869999999999976°S`, and it ignores `{:.3}`.

<details>
<summary>📋 Show the numbers</summary>

The same position, typed in and after a trip through ECEF and back with geo3d (as a computed one
would be), printed plainly and with three decimals.

| Position and format | geo3d | sguaba |
|---|---|---|
| typed in, `{}` | `-33.87°, 151.21°, 20 m` | `33.869999999999976°S, 151.21°E, 20m` |
| typed in, `{:.3}` | `-33.870°, 151.210°, 20 m` | `33.869999999999976°S, 151.21°E, 20m` |
| computed, `{}` | `-33.87°, 151.21°, 19.99999999992662 m` | `33.869999999999976°S, 151.21°E, 19.99999999992662m` |
| computed, `{:.3}` | `-33.870°, 151.210°, 20 m` | `33.869999999999976°S, 151.21°E, 19.99999999992662m` |

</details>

## 📡 Look angles

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/look-angles-dark.svg">
  <img alt="Worst look-angle error on hard cases: geo3d 25 nm, jord 13 cm, sguaba 13 cm" src="charts/look-angles-light.svg">
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
| geo3d | 0 | 0.002 | 0.004 | 1.7e-11 |
| map_3d | 0.002 | 0 | 0.004 | 1.7e-11 |
| nav-types | 0.001 | 0.002 | — | — |
| sguaba | 0.001 | 0.002 | 1 745 | 1.0e-7 |
| jord | 0.004 | 0.004 | 0.004 | 2.6e-11 |

2 000 pairs straight up and 1e-7° to 1e-3° off it, from 1 km to geostationary range.

| Crate | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing error (°) |
|---|---|---|---|---|
| geo3d | 0 | 0.015 | 0.025 | 1.6e-10 |
| map_3d | 0.015 | 0 | 0.025 | 1.6e-10 |
| nav-types | 0.015 | 0.019 | — | — |
| sguaba | 0.026 | 0.028 | 129 550 | 4.8e-7 |
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
  <img alt="Feature matrix: geo3d is the only crate with ground strike, 3D volumes, closest approach, covariance and datum shifts" src="charts/features-light.svg">
</picture>

geo3d covers everything from a GPS fix to a decision: where to point, where a line of sight lands,
whether something is inside a zone, when two objects come closest, how sure you can be, and how to
shift datums. It does
this with no dependencies and no `unsafe` code.

<details>
<summary>📋 Show the numbers</summary>

Checked in each crate's source at the versions above. Dependencies are counted with default
features, from `cargo tree -e normal`.

| | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo | utm |
|---|---|---|---|---|---|---|---|---|
| Direct dependencies | 0 | 0 | 1 | 5 | 0 | 2 | 13 | 0 |
| All crates pulled in | 0 | 0 | 19 | 21 | 0 | 5 | 41 | 0 |
| `unsafe` | `#![forbid(unsafe_code)]` | none, not forbidden | none, not forbidden | marks frame assumptions: users write `unsafe` to build ECEF ↔ ENU/NED transforms; no memory-unsafe code | `#![deny(unsafe_code)]` | none, not forbidden | 2 blocks | none, not forbidden |
| MSRV | 1.85 | not declared | not declared | 1.89 | 1.85 | 1.70 | 1.88 | not declared |
| Without std (no_std) | yes, with `libm` | no | no | yes, with `libm` | no | no | no | yes, with its `no_std` feature |
| Ellipsoids | 8 named, or any | 10 named | WGS84 | WGS84 | 4 named, or any; spheres | any | WGS84 (`Geodesic`); spheres | WGS84 |
| Look angles | yes | yes | no (ENU/NED only) | yes | yes | — | — | — |
| Ray–ellipsoid ground strike | yes | no | no | no | no | — | — | — |
| 3D volumes | sphere, cone, cylinder, altitude band | no | no | no | no (caps and polygons on a sphere) | — | no (2D shapes) | — |
| Closest approach | 3D, constant velocity | no | no | no | on a sphere (surface tracks) | — | no | — |
| Datum shift (Helmert) | yes | no | no | no | no | — | no | — |
| Inertial frame | TEME (GMST rotation) | GMST rotation (`ecef2eci`) | no | no | no | — | no | — |
| Geodesics | Karney | haversine only | haversine only | haversine only | great circles on a sphere | Karney | Karney (geographiclib-rs) | — |
| Polygon area | geodesic | no | no | no | on a sphere | geodesic | geodesic (geographiclib-rs); on a sphere | — |
| Covariance through frames | ECEF, ENU, NED, AER, geodetic | no | no | no | no | — | no | — |
| n-vector type | no (used inside) | no | yes | no | yes | no | no | — |
| Route tools | cross-track, along-track, crossings, centre, interpolation | no | no | no | yes | no | cross-track, interpolation | — |
| Attitude (yaw, pitch, roll) | yes, exact at ±90° pitch | no | no | yes | yes | no | no | — |
| UTM/UPS grid | UTM and UPS, Krüger to the 8th order | — | — | — | — | — | — | UTM only (USGS series) |
| Prints positions | yes: lossless, `{:.N}`, `{:#}` for N/S/E/W | no | no | yes, ignoring `{:.N}` | no | no | no | no |

</details>

## ⚡ Speed

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="charts/speed-dark.svg">
  <img alt="Time per call for 13 operations; geo3d is fastest in 3 and tied in 3" src="charts/speed-light.svg">
</picture>

geo3d is the fastest, or level with the fastest, in six of the thirteen tests. It's slower in the
other seven, and we know why:

- **x, y, z back to latitude and longitude** (two tests): about 10 % behind map_3d and jord, and
  further behind nav-types and sguaba. geo3d spends that time on a final check that makes it the
  most accurate of the five, at every height.
- **One-off ENU conversions:** a `LocalFrame` costs a little to build, and pays off from the second
  target onwards.
- **Quick great-circle distance:** geo3d uses a formula that stays exact on the far side of the
  Earth. The faster formula most crates use can be almost a metre out there.
- **Distance along the ground:** about 4 % behind geographiclib-rs. geo is faster again because it
  skips the compass bearings.
- **UTM, both ways:** the utm crate is quicker (1.5× to the grid, 4× back). geo3d finds the
  latitude exactly on the way back, where the utm crate's formulas can be 77 cm out.

<details>
<summary>📋 Show the numbers</summary>

Measured on:

- AMD Ryzen 9 5950X (16 cores, 32 threads), Ubuntu 24.04.5 LTS, Linux 7.0.0
- rustc 1.98.0 (88d9e12ae 2026-08-18), default `x86_64-unknown-linux-gnu` target
- criterion 0.5.1: 300 ms warm-up, 1 s measurement, 50 samples

Median nanoseconds per call. In bold: the fastest, and anything within 3 % of it.

| Operation | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo | utm |
|---|---|---|---|---|---|---|---|---|
| Geodetic → ECEF | **21** | 30 | 36 | 53 | 28 | — | — | — |
| ECEF → geodetic, −500 m to 12 km | 119 | 109 | 87 | **69** | 104 | — | — | — |
| ECEF → geodetic, 400 km to 35 786 km | 121 | 108 | **85** | 98 | 102 | — | — | — |
| ENU, one observer (frame reused) | **28** | 49 | 146 | 59 | 34 | — | — | — |
| ENU, single call (frame built each time) | 53 | **49** | 185 | 132 | 90 | — | — | — |
| Look angles, one observer (frame reused) | **67** | 94 | — | 98 | **67** | — | — | — |
| Look angles, single call | **91** | **93** | — | 170 | 125 | — | — | — |
| Geodesic inverse | 861 | — | — | — | — | 830 | **768 (distance only)** | — |
| Geodesic direct | **275** | — | — | — | — | 316 | 286 (position only) | — |
| Great-circle distance | 60 | 46 | **38** | 53 | 67 | — | **38** | — |
| Polygon area, 10 vertices | **4 985** | — | — | — | — | **4 960** | 5 174 | — |
| Latitude/longitude → UTM | 198 | — | — | — | — | — | — | **129** |
| UTM → latitude/longitude | 300 | — | — | — | — | — | — | **75** |

Inputs: 1 000 points anywhere on Earth from −500 m to 12 km high (and from 400 km to 35 786 km for
the second ECEF → geodetic row); a radar near Sydney with 1 000 aircraft within a degree of it, up
to 12 km high; 1 000 random point pairs; 1 000 lines on random bearings up to 20 000 km long; 1 000
polygons of 10 vertices, 1 km to 100 km across; 1 000 points between 80° S and 84° N for UTM.
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
| [utm](https://crates.io/crates/utm) | 0.1.6 | UTM coordinates on WGS84, by the USGS series |
| [geo](https://crates.io/crates/geo) | 0.33.1 | 2D geometry; only its `Haversine` and `Geodesic` distances and its `GeodesicArea` and `ChamberlainDuquetteArea` areas and its `CrossTrackDistance` are used |

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
