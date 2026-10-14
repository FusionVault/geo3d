## Round trip: geodetic → ECEF → geodetic

66960 points per height: a 1° grid plus latitudes ±89.9°, ±89.99°, ±89.999°. Each library uses its own forward and inverse conversion.

### Max latitude error (degrees)

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| -11 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 3.6e-14 | 2.1e-14 |
| 0 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 10 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 400 km | 2.1e-14 | 2.5e-9 | 1.4e-14 | 4.3e-14 | 2.1e-14 |
| 20200 km | 1.4e-14 | 7.0e-5 | 2.1e-14 | 5.0e-14 | 2.1e-14 |
| 35786 km | 1.8e-14 | 1.1e-4 | 1.4e-14 | 4.3e-14 | 2.1e-14 |

### Max longitude error (degrees)

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| -11 km | 1.4e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 0 km | 1.4e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.1e-14 |
| 10 km | 1.4e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 400 km | 1.4e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 20200 km | 1.4e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 35786 km | 1.4e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |

### Max height error (metres)

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| -11 km | 1.9e-9 | 2.0e-9 | 3.0e-9 | 1.9e-7 | 2.7e-9 |
| 0 km | 1.9e-9 | 2.1e-9 | 1.4e-9 | 1.2e-9 | 2.1e-9 |
| 10 km | 1.9e-9 | 2.2e-9 | 2.5e-9 | 2.6e-9 | 3.6e-9 |
| 400 km | 2.8e-9 | 2.2e-9 | 2.7e-9 | 3.2e-7 | 3.2e-9 |
| 20200 km | 1.1e-8 | 6.3e-6 | 1.1e-8 | 1.2e-6 | 1.1e-8 |
| 35786 km | 1.5e-8 | 1.4e-5 | 1.5e-8 | 2.2e-6 | 1.5e-8 |

## Geodesics, against geographiclib-rs

100000 random point pairs (uniform over the sphere) and 10000 nearly antipodal pairs (the second point within 0.5° of the first one's antipode). A failure is `None`, NaN or a panic. geo's `Geodesic` calls geographiclib-rs, so it matches it exactly.

### Inverse: distance and azimuths between two points

| Library | Pairs | Count | Failures | Max distance error (m) | Mean distance error (m) | Max azimuth error (°) |
|---|---|---|---|---|---|---|
| geo3d | random | 100000 | 0 | 3.7e-9 | 4.7e-13 | 5.7e-14 |
| geo (`Geodesic`) | random | 100000 | 0 | 0 | 0 | — |
| geo3d | nearly antipodal | 10000 | 0 | 3.7e-9 | 3.7e-13 | 4.5e-13 |
| geo (`Geodesic`) | nearly antipodal | 10000 | 0 | 0 | 0 | — |

### Direct: travel the reference distance on the reference bearing

| Library | Lines | Count | Failures | Max position error (m) | Mean position error (m) | Max final azimuth error (°) |
|---|---|---|---|---|---|---|
| geo3d | random | 100000 | 0 | 0 | 0 | 0 |
| geo (`Geodesic`) | random | 100000 | 0 | 0 | 0 | — |
| geo3d | nearly antipodal | 10000 | 0 | 0 | 0 | 0 |
| geo (`Geodesic`) | nearly antipodal | 10000 | 0 | 0 | 0 | — |

### Great-circle (spherical) distance against the ellipsoidal geodesic, all pairs

| Library | Sphere radius | Max error | Mean error | Max error (km) | Failures |
|---|---|---|---|---|---|
| geo3d | 6371008.771 m | 0.56 % | 0.13 % | 37.8 km | 0 |
| map_3d | 6371000 m | 0.56 % | 0.13 % | 37.8 km | 0 |
| nav-types | 6378137 m | 0.67 % | 0.16 % | 50.7 km | 0 |
| sguaba | 6378137 m | 0.67 % | 0.16 % | 50.7 km | 0 |
| jord | 6371000.8 m | 0.56 % | 0.13 % | 37.8 km | 0 |
| geo | 6371008.8 m | 0.56 % | 0.13 % | 37.8 km | 0 |

## Local frames and look angles

Every library is compared with geo3d and with map_3d. For targets built from a known pointing, that pointing is the reference too. Disagreements are in micrometres. For look angles, each answer (azimuth, elevation, range) is turned back into a point and the distance between the two points is reported, so a pointing error counts by how far it moves the target.

### 15000 observer–target pairs: nearby aircraft, and targets anywhere up to geostationary height

| Library | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs geo3d (µm) | Look angles vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing vs the built pointing (°) |
|---|---|---|---|---|---|---|
| geo3d | 0.000 | 0.017 | 0.000 | 0.077 | — | — |
| map_3d | 0.017 | 0.000 | 0.077 | 0.000 | — | — |
| nav-types | 0.018 | 0.017 | — | — | — | — |
| sguaba | 0.037 | 0.038 | 8.237 | 8.237 | — | — |
| jord | 0.022 | 0.020 | 0.081 | 0.060 | — | — |

### 2400 pairs: on the horizon, on and just off (1e-7°, 1e-5°) due north, east, south and west, at 10 km and 1 000 km

| Library | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs geo3d (µm) | Look angles vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing vs the built pointing (°) |
|---|---|---|---|---|---|---|
| geo3d | 0.000 | 0.002 | 0.000 | 0.002 | 0.006 | 3.3e-11 |
| map_3d | 0.002 | 0.000 | 0.002 | 0.000 | 0.006 | 3.7e-11 |
| nav-types | 0.001 | 0.002 | — | — | — | — |
| sguaba | 0.001 | 0.002 | 1745 | 1745 | 1745 | 1.0e-7 |
| jord | 0.004 | 0.004 | 0.004 | 0.004 | 0.006 | 3.4e-11 |

### 2000 pairs: straight up and just off it (1e-7° to 1e-3°), 1 km to geostationary range

| Library | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs geo3d (µm) | Look angles vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing vs the built pointing (°) |
|---|---|---|---|---|---|---|
| geo3d | 0.000 | 0.015 | 0.000 | 0.016 | 0.030 | 2.1e-10 |
| map_3d | 0.015 | 0.000 | 0.016 | 0.000 | 0.030 | 2.1e-10 |
| nav-types | 0.015 | 0.017 | — | — | — | — |
| sguaba | 0.027 | 0.030 | 129550 | 129550 | 129550 | 4.8e-7 |
| jord | 0.023 | 0.023 | 129550 | 129550 | 129550 | 4.8e-7 |

## Robustness: edge cases

"ok" means within 1 mm and 1e-8°. For NaN input, NaN out or a rejection (`None`) are both fine; a panic is not.

### What each library returns

| Case | geo3d | map_3d | nav-types | sguaba | jord | geographiclib-rs | geo |
|---|---|---|---|---|---|---|---|
| ECEF → geodetic at the Earth's centre (0, 0, 0); want ±90°, h = −b | ok | NaN | partly NaN (NaN°, 0.000000°, -6292741.6546 m) | partly NaN (90.000000°, NaN°, -6356752.3142 m) | NaN | — | — |
| ECEF → geodetic exactly above the North Pole (0, 0, b + 100 m) | ok | ok | partly NaN (NaN°, 0.000000°, NaN m) | ok | ok | — | — |
| Round trip exactly at the South Pole (−90°, 0°, 100 m) | ok | ok | partly NaN (NaN°, 0.000000°, NaN m) | ok | ok | — | — |
| ENU across the antimeridian, 179.9999° → −179.9999° (want 22.264 m east) | ok | ok | ok | ok | ok | — | — |
| Round trip at 1e8 m altitude (45°, 45°, 100 000 km) | ok | wrong (45.000126°, 45.000000°, 100000000.0000 m) | ok | ok | ok | — | — |
| Geodetic → ECEF with a NaN latitude | NaN | NaN | rejected (`None`) | rejected (`None`) | NaN | — | — |
| ECEF → geodetic with a NaN x | NaN | NaN | NaN | panic | NaN | — | — |
| Geodesic between exact antipodes (0°, 0°) → (0°, 180°); want 20 003 931.459 m | ok | — | — | — | — | ok | ok |
| Geodesic with a NaN latitude | rejected (`None`) | — | — | — | — | NaN | NaN |

