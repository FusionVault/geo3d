## Round trip: geodetic → ECEF → geodetic

66960 points per height: a 1° grid plus latitudes ±89.9°, ±89.99°, ±89.999°. Each library uses its own forward and inverse conversion.

### Max latitude error (degrees)

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| -11 km | 1.4e-14 | 2.1e-14 | 2.1e-14 | 3.6e-14 | 2.1e-14 |
| 0 km | 1.4e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 10 km | 2.1e-14 | 2.1e-14 | 2.1e-14 | 4.3e-14 | 2.1e-14 |
| 400 km | 1.4e-14 | 2.5e-9 | 1.4e-14 | 4.3e-14 | 2.1e-14 |
| 20200 km | 1.4e-14 | 7.0e-5 | 2.1e-14 | 5.0e-14 | 2.1e-14 |
| 35786 km | 2.1e-14 | 1.1e-4 | 1.4e-14 | 4.3e-14 | 2.1e-14 |

### Max longitude error (degrees)

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| -11 km | 2.8e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 0 km | 2.8e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.1e-14 |
| 10 km | 2.8e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 400 km | 2.8e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 20200 km | 2.8e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |
| 35786 km | 2.8e-14 | 1.4e-14 | 1.4e-14 | 5.7e-14 | 2.8e-14 |

### Max height error (metres)

| Height | geo3d | map_3d | nav-types | sguaba | jord |
|---|---|---|---|---|---|
| -11 km | 2.1e-9 | 2.0e-9 | 3.0e-9 | 1.9e-7 | 2.7e-9 |
| 0 km | 2.1e-9 | 2.1e-9 | 1.4e-9 | 1.2e-9 | 2.1e-9 |
| 10 km | 2.1e-9 | 2.2e-9 | 2.5e-9 | 2.6e-9 | 3.6e-9 |
| 400 km | 2.0e-9 | 2.2e-9 | 2.7e-9 | 3.2e-7 | 3.2e-9 |
| 20200 km | 7.5e-9 | 6.3e-6 | 1.1e-8 | 1.2e-6 | 1.1e-8 |
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

## Polygon areas, against geographiclib-rs

Simple polygons of 3 to 32 vertices, counter-clockwise, with longitudes in [-180°, 180°) as GeoJSON has them. The rings round a pole have vertices on one parallel. Errors are in the size of the area, whatever its sign. A failure is `None`, NaN or a panic. geo's `GeodesicArea` calls geographiclib-rs, so it matches it exactly.

### Area of a polygon

| Library | Polygons | Count | Failures | Max error (m²) | Max error (share of the area) |
|---|---|---|---|---|---|
| geo3d | fields, 10 m to 1 km across | 1000 | 0 | 0 | 0 % |
| geo3d | regions, 1 km to 100 km across | 1000 | 0 | 0 | 0 % |
| geo3d | across the antimeridian, 10 km to 100 km | 1000 | 0 | 0 | 0 % |
| geo3d | around a pole, 60° to 89.9° | 1000 | 0 | 1.2e-1 | 9.5e-9 % |
| geo (`GeodesicArea`) | fields, 10 m to 1 km across | 1000 | 0 | 0 | 0 % |
| geo (`GeodesicArea`) | regions, 1 km to 100 km across | 1000 | 0 | 0 | 0 % |
| geo (`GeodesicArea`) | across the antimeridian, 10 km to 100 km | 1000 | 0 | 0 | 0 % |
| geo (`GeodesicArea`) | around a pole, 60° to 89.9° | 1000 | 0 | 0 | 0 % |
| geo (`ChamberlainDuquetteArea`, sphere) | fields, 10 m to 1 km across | 1000 | 0 | 3.0e3 | 0.67 % |
| geo (`ChamberlainDuquetteArea`, sphere) | regions, 1 km to 100 km across | 1000 | 0 | 1.8e7 | 0.71 % |
| geo (`ChamberlainDuquetteArea`, sphere) | across the antimeridian, 10 km to 100 km | 1000 | 0 | 3.0e12 | 9.4e5 % |
| geo (`ChamberlainDuquetteArea`, sphere) | around a pole, 60° to 89.9° | 1000 | 0 | 3.4e13 | 100.00 % |
| jord (spherical excess, sphere) | fields, 10 m to 1 km across | 1000 | 0 | 2.7e3 | 1.35 % |
| jord (spherical excess, sphere) | regions, 1 km to 100 km across | 1000 | 0 | 2.1e7 | 0.74 % |
| jord (spherical excess, sphere) | across the antimeridian, 10 km to 100 km | 1000 | 0 | 2.1e7 | 0.73 % |
| jord (spherical excess, sphere) | around a pole, 60° to 89.9° | 1000 | 0 | 2.5e11 | 0.89 % |

## Cross-track distance, against the exact answer

Paths due north along a meridian at random places, and a point 1 m to 1 km to one side of the middle. On a sphere the exact distance is R·asin(cos φ·sin Δλ); each crate is checked on its own sphere radius. A failure is `None`, NaN or a panic.

### Distance from a point to a path

| Library | Paths | Count | Failures | Max error (m) |
|---|---|---|---|---|
| geo3d | short paths, 1 µm to 1 km | 10000 | 0 | 1.6e-9 |
| jord | short paths, 1 µm to 1 km | 10000 | 0 | 7.3e-5 |
| geo | short paths, 1 µm to 1 km | 10000 | 0 | 3.4e-9 |
| textbook n1 × n2 (NavLab) | short paths, 1 µm to 1 km | 10000 | 0 | 2.0e3 |
| geo3d | long paths, 1 km to 1 000 km | 10000 | 0 | 1.6e-9 |
| jord | long paths, 1 km to 1 000 km | 10000 | 0 | 5.1e-9 |
| geo | long paths, 1 km to 1 000 km | 10000 | 0 | 3.8e-9 |
| textbook n1 × n2 (NavLab) | long paths, 1 km to 1 000 km | 10000 | 0 | 1.8e-6 |

## Attitude: yaw, pitch and roll there and back

Each crate turns yaw, pitch and roll into its own rotation and back. "Via ECEF" instead reads the angles out of an attitude that has been turned into ECEF and back, as one from an inertial navigation system would be, so every entry carries ordinary rounding (sguaba has no matrix-to-angles function to call directly). The pitch error is in degrees. The turn error is how far the rebuilt attitude points from the original, as a distance at 1 000 km: near ±90° pitch, yaw and roll can't be told apart, so only the rotation they rebuild can be checked.

### Yaw, pitch and roll → rotation → yaw, pitch and roll

| Library | Attitudes | Count | Failures | Max pitch error (°) | Max turn error at 1 000 km (m) |
|---|---|---|---|---|---|
| geo3d | any pitch | 10000 | 0 | 1.4e-14 | 6.1e-10 |
| sguaba | any pitch | 10000 | 0 | 1.4e-12 | 2.4e-8 |
| jord | any pitch | 10000 | 0 | 2.1e-14 | 1.3e-9 |
| geo3d | within 1e-15° to 0.1° of ±90° pitch | 10000 | 0 | 0 | 5.6e-10 |
| sguaba | within 1e-15° to 0.1° of ±90° pitch | 10000 | 0 | 1.9e-6 | 1.9e6 |
| jord | within 1e-15° to 0.1° of ±90° pitch | 10000 | 0 | 1.4e-14 | 1.4e-9 |
| geo3d | via ECEF, within 1e-15° to 0.1° of ±90° pitch | 10000 | 0 | 1.4e-14 | 8.9e-10 |
| jord | via ECEF, within 1e-15° to 0.1° of ±90° pitch | 10000 | 0 | 1.4e-14 | 2.0e6 |

## UTM grid coordinates

geo3d and the `utm` crate, against GeographicLib's GeoConvert at nanometre precision, and against each other on 10 000 random points per set (geo3d agrees with GeographicLib to a few nanometres). "Grid and back" is each crate's own round trip, as a distance on the ground. A wrong zone is a different zone from GeographicLib's, or a refusal.

### UTM easting and northing

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

### UTM at the edges of its rules

| Case | geo3d | utm |
|---|---|---|
| Exactly on the equator (0°, 3°) | zone 31, 500000 E, 0 N | zone 31, 500000 E, 10000000 N |
| On the 180° meridian (10°, 180°) | zone 1, 171071 E, 1106909 N | zone 61, 171071 E, 1106909 N |
| Beyond 84° N, where UTM stops (85°, 10°) | UPS (not UTM) | zone 32, 509730 E, 9439817 N |

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
| geo3d | 0.000 | 0.002 | 0.000 | 0.002 | 0.004 | 1.7e-11 |
| map_3d | 0.002 | 0.000 | 0.002 | 0.000 | 0.004 | 1.7e-11 |
| nav-types | 0.001 | 0.002 | — | — | — | — |
| sguaba | 0.001 | 0.002 | 1745 | 1745 | 1745 | 1.0e-7 |
| jord | 0.004 | 0.004 | 0.004 | 0.004 | 0.004 | 2.6e-11 |

### 2000 pairs: straight up and just off it (1e-7° to 1e-3°), 1 km to geostationary range

| Library | ENU vs geo3d (µm) | ENU vs map_3d (µm) | Look angles vs geo3d (µm) | Look angles vs map_3d (µm) | Look angles vs the built pointing (µm) | Pointing vs the built pointing (°) |
|---|---|---|---|---|---|---|
| geo3d | 0.000 | 0.015 | 0.000 | 0.016 | 0.025 | 1.6e-10 |
| map_3d | 0.015 | 0.000 | 0.016 | 0.000 | 0.025 | 1.6e-10 |
| nav-types | 0.015 | 0.019 | — | — | — | — |
| sguaba | 0.026 | 0.028 | 129550 | 129550 | 129550 | 4.8e-7 |
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

