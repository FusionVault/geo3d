# geo3d

[![crates.io](https://img.shields.io/crates/v/geo3d.svg)](https://crates.io/crates/geo3d)
[![docs.rs](https://docs.rs/geo3d/badge.svg)](https://docs.rs/geo3d)
[![ci](https://github.com/FusionVault/geo3d/actions/workflows/ci.yml/badge.svg)](https://github.com/FusionVault/geo3d/actions/workflows/ci.yml)
[![MSRV](https://img.shields.io/crates/msrv/geo3d.svg)](Cargo.toml)

**3D positions on Earth, made simple.**

Turn GPS coordinates into 3D space and answer everyday questions: *Where is that plane relative to
me? Where do I point the antenna? How far apart are they? Can they see each other? Will they
collide?* Pure Rust, no dependencies.

```rust
use geo3d::{Geodetic, LocalFrame};

let me = Geodetic::new(51.5, -0.12, 30.0);          // latitude°, longitude°, height in metres
let plane = Geodetic::new(51.6, -0.12, 3_000.0);

let look = LocalFrame::new(me).look_angles(plane);   // where to point
// azimuth 0° (due north), elevation 14.9°, range 11 518 m
assert!(look.elevation_deg > 14.0 && look.range_m > 11_000.0);
```

```bash
cargo add geo3d
```

## Why geo3d?

- **Accurate.** Exact conversions, and geodesics accurate to 15 nanometres. Every function states
  its accuracy, and the tests pin it to GeographicLib and two independent implementations.
- **Reliable.** Geodesics work for every pair of points, including ones on opposite sides of the
  Earth, where Vincenty-based libraries fail. No `unsafe`, no dependencies, and every example in this
  README runs as a test.
- **Simple.** Plain, typed coordinates with named fields. Degrees and metres everywhere, latitude
  always first, so nothing gets mixed up.
- **Fast.** Nothing allocates. A `LocalFrame` does its trigonometry once; after that, each ECEF
  point costs a few multiplications. Rust 1.85+, with optional bridges to `glam`, `nalgebra`, `mint`,
  `geo-types` and `serde`.

## The coordinate types

Latitude and longitude are angles, so you cannot subtract two positions and get a distance or a
direction. geo3d converts them into frames where ordinary 3D maths works.

| Type | What it is | Units |
|---|---|---|
| `Geodetic` | latitude, longitude, height — what GPS gives you | degrees, metres |
| `Ecef` | x, y, z from the Earth's centre, turning with the Earth | metres |
| `Enu` / `Ned` | east-north-up / north-east-down from an observer | metres |
| `Aer` | azimuth, elevation, range — where to point | degrees, metres |
| `Eci` | x, y, z from the Earth's centre, not turning — for satellites (TEME) | metres |
| `Vec3` | a plain 3D vector — a direction or a velocity | metres, m/s |

Conventions, the same everywhere:

- **Degrees and metres.** Latitude comes first: `Geodetic::new(lat, lon, height)`.
- **Azimuths and bearings** run clockwise from north, 0° to 360°. Longitudes come out between −180°
  and 180°.
- **Heights are above the WGS84 ellipsoid** (what GPS reports), not above sea level.
- **WGS84 is the default.** Every calculation is also a method on `Ellipsoid`, so `GRS80` and the
  others work the same way.

New to the terms? See the [glossary](#glossary).

## Recipes

Every example below is compiled and run as a test, so you can copy it.

### Convert latitude/longitude/height to x, y, z

```rust
use geo3d::{ecef_to_geodetic, geodetic_to_ecef, Geodetic, GRS80};

let sydney = Geodetic::new(-33.87, 151.21, 20.0);

let xyz = geodetic_to_ecef(sydney);    // x −4 646 032 m, y 2 553 122 m, z −3 534 494 m
let back = ecef_to_geodetic(xyz);      // Sydney again
assert!((back.lat_deg - sydney.lat_deg).abs() < 1e-11);

let on_grs80 = GRS80.to_ecef(sydney);  // any other ellipsoid works the same way
```

### Where is it relative to me?

```rust
use geo3d::{Geodetic, LocalFrame};

let radar = LocalFrame::new(Geodetic::new(51.5, -0.12, 30.0));
let target = Geodetic::new(51.6, -0.12, 3_000.0);

let enu = radar.enu_of(target);        // east 0 m, north 11 131 m, up 2 960 m
let ned = radar.ned_of(target);        // the same as north-east-down
let back = radar.geodetic_of(enu);     // and back to latitude/longitude/height
assert!(enu.north > 11_000.0 && ned.down < 0.0);
```

Make one `LocalFrame` per observer and reuse it: building it is the expensive part.

### Where do I point?

```rust
use geo3d::{Geodetic, LocalFrame};

let station = LocalFrame::new(Geodetic::new(-33.94, 151.18, 21.0));
let aircraft = Geodetic::new(-33.80, 151.30, 2_500.0);

let aer = station.look_angles(aircraft);  // azimuth 35.6°, elevation 7.3°, range 19 254 m
let again = station.destination(aer);     // and back: the aircraft's position
assert!(aer.is_above_horizon());
```

### Where does my line of sight hit the ground?

```rust
use geo3d::{ground_strike, Geodetic};

let drone = Geodetic::new(-33.94, 151.18, 120.0);    // 120 m up

let spot = ground_strike(drone, 90.0, -30.0);        // looking east, 30° down: 208 m east
assert!(spot.is_some());
assert!(ground_strike(drone, 90.0, 10.0).is_none()); // looking up: never hits
```

### How far apart?

There are three different distances. Pick the one you mean:

| Distance | Function | Accuracy |
|---|---|---|
| Straight line through space (slant range) | `distance3d_m(a, b)` | exact |
| Along the surface | `WGS84.geodesic_inverse(a, b)` | about 15 nm |
| Along the surface, quick estimate | `WGS84.great_circle_distance_m(a, b)` | within 0.6 % |

```rust
use geo3d::{distance3d_m, Geodetic, WGS84};

let lax = Geodetic::new(33.9425, -118.4081, 0.0);
let jfk = Geodetic::new(40.6413, -73.7781, 0.0);

let straight = distance3d_m(lax, jfk);                    // 3 918 745 m, through the Earth
let surface = WGS84.geodesic_inverse(lax, jfk).unwrap();  // 3 983 005 m, starting on bearing 65.9°
let quick = WGS84.great_circle_distance_m(lax, jfk);      // 3 974 267 m

// The other way round: go that far on that bearing and you arrive at JFK.
let (arrive, _) = WGS84.geodesic_direct(lax, surface.initial_bearing_deg, surface.distance_m);
assert!((arrive.lat_deg - jfk.lat_deg).abs() < 1e-8);
```

The surface distance works for every pair of points, including ones on opposite sides of the Earth
(where the older Vincenty method used by many libraries fails).

### Is it inside?

Every shape implements `Volume`, so they all answer `contains`:

```rust
use geo3d::{AltitudeBand, Cone, Cylinder, Geodetic, Sphere, Volume};

let site = Geodetic::new(0.0, 0.0, 0.0);
let p = Geodetic::new(0.0, 0.02, 500.0);  // 2.2 km east, 500 m up

let bubble = Sphere::new(site, 10_000.0);                             // within 10 km
let camera = Cone::new(site, 90.0, 0.0, 20.0, 5_000.0);               // looking east, ±20°, 5 km reach
let zone = Cylinder::new(site, 8_000.0, AltitudeBand::below(2_000.0)); // 8 km radius, up to 2 km

assert!(bubble.contains(p) && camera.contains(p) && zone.contains(p));
```

`AltitudeBand` adds a height range to any 2D footprint: `band.within_prism(inside_polygon, height)`.

### Will they collide?

```rust
use geo3d::{Course, Geodetic, LocalFrame, Track, WGS84};

let a_pos = Geodetic::new(-33.80, 151.30, 2_500.0);
let b_pos = Geodetic::new(-33.90, 151.25, 2_400.0);

// Heading°, ground speed m/s, climb m/s → a 3D velocity at each position.
let a_vel = LocalFrame::new(a_pos).velocity_of_course(Course::new(225.0, 120.0, 0.0)); // south-west
let b_vel = LocalFrame::new(b_pos).velocity_of_course(Course::new(0.0, 100.0, 0.0));   // north

let a = Track::new(WGS84.to_ecef(a_pos), a_vel);
let b = Track::new(WGS84.to_ecef(b_pos), b_vel);

let cpa = a.closest_approach(b);      // closest in 59 s, 434 m apart
let alert = a.conflict(b, 1_000.0);   // Some(cpa): they come within 1 km
let soon = a.at(60.0);                // where a will be in a minute
assert!(alert.is_some());
```

### Can they see each other?

```rust
use geo3d::{Geodetic, WGS84};

let plane = WGS84.to_ecef(Geodetic::new(0.0, 0.0, 10_000.0));  // 10 km up
let near = WGS84.to_ecef(Geodetic::new(0.0, 1.0, 10_000.0));   // 111 km away
let far = WGS84.to_ecef(Geodetic::new(0.0, 10.0, 10_000.0));   // 1 113 km away

assert!(WGS84.line_of_sight(plane, near));
assert!(!WGS84.line_of_sight(plane, far));       // the Earth is in the way

let horizon = WGS84.horizon_distance_m(10_000.0); // 357 km
```

The Earth here is a smooth ellipsoid: no terrain, no atmospheric refraction.

### Satellites: the inertial frame

```rust
use geo3d::{eci, Ecef};

let angle = eci::gmst_from_unix(1_800_000_000.0);  // how far the Earth has turned at that moment

let fixed = Ecef::new(7_000_000.0, 0.0, 0.0);
let inertial = fixed.to_eci(angle);
let back = inertial.to_ecef(angle);
assert!((back.x - fixed.x).abs() < 1e-6);
```

This inertial frame is **TEME**, the one SGP4 / TLE propagators output, so their positions plug
straight in. It is not J2000 / GCRF, which differs by about 0.4°. Unix time is UTC, which can be up
to 0.9 s off the UT1 the Earth's rotation needs (up to about 400 m at the equator); pass a UT1
Julian date to `eci::gmst_rad` for better.

### Shift between datums

```rust
use geo3d::{Ecef, Helmert7};

// WGS84 → OSGB36 (Great Britain), as published: shifts in metres, rotations in arc-seconds, scale in ppm.
let h = Helmert7::position_vector(-446.448, 125.157, -542.060, -0.1502, -0.2470, -0.8421, 20.4894);

let p = Ecef::new(3_874_938.0, 116_218.0, 5_047_168.0);
let shifted = h.apply(p);
let back = h.apply_inverse(shifted);  // the exact reverse
assert!(back.distance_to(p) < 1e-6);
```

Publishers use one of two rotation-sign conventions. For the other one, use
`Helmert7::coordinate_frame` with the same arguments.

## Coming from MATLAB, pymap3d or map_3d?

The same operations, by their familiar names. Build the `LocalFrame` once per observer:
`let frame = LocalFrame::new(observer);`

| MATLAB / pymap3d / map_3d | geo3d |
|---|---|
| `geodetic2ecef` / `ecef2geodetic` | `geodetic_to_ecef(g)` / `ecef_to_geodetic(e)` |
| `geodetic2enu` / `enu2geodetic` | `frame.enu_of(g)` / `frame.geodetic_of(enu)` |
| `geodetic2ned` / `ned2geodetic` | `frame.ned_of(g)` / `frame.geodetic_of_ned(ned)` |
| `ecef2enu` / `enu2ecef` | `frame.to_enu(e)` / `frame.to_ecef(enu)` |
| `geodetic2aer` / `aer2geodetic` | `frame.look_angles(g)` / `frame.destination(aer)` |
| `ecef2aer` / `aer2ecef` | `frame.look_angles_ecef(e)` / `frame.ecef_of_aer(aer)` |
| `enu2aer` / `aer2enu` | `enu.to_aer()` / `aer.to_enu()` |
| `ecef2enuv` / `enu2ecefv` (directions, velocities) | `frame.dir_to_enu(v)` / `frame.dir_to_ecef(enu)` |
| `lookAtSpheroid` | `frame.ground_strike(azimuth, elevation)` |
| `vdist` / `distance` | `WGS84.geodesic_inverse(a, b)` |
| `vreckon` / `reckon` | `WGS84.geodesic_direct(a, bearing, distance)` |
| `ecef2eci` / `eci2ecef` | `e.to_eci(gmst)` / `i.to_ecef(gmst)` (TEME; see above) |

Angles go in and come out in degrees; there is no `deg` flag to remember.

## Accuracy

| Calculation | Method | Accuracy |
|---|---|---|
| latitude/longitude/height → x, y, z | closed form | exact (to rounding) |
| x, y, z → latitude/longitude/height | Bowring, iterated | 1e-11° and 0.01 mm, from 11 km underground to geostationary orbit |
| local frames, look angles, ground strike, volumes | closed form | exact (to rounding), including straight up |
| surface distance and destination | Karney (as in GeographicLib) | about 15 nm, for every pair of points |
| great-circle helpers, horizon distance | sphere of mean radius | within 0.6 % |
| closest approach, tracks | constant velocity | exact for straight-line motion |
| datum shift (`Helmert7`) | small-angle Helmert | exact for the given parameters |
| inertial frame (`eci`) | sidereal time (TEME) | rotation exact; UTC-as-UT1 up to ~400 m |

The tests pin these to GeographicLib reference values, cross-check the conversions against an
independent crate ([`map_3d`](https://crates.io/crates/map_3d)) over a global grid, and cross-check
the geodesics against [`geographiclib-rs`](https://crates.io/crates/geographiclib-rs) on 25 000
pairs, 5 000 of them nearly antipodal, on four ellipsoids. They agree to within 4 nm.

## Optional features

All off by default, so the default build has no dependencies.

| Feature | Adds |
|---|---|
| `serde` | `Serialize` / `Deserialize` on every type |
| `mint` | `From` / `Into` for `mint` vectors and points (and through it `cgmath` and others) |
| `glam` | `From` / `Into` for `glam::DVec3` |
| `nalgebra` | `From` / `Into` for `nalgebra::Vector3` and `Point3` |
| `geo-types` | `Geodetic` ↔ `geo_types::Point` / `Coord` (x = longitude, y = latitude, height dropped) |

```bash
cargo add geo3d --features glam,serde
```

```rust,ignore
let v: glam::DVec3 = geo3d::Ecef::new(1.0, 2.0, 3.0).into();
let p: geo_types::Point = geo3d::Geodetic::new(-33.87, 151.21, 20.0).into(); // x = lon, y = lat
```

## Performance

Nothing allocates. `LocalFrame`, `Sphere`, `Cone` and `Cylinder` do their trigonometry once when
you build them, so build one per observer or shape and reuse it. `cargo bench` measures the hot
paths. Build with `RUSTFLAGS="-C target-cpu=native"` to use the CPU's fused multiply-add.

## Not included

- 2D geometry such as polygons and buffers: use [`geo`](https://crates.io/crates/geo).
- Map grids and projections such as UTM, MGRS or Web Mercator: use [`utm`](https://crates.io/crates/utm),
  [`proj`](https://crates.io/crates/proj) or [`geodesy`](https://crates.io/crates/geodesy).
- Height above sea level (geoid models), terrain, and atmospheric refraction.
- Precise Earth-orientation models (precession, nutation, polar motion).

## Glossary

- **Ellipsoid**: the smooth, slightly squashed sphere that positions are measured on. WGS84 is the
  one GPS uses.
- **Geodetic**: latitude, longitude and height above the ellipsoid.
- **ECEF** (Earth-Centred, Earth-Fixed): x, y, z in metres from the Earth's centre. The axes turn
  with the Earth.
- **ENU / NED** (East-North-Up / North-East-Down): offsets from an observer, on a flat plane that
  touches the Earth where they stand.
- **AER** (Azimuth, Elevation, Range): compass direction, angle above the horizon, and distance.
- **ECI** (Earth-Centred Inertial): like ECEF, but the axes don't turn. Used for satellites.
- **TEME**: the inertial frame SGP4 / TLE satellite propagators use; geo3d's `Eci`.
- **Geodesic**: the shortest path along the ellipsoid's surface.
- **Great circle**: the shortest path on a sphere. A quick approximation of a geodesic.
- **Datum**: a reference frame for coordinates. Helmert parameters convert between two of them.

## Develop

```bash
cargo test --all-features       # unit, cross-check and doc tests, including this README
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --check
cargo run --example tracking    # a ground station tracking aircraft
cargo bench                     # micro-benchmarks
```

CI runs these, plus a build on the minimum supported Rust version (1.85).

## License

Apache-2.0. See [LICENSE](LICENSE). Release notes are in [CHANGELOG.md](CHANGELOG.md).
