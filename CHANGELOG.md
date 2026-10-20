# Changelog

All notable changes to geo3d. We use [semantic versioning](https://semver.org), and each
release is tagged `vX.Y.Z`.

## Unreleased

- 🗺️ **A point-in-polygon test for a polygon used once** ([#68](https://github.com/FusionVault/geo3d/issues/68)): `polygon_containment` and `multi_polygon_containment`, behind the `footprint` feature.
  - It tests edge by edge with no index, and allocates nothing. A ring is any iterator of positions, so a polygon can be tested straight from where it lies, such as one row of a column of WKB.
  - A small polygon takes about 13 ns a test, level with `geo::Contains` on a polygon already built, and 16 times faster than preparing a `Footprint` for the one test.
  - It gives a `Footprint`'s answer for the same rings, valid or not: the two share one edge rule and one rule for combining rings. `Footprint` is as fast as before; it compiles to the same instructions.

## 0.3.2 — 2026-10-20

- 🎯 **Exact heights on the polar axis and in the equatorial plane** ([#63](https://github.com/FusionVault/geo3d/issues/63)), as `to_geodetic` promises. There the height is `|z| − b` or `r − a`, rounded once. Before, it came back up to a nanometre off: `(0°, 0°, 100 m)` round-tripped to `99.99999999969262` m. Every other position gives the same answer as 0.3.1, to the bit (checked on a million).
- ⚡ **`Footprint::containment` is inlined** ([#65](https://github.com/FusionVault/geo3d/issues/65)), as geo3d's other hot functions are, so a caller's loop tests points without a call into the crate. That makes it 5% to 6% faster per point in a caller's loop, measured.

## 0.3.1 — 2026-10-20

- 🗺️ **Fast point-in-polygon: `Footprint` and `Prism`** ([#60](https://github.com/FusionVault/geo3d/issues/60)), behind the new `footprint` feature:
  - **`Footprint`** prepares a polygon once, with holes and several parts. Then `containment` places a point as inside, on the boundary or outside.
  - **Each ring's edges are indexed by latitude,** so a point meets only the edges that cross its latitude, after a bounding-box reject. A point takes about 25 ns on a 50-vertex ring, against 55 ns for `geo::Contains`. On a 5,000-vertex coastline it takes 57 ns, against 3.1 µs.
  - **The answer is exact,** with Shewchuk's orientation test (the `robust` crate, now an optional dependency), and matches `geo`'s `coordinate_position` point for point.
  - **`Prism`** extrudes a footprint through an `AltitudeBand`, as a `Volume`.
  - With `geo-types`, a `Polygon` or `MultiPolygon` converts into a `Footprint`.
  - The feature needs an allocator, not the standard library.
- 🧭 **`Geodetic::from_lon_lat`** builds a position longitude first, the order of GeoJSON, WKT and GeoArrow, so the two cannot be swapped at those boundaries.
- 🔍 **Accuracy pass.** A review of every formula for precision traps, with each finding checked
  against exact answers. Fixed:
  - `line_of_sight_with_refraction` bent the whole path, even far above the air, so it called a
    geostationary satellite 18° below the horizon visible. The air now bends only the part of the
    path below 7.35 km (ITU-R P.453's refractivity scale height). Links through the air are
    unchanged, and a distant satellite shows down to 0.8° below the horizon (ITU-R P.834: 0.76°).
    `horizon_distance_with_refraction_m` uses the same model above 7.35 km;
  - `Cylinder` measured its radius flat in the centre's frame, which shrank it with height: 107 m
    short at 13.7 km up on a 50 km zone. The radius is now along the ground (the geodesic), at any
    height;
  - `geodesic_inverse` gave a distance of 0 m, not `None`, for a NaN or out-of-range latitude from
    some starting points (the equator at 0°, for one), and `geodesic_direct` gave a finite latitude
    with a NaN longitude. Both now check their input;
  - `initial_bearing_deg` lost precision on short paths: 8e-5° out at 1 mm, and 0.02° at 1 µm. It
    now keeps it at any length;
  - `ray_intersect`, `ground_strike` and `line_of_sight` work their terms to 32 digits, so rays that
    graze the Earth from far away, and points just above the surface, keep full precision;
  - `Cone` compares angles, not cosines, so beams narrower than about a microradian keep their
    shape;
  - `Geodetic::normalized` left a longitude already in range very slightly changed (−0.1° became
    −0.10000000000002274°); it's now exact;
  - `Ellipsoid::ep2` is 300 times more precise.

  The comparison results are unchanged.
- 🌐 **WebAssembly** ([#19](https://github.com/FusionVault/geo3d/issues/19)). CI builds for the
  browser (`wasm32-unknown-unknown`, with and without std, and with every feature), and runs the
  whole test suite under WebAssembly (`wasm32-wasip1`, in wasmtime). A geoid test that assumed a
  64-bit `usize` now holds on 32-bit too.
- 📚 **docs.rs.** Items that need an optional feature now say so ("Available on crate feature
  serde only"), and a nightly CI job builds the docs as docs.rs does. The examples show on the
  items they use again: the rewrite had stopped docs.rs from scraping them.
- 🌐 **Every EPSG ellipsoid** ([#18](https://github.com/FusionVault/geo3d/issues/18)):
  - `Ellipsoid::from_epsg(code)` gives all 50 ellipsoids in the EPSG registry, each from the
    numbers EPSG defines it by, converted to metres where EPSG gives feet;
  - new named ones: `CGCS2000` (China, BeiDou), `KRASSOVSKY_1940` (Pulkovo 1942, Beijing 1954),
    `GSK_2011` (Russia), `IAG_1975` (Xian 1980), `GRS67`, `AUSTRALIAN_NATIONAL` (AGD66 and AGD84),
    `WGS66`, `WGS60`, `EVEREST_1830` and `EVEREST_1830_MODIFIED`. Every named ellipsoid now gives
    its EPSG code and where it's used;
  - `Ellipsoid::from_axes(a, b)`, for ellipsoids defined by their two radii.
- 🐛 `CLARKE_1866` now uses its defining radii. Its polar radius was 1 µm short.
- 🚀 **Examples, rewritten.** Eight short programs, one real job each: `radar`, `air_traffic`,
  `satellite`, `drone_camera`, `field_survey`, `search_and_rescue`, `flight_route` and `sea_level`.
  They replace `tracking`, and CI runs them all. The README lists them.
- 🎯 **Range rate and intercepts** ([#17](https://github.com/FusionVault/geo3d/issues/17)):
  - `Track::range_rate_mps` and `LocalFrame::range_rate_mps`: how fast the distance is changing,
    negative while closing, for radar Doppler. The docs give the Doppler shift for radar and
    one-way signals;
  - `Track::intercept_from` and `Track::slowest_intercept_from`: when and where a pursuer flying
    straight at a steady speed meets a target, which way to fly, and the slowest speed that works,
    as an `Intercept`. The equation's terms are worked out to 32 digits, so the answer is within a
    few nanometres even at equal speeds, where the textbook formula divides by zero (it fails in
    half of those cases). jord's sphere solver fails in 1 % to 16 % of its cases;
  - `Track::from_eci`: an Earth-fixed track from an inertial (TEME) position and velocity, as SGP4
    gives them. It takes the Earth's spin off at the exact rate of geo3d's sidereal time.
    `examples/satellite.rs` now prints the Doppler shift on the ISS's radio.
- 🌊 **Heights above sea level** ([#16](https://github.com/FusionVault/geo3d/issues/16)):
  - the `Geoid` trait: `undulation_m`, the geoid's height above the ellipsoid, then
    `sea_level_height_m` and `geodetic_from_sea_level` to convert. Any `Fn(lat, lon) -> metres` is
    a `Geoid`, so national models plug straight in;
  - `GeoidGrid`, which reads GeographicLib's EGM84, EGM96 and EGM2008 grids (`.pgm`) straight from
    their bytes: from a file, from `include_bytes!`, or memory-mapped. It never copies or allocates,
    and works without the standard library.

  It interpolates as GeographicLib does (a cubic on 12 points by default, or bilinear), and gives
  the same answers as GeographicLib's GeoidEval. Against the gravity models themselves, at
  GeographicLib's 500 000 test points, it is within 2.6 mm of EGM96 and 2.2 mm of EGM2008. The
  egm96 crate's default is up to 1.9 m out, and egm2008's up to 24 m.
- 🖨️ **Printing** ([#15](https://github.com/FusionVault/geo3d/issues/15)). `Display` for `Geodetic`,
  `Aer`, `Ecef`, `Eci`, `Enu`, `Ned`, `Frd`, `Vec3`, `Attitude`, `Quaternion` and `UtmUps`:
  - plain `{}` keeps every digit, so a printed position reads back exactly;
  - `{:.N}` gives N decimals (metres get four fewer where a type mixes them with degrees);
  - `{:#}` writes a position's hemispheres as letters (`33.87° S, 151.21° E, 20 m`).

  It works without the standard library too.
- 🗺️ **UTM and UPS** ([#14](https://github.com/FusionVault/geo3d/issues/14)).
  `Ellipsoid::to_utm_ups`, `to_utm_ups_zone` (a chosen zone) and `from_utm_ups`. Standard zones
  include the Norway and Svalbard exceptions, and UPS covers the poles. Built on two new projections
  for any ellipsoid:
  - `TransverseMercator`, by Krüger's series to the eighth order (Karney 2011), with grid
    convergence and scale;
  - `PolarStereographic`, exact.

  Within 2 nm of GeographicLib on its own reference points, and 4 nm grid-and-back anywhere. The
  `utm` crate is up to 77 cm out on a round trip.
- 🎯 Sine and cosine of whole multiples of 30° and 45° are now correctly rounded, as in current
  GeographicLib. Results at those exact angles can change in the last digit.
- 🔄 **Quaternions.** `Quaternion { w, x, y, z }` in the Hamilton convention, scalar first, turning
  body vectors into NED, as PX4 and MAVLink store attitude. `Attitude::to_quaternion` and
  `from_quaternion` (exact pointing straight up), `Quaternion::to_rotation` and `from_rotation`
  (Shepperd's method, precise for every rotation, half turns included), `*` to combine turns, and
  `conjugate` to reverse one. Always the same sign for the same turn. Converts to `glam::DQuat`,
  `nalgebra::Quaternion` and `mint::Quaternion`; glam's and nalgebra's own yaw, pitch and roll agree
  with it.
- ✈️ **Yaw, pitch and roll** ([#13](https://github.com/FusionVault/geo3d/issues/13)). `Attitude`
  (aerospace yaw, pitch and roll, from North-East-Down) and `Frd` (forward, right, down in a
  vehicle's own frame):
  - `Attitude::body_to_ned`, `ned_to_body`, `to_rotation` and `from_rotation`;
  - `LocalFrame::geodetic_of_body` and `body_of`: positions around a vehicle;
  - `LocalFrame::body_to_ecef_rotation` and `attitude_from_ecef_rotation`, the rotation into
    ECEF and back (Cesium's `headingPitchRollToFixedFrame`);
  - `LocalFrame::enu_to_ecef_rotation` and `ned_to_ecef_rotation` expose the frame's own rotation.

  Exact pointing straight up or down: the angles always rebuild the same orientation, where
  sguaba (and jord, after a trip through ECEF) can come back pointing the wrong way.
- 🔌 **Runs without the standard library** ([#12](https://github.com/FusionVault/geo3d/issues/12)).
  On microcontrollers and flight controllers, `default-features = false, features = ["libm"]` takes
  the maths from `libm`. Everything else is the same, every optional feature works too, and CI
  builds it for an Arm Cortex-M4F. Nothing changes for anyone else: the new `std` feature is on by
  default, and the default build still has no dependencies.
  - ⚠️ If you already set `default-features = false`, add `features = ["std"]` (or `libm`), or
    the build will stop with a message saying so.
  - Keep `std` where you have it: its sine, cosine and atan2 round the last digit more carefully
    than libm's (worst 0.5 against 0.7–1.2 units in the last place). libm is 10–20 % faster.
- 🧭 **Route tools** ([#7](https://github.com/FusionVault/geo3d/issues/7)), built on n-vectors so
  the poles and the antimeridian need no special care:
  - `cross_track_distance_m` and `along_track_distance_m`: how far off a route, and how far along;
  - `great_circle_intersection`: where two routes cross;
  - `mean_position`: the centre of a group of positions;
  - `great_circle_interpolate`: a point part-way along a route, at constant speed.

  They stay exact on short paths. The textbook n-vector formula can be kilometres out on a path a
  few micrometres long; geo3d builds the path from the latitude and longitude differences instead.
- 🐛 **`great_circle_destination` near the poles.** It computed the latitude with `asin`, which
  loses precision near ±90°: up to 11 cm out arriving close to a pole. It now uses n-vectors and
  `atan2`, and is exact everywhere.
- 📚 **The NavLab examples** (`geo3d::navlab`): all ten of FFI's classic n-vector problems, solved
  with geo3d and checked against the published answers.
- 🎯 **Uncertainty between frames** ([#6](https://github.com/FusionVault/geo3d/issues/6)). `Cov3`
  is a 3×3 covariance. It moves between frames with the position:
  - ENU ↔ ECEF (`LocalFrame::enu_cov_to_ecef`, `ecef_cov_to_enu`) and ENU ↔ NED, exactly;
  - AER ↔ ENU (`Aer::cov_to_enu`, `Enu::cov_to_aer`) and geodetic ↔ ECEF
    (`Ellipsoid::geodetic_cov_to_ecef`, `ecef_cov_to_geodetic`), to first order.

  The Jacobians are checked against finite differences, and the results against a Monte Carlo
  run. `Cov3` converts to `nalgebra::Matrix3`, `glam::DMat3` and the `mint` matrices.
- 🗺️ **Polygon area** ([#4](https://github.com/FusionVault/geo3d/issues/4)).
  `Ellipsoid::polygon_area` gives the area and perimeter of a polygon whose edges are geodesics. It
  is a port of GeographicLib's method, and matches it to within 0.2 m² on any polygon, including
  ones around a pole or across the antimeridian. Counter-clockwise polygons have a positive area.
- 🛰️ **Satellite time** ([#8](https://github.com/FusionVault/geo3d/issues/8)):
  - `eci::gmst_from_unix_with_dut1` corrects UTC to UT1, which removes up to 400 m of error.
  - `eci::era_rad` and `eci::era_from_unix_with_dut1` give the Earth Rotation Angle (IAU 2000).
  - Sidereal time from Unix time now resolves 0.2 µs instead of 40 µs.
  - Both formulas match SOFA's reference values to 1e-12 rad.
  - A new example, `examples/satellite.rs`, tracks the ISS with the `sgp4` crate.
- 📡 **Refraction for radar and radio** ([#9](https://github.com/FusionVault/geo3d/issues/9)).
  `horizon_distance_with_refraction_m` and `line_of_sight_with_refraction` take an
  effective-Earth-radius factor `k`. `STANDARD_REFRACTION_K` is 4/3, the ITU-R standard for radio,
  which reaches about 15 % further than a straight line.
- 🐛 **`line_of_sight` between points exactly on the surface.** It used to say two points on the
  ground could see each other at any distance, even when the line between them went hundreds of
  metres underground. Now it checks the deepest point of the line, so they only see each other if
  the line stays above the surface.
- 🛡️ **An edge-case test suite** ([#10](https://github.com/FusionVault/geo3d/issues/10)). It runs
  every public conversion at the exact poles, the antimeridian, the Earth's centre, 100 000 km up,
  deep underground, exact antipodes, and NaN or infinite input. It found and fixed two bugs:
  - Bearings and azimuths could come out as exactly 360° (from a tiny negative angle). They are
    now always in `[0, 360)`. This affects look angles, `Course`, `initial_bearing_deg` and
    geodesic bearings.
  - `ground_strike` returned `Some(NaN position)` for a NaN latitude or longitude at height 0.
    It now returns `None`.
- 🧹 CI now runs [cargo-semver-checks](https://github.com/obi1kenobi/cargo-semver-checks) on every
  pull request, so accidental breaking changes are caught before release
  ([#11](https://github.com/FusionVault/geo3d/issues/11)).
- 🎯 **Faster, more accurate x, y, z → latitude/longitude/height**
  ([#28](https://github.com/FusionVault/geo3d/issues/28)). `to_geodetic` now uses an exact formula
  (Vermeille's, refined by Karney for GeographicLib) and one correction step. It no longer
  iterates.
  - It is about 40 % faster.
  - It is at least as accurate at every height: within 3 nm near the Earth and 20 nm at satellite
    heights. This was measured at 50-digit precision on 5 000 points.
  - Deep inside the Earth, the old method could be up to 53 m out. The new one stays within 3 nm.
  - A point exactly on the polar axis now always reports longitude 0.

## 0.3.0 — 2026-10-14

An accuracy release. Your code keeps compiling: no public type, method or signature has changed.
Some results are now more accurate, as listed below.

### 🎯 More accurate

- **Distance along the ground now works for every pair of points**
  ([#3](https://github.com/FusionVault/geo3d/issues/3)). `geodesic_inverse` and `geodesic_direct`
  now use Karney's algorithm, the one GeographicLib uses, instead of Vincenty's.
  - It is accurate to about 15 nm, up from about 0.5 mm.
  - `geodesic_inverse` used to return `None` for points on nearly opposite sides of the Earth. Now
    it only returns `None` for NaN or infinite input, or a latitude beyond ±90°.
  - It matches geographiclib-rs to within 4 nm on 25 000 test pairs.
  - The code is ported from GeographicLib and geographiclib-rs, which are both MIT-licensed. See
    NOTICE.
- **Four precision fixes** ([#5](https://github.com/FusionVault/geo3d/issues/5)). Each has a test
  that failed before the fix:
  - **Elevation near straight up.** It used to round to exactly 90°, up to 0.6 m out at
    geostationary range. This affects every look angle.
  - **Great-circle distance.** It was up to 0.9 m short near the far side of the Earth. It is now
    accurate to about a nanometre everywhere.
  - **`Vec3::angle_to`.** Tiny angles used to come out as 0.
  - **Ray–ground intersections.** These lost precision when starting close to the surface. This
    affects `ray_intersect`, `ground_strike` and `line_of_sight`.

### 📖 Clearer docs

- The README and the docs.rs front page were rewritten in plain language. They now have recipes, a
  "find what you need" table, an accuracy table, a glossary, and a lookup table for MATLAB, pymap3d
  and map_3d users.
- There is a new [comparison](comparisons/README.md) against five other Rust crates, with charts.
  It covers accuracy, tricky cases, features and speed, including the places where geo3d is
  slower.
- There is a new [contributing guide](CONTRIBUTING.md).
- A correction: the `eci` frame is TEME, the one SGP4 satellite models use. It is not J2000, which
  differs by about 0.4°. Treating UTC as UT1 can add up to 400 m. The old docs promised "tens of
  metres", which wasn't true.
- docs.rs now shows the optional features. Searching for familiar names like `geodetic2aer`,
  `vdist` or `haversine` now finds the matching geo3d function.

### ⚡ Speed

- `geodesic_direct` is a little faster than before, and about 13 % faster than geographiclib-rs.
- `geodesic_inverse` takes about twice as long as before (about 0.86 µs per call). That extra work
  is what makes it exact everywhere. It is within 4 % of geographiclib-rs.
- Look angles and great-circle distance are 15–20 ns slower per call, because of the precision
  fixes.

### 🧹 Housekeeping

- The code is split into one module per topic.
- The crates.io description is simpler. The `no-std::no-alloc` category was removed, because geo3d
  needs `std` today ([#12](https://github.com/FusionVault/geo3d/issues/12) will change that).
- CI now checks that every code change comes with a CHANGELOG line.

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
