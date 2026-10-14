//! Accuracy and robustness of geo3d and the other crates, printed as Markdown tables.
//!
//! `cargo run --release --bin accuracy`. Deterministic: the same points every run.

use std::panic::{catch_unwind, AssertUnwindSafe};

use geo3d::{Aer, Geodetic, LocalFrame};
use geo3d_comparisons::*;

type Conv = fn([f64; 3]) -> Option<[f64; 3]>;
type Pair = fn(Lla, Lla) -> Option<[f64; 3]>;

/// (name, geodetic → ECEF, ECEF → geodetic)
const CONVERSIONS: [(&str, Conv, Conv); 5] = [
    ("geo3d", with_geo3d::to_ecef, with_geo3d::to_geodetic),
    ("map_3d", with_map_3d::to_ecef, with_map_3d::to_geodetic),
    ("nav-types", with_nav_types::to_ecef, with_nav_types::to_geodetic),
    ("sguaba", with_sguaba::to_ecef, with_sguaba::to_geodetic),
    ("jord", with_jord::to_ecef, with_jord::to_geodetic),
];

/// (name, ENU of a target, AER of a target), both including the frame construction.
const FRAMES: [(&str, Pair, Option<Pair>); 5] = [
    (
        "geo3d",
        |o, t| with_geo3d::Frame::new(o)?.enu(t),
        Some(|o, t| with_geo3d::Frame::new(o)?.aer(t)),
    ),
    (
        "map_3d",
        |o, t| with_map_3d::Frame::new(o)?.enu(t),
        Some(|o, t| with_map_3d::Frame::new(o)?.aer(t)),
    ),
    ("nav-types", |o, t| with_nav_types::Frame::new(o)?.enu(t), None),
    (
        "sguaba",
        |o, t| with_sguaba::Frame::new(o)?.enu(t),
        Some(|o, t| with_sguaba::Frame::new(o)?.aer(t)),
    ),
    (
        "jord",
        |o, t| with_jord::Frame::new(o)?.enu(t),
        Some(|o, t| with_jord::Frame::new(o)?.aer(t)),
    ),
];

const WGS84_B: f64 = 6_356_752.314_245_179;

fn main() {
    // Panics are caught and counted below; keep their messages out of the output.
    std::panic::set_hook(Box::new(|_| {}));
    round_trip();
    geodesics();
    areas();
    look_angles();
    robustness();
}

/// `None` if the call panicked.
fn guarded<T>(f: impl FnOnce() -> Option<T>) -> Result<Option<T>, ()> {
    catch_unwind(AssertUnwindSafe(f)).map_err(|_| ())
}

/// A result that is there, finite and not refused; anything else counts as a failure.
fn finite<const N: usize>(r: Result<Option<[f64; N]>, ()>) -> Option<[f64; N]> {
    r.ok().flatten().filter(|v| v.iter().all(|x| x.is_finite()))
}

/// The difference between two angles on the circle, degrees.
fn dang(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(360.0);
    d.min(360.0 - d)
}

/// A longitude into `[-180, 180)`.
fn wrap(lon: f64) -> f64 {
    (lon + 180.0).rem_euclid(360.0) - 180.0
}

fn sci(x: f64) -> String {
    if x == 0.0 {
        "0".into()
    } else {
        format!("{x:.1e}")
    }
}

fn table(title: &str, header: &[&str], rows: &[Vec<String>]) {
    println!("### {title}\n");
    println!("| {} |", header.join(" | "));
    println!("|{}", "---|".repeat(header.len()));
    for r in rows {
        println!("| {} |", r.join(" | "));
    }
    println!();
}

// ── 1. geodetic → ECEF → geodetic, each library with its own pair ────────────────────────────────

fn round_trip() {
    let heights = [-11_000.0, 0.0, 10_000.0, 400_000.0, 20_200_000.0, 35_786_000.0];
    // A 1° global grid, plus the last thousandth of a degree before each pole.
    let mut lats: Vec<f64> = (0..180).map(|i| -89.5 + i as f64).collect();
    lats.extend([89.9, 89.99, 89.999, -89.9, -89.99, -89.999]);
    let lons: Vec<f64> = (0..360).map(|i| -179.5 + i as f64).collect();

    let mut lat_rows = Vec::new();
    let mut lon_rows = Vec::new();
    let mut h_rows = Vec::new();
    for h in heights {
        let label = format!("{} km", h / 1000.0);
        let (mut lat_row, mut lon_row, mut h_row) = (vec![label.clone()], vec![label.clone()], vec![label]);
        for (_, to_ecef, to_geodetic) in CONVERSIONS {
            let (mut dlat, mut dlon, mut dh, mut failed) = (0.0f64, 0.0f64, 0.0f64, 0);
            for &lat in &lats {
                for &lon in &lons {
                    match finite(guarded(|| to_geodetic(to_ecef([lat, lon, h])?))) {
                        Some(r) => {
                            dlat = dlat.max((r[0] - lat).abs());
                            dlon = dlon.max(dang(r[1], lon));
                            dh = dh.max((r[2] - h).abs());
                        }
                        None => failed += 1,
                    }
                }
            }
            let note = if failed > 0 {
                format!(" ({failed} failed)")
            } else {
                String::new()
            };
            lat_row.push(sci(dlat) + &note);
            lon_row.push(sci(dlon) + &note);
            h_row.push(sci(dh) + &note);
        }
        lat_rows.push(lat_row);
        lon_rows.push(lon_row);
        h_rows.push(h_row);
    }
    let points = lats.len() * lons.len();
    println!("## Round trip: geodetic → ECEF → geodetic\n");
    println!("{points} points per height: a 1° grid plus latitudes ±89.9°, ±89.99°, ±89.999°. Each library uses its own forward and inverse conversion.\n");
    let header: Vec<&str> = ["Height"].into_iter().chain(CONVERSIONS.map(|c| c.0)).collect();
    table("Max latitude error (degrees)", &header, &lat_rows);
    table("Max longitude error (degrees)", &header, &lon_rows);
    table("Max height error (metres)", &header, &h_rows);
}

// ── 2. Geodesics against geographiclib-rs ───────────────────────────────────────────────────────

#[derive(Default)]
struct Stats {
    n: usize,
    failed: usize,
    max: f64,
    sum: f64,
    max_az: f64,
}

impl Stats {
    fn add(&mut self, err: f64, az_err: f64) {
        self.n += 1;
        self.max = self.max.max(err);
        self.sum += err;
        self.max_az = self.max_az.max(az_err);
    }

    fn row(&self, name: &str, set: &str, total: usize, with_az: bool) -> Vec<String> {
        vec![
            name.into(),
            set.into(),
            total.to_string(),
            self.failed.to_string(),
            sci(self.max),
            sci(self.sum / self.n.max(1) as f64),
            if with_az { sci(self.max_az) } else { "—".into() },
        ]
    }
}

fn geodesics() {
    let mut rng = Rng::new(0x9e37_79b9_7f4a_7c15);
    let random: Vec<(LatLon, LatLon)> = (0..100_000).map(|_| (rng.lat_lon(), rng.lat_lon())).collect();
    let antipodal: Vec<(LatLon, LatLon)> = (0..10_000)
        .map(|_| {
            let a = rng.lat_lon();
            let lat = (-a[0] + rng.range(-0.5, 0.5)).clamp(-90.0, 90.0);
            let lon = wrap(a[1] + 180.0 + rng.range(-0.5, 0.5));
            (a, [lat, lon])
        })
        .collect();

    let mut inverse_rows = Vec::new();
    let mut direct_rows = Vec::new();
    let mut gc: Vec<(&str, f64, Stats)> = vec![
        ("geo3d", 6_371_008.771, Stats::default()),
        ("map_3d", 6_371_000.0, Stats::default()),
        ("nav-types", 6_378_137.0, Stats::default()),
        ("sguaba", 6_378_137.0, Stats::default()),
        ("jord", 6_371_000.8, Stats::default()),
        ("geo", 6_371_008.8, Stats::default()),
    ];
    let gc_fns: [fn(LatLon, LatLon) -> Option<f64>; 6] = [
        with_geo3d::great_circle_m,
        with_map_3d::great_circle_m,
        with_nav_types::great_circle_m,
        with_sguaba::great_circle_m,
        with_jord::great_circle_m,
        with_geo::great_circle_m,
    ];

    for (set, pairs) in [("random", &random), ("nearly antipodal", &antipodal)] {
        let (mut g3_inv, mut geo_inv, mut g3_dir, mut geo_dir) =
            (Stats::default(), Stats::default(), Stats::default(), Stats::default());
        for &(a, b) in pairs.iter() {
            let reference = with_geographiclib::inverse(a, b).unwrap();
            let [s, azi1, azi2] = reference;

            match finite(guarded(|| with_geo3d::inverse(a, b))) {
                Some(r) => g3_inv.add((r[0] - s).abs(), dang(r[1], azi1).max(dang(r[2], azi2))),
                None => g3_inv.failed += 1,
            }
            match finite(guarded(|| with_geo::geodesic_distance_m(a, b).map(|d| [d]))) {
                Some([d]) => geo_inv.add((d - s).abs(), 0.0),
                None => geo_inv.failed += 1,
            }

            // Direct: travel the reference distance on the reference bearing; compare the arrival
            // with geographiclib's own direct solution, in metres along the surface.
            let [lat2, lon2, azi2_ref] = with_geographiclib::direct(a, azi1, s).unwrap();
            let miss = |p: LatLon| with_geographiclib::inverse([lat2, lon2], p).unwrap()[0];
            match finite(guarded(|| with_geo3d::direct(a, azi1, s))) {
                Some(r) => g3_dir.add(miss([r[0], r[1]]), dang(r[2], azi2_ref)),
                None => g3_dir.failed += 1,
            }
            match finite(guarded(|| with_geo::geodesic_destination(a, azi1, s))) {
                Some(p) => geo_dir.add(miss(p), 0.0),
                None => geo_dir.failed += 1,
            }

            for ((_, _, st), f) in gc.iter_mut().zip(gc_fns) {
                match finite(guarded(|| f(a, b).map(|d| [d]))) {
                    // The second slot (`max_az`) carries the absolute error here, in metres.
                    Some([d]) => st.add((d - s).abs() / s, (d - s).abs()),
                    None => st.failed += 1,
                }
            }
        }
        inverse_rows.push(g3_inv.row("geo3d", set, pairs.len(), true));
        inverse_rows.push(geo_inv.row("geo (`Geodesic`)", set, pairs.len(), false));
        direct_rows.push(g3_dir.row("geo3d", set, pairs.len(), true));
        direct_rows.push(geo_dir.row("geo (`Geodesic`)", set, pairs.len(), false));
    }

    println!("## Geodesics, against geographiclib-rs\n");
    println!("{} random point pairs (uniform over the sphere) and {} nearly antipodal pairs (the second point within 0.5° of the first one's antipode). A failure is `None`, NaN or a panic. geo's `Geodesic` calls geographiclib-rs, so it matches it exactly.\n", random.len(), antipodal.len());
    let header = [
        "Library",
        "Pairs",
        "Count",
        "Failures",
        "Max distance error (m)",
        "Mean distance error (m)",
        "Max azimuth error (°)",
    ];
    table(
        "Inverse: distance and azimuths between two points",
        &header,
        &inverse_rows,
    );
    let header = [
        "Library",
        "Lines",
        "Count",
        "Failures",
        "Max position error (m)",
        "Mean position error (m)",
        "Max final azimuth error (°)",
    ];
    table(
        "Direct: travel the reference distance on the reference bearing",
        &header,
        &direct_rows,
    );

    let rows: Vec<Vec<String>> = gc
        .iter()
        .map(|(name, r, st)| {
            vec![
                name.to_string(),
                format!("{r} m"),
                format!("{:.2} %", st.max * 100.0),
                format!("{:.2} %", st.sum / st.n.max(1) as f64 * 100.0),
                format!("{:.1} km", st.max_az / 1000.0),
                st.failed.to_string(),
            ]
        })
        .collect();
    table(
        "Great-circle (spherical) distance against the ellipsoidal geodesic, all pairs",
        &[
            "Library",
            "Sphere radius",
            "Max error",
            "Mean error",
            "Max error (km)",
            "Failures",
        ],
        &rows,
    );
}

// ── 3. Polygon areas against geographiclib-rs ──────────────────────────────────────────────────

type Area = fn(&[LatLon]) -> Option<f64>;

/// A percentage, readable at any size: scientific notation below 0.01 % and above 1 000 %.
fn pct(x: f64) -> String {
    if (1e-4..10.0).contains(&x) {
        format!("{:.2} %", x * 100.0)
    } else {
        format!("{} %", sci(x * 100.0))
    }
}

fn areas() {
    let mut rng = Rng::new(0x5851_f42d_4c95_7f2d);
    let mut sets: Vec<(&str, Vec<Vec<LatLon>>)> = Vec::new();
    // Up to 100 km across, so the edges' path barely matters: a geodesic, a great circle and a
    // straight line on an equal-area map (geo's spherical method) enclose almost the same area.
    let stars = |rng: &mut Rng, lo: f64, hi: f64, lon: Option<f64>| -> Vec<Vec<LatLon>> {
        (0..1_000)
            .map(|_| {
                let centre = [rng.range(-70.0, 70.0), lon.unwrap_or_else(|| rng.range(-175.0, 175.0))];
                let (size, n) = (10f64.powf(rng.range(lo, hi)), 3 + (rng.f64() * 30.0) as usize);
                star(rng, centre, size, n)
            })
            .collect()
    };
    sets.push(("fields, 10 m to 1 km across", stars(&mut rng, -4.0, -2.0, None)));
    sets.push(("regions, 1 km to 100 km across", stars(&mut rng, -2.0, -0.05, None)));
    sets.push((
        "across the antimeridian, 10 km to 100 km",
        stars(&mut rng, -1.05, -0.05, Some(180.0)),
    ));
    let rings = (0..1_000)
        .map(|_| {
            let (north, lat) = (rng.f64() < 0.5, rng.range(60.0, 89.9));
            let (n, start) = (3 + (rng.f64() * 30.0) as usize, rng.range(-180.0, 180.0));
            // Counter-clockwise seen from above: eastward round the north pole, westward round the south.
            let step = if north { 360.0 } else { -360.0 } / n as f64;
            (0..n)
                .map(|i| [if north { lat } else { -lat }, wrap(start + step * i as f64)])
                .collect()
        })
        .collect();
    sets.push(("around a pole, 60° to 89.9°", rings));

    let libs: [(&str, Area); 4] = [
        ("geo3d", with_geo3d::polygon_area_m2),
        ("geo (`GeodesicArea`)", with_geo::geodesic_area_m2),
        ("geo (`ChamberlainDuquetteArea`, sphere)", with_geo::spherical_area_m2),
        ("jord (spherical excess, sphere)", with_jord::polygon_area_m2),
    ];
    let mut rows = Vec::new();
    for (name, f) in libs {
        for (set, polys) in &sets {
            let (mut failed, mut max_abs, mut max_rel) = (0, 0.0f64, 0.0f64);
            for poly in polys {
                let reference = with_geographiclib::polygon_area_m2(poly).unwrap();
                match finite(guarded(|| f(poly).map(|a| [a]))) {
                    Some([a]) => {
                        let err = (a.abs() - reference.abs()).abs();
                        max_abs = max_abs.max(err);
                        max_rel = max_rel.max(err / reference.abs());
                    }
                    None => failed += 1,
                }
            }
            rows.push(vec![
                name.to_string(),
                set.to_string(),
                polys.len().to_string(),
                failed.to_string(),
                sci(max_abs),
                pct(max_rel),
            ]);
        }
    }
    println!("## Polygon areas, against geographiclib-rs\n");
    println!("Simple polygons of 3 to 32 vertices, counter-clockwise, with longitudes in [-180°, 180°) as GeoJSON has them. The rings round a pole have vertices on one parallel. Errors are in the size of the area, whatever its sign. A failure is `None`, NaN or a panic. geo's `GeodesicArea` calls geographiclib-rs, so it matches it exactly.\n");
    table(
        "Area of a polygon",
        &[
            "Library",
            "Polygons",
            "Count",
            "Failures",
            "Max error (m²)",
            "Max error (share of the area)",
        ],
        &rows,
    );
}

// ── 4. Local frames and look angles, compared pairwise ──────────────────────────────────────────

fn aer_to_enu(a: AerV) -> EnuV {
    let (sa, ca) = a[0].to_radians().sin_cos();
    let (se, ce) = a[1].to_radians().sin_cos();
    [a[2] * ce * sa, a[2] * ce * ca, a[2] * se]
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Micrometres, to three significant figures or so.
fn um(metres: f64) -> String {
    let x = metres * 1e6;
    match x {
        x if x < 10.0 => format!("{x:.3}"),
        x if x < 1_000.0 => format!("{x:.1}"),
        x => format!("{x:.0}"),
    }
}

/// Observer, target, and for constructed targets the pointing (azimuth, elevation, range) that
/// built them.
type Case = (Lla, Lla, Option<AerV>);

fn look_angles() {
    let mut rng = Rng::new(42);
    let (mut general, mut level, mut zenith): (Vec<Case>, Vec<Case>, Vec<Case>) = Default::default();
    for _ in 0..100 {
        let [lat, lon] = rng.lat_lon();
        let obs = [lat, lon, rng.range(0.0, 3_000.0)];
        for _ in 0..100 {
            // Nearby aircraft: within half a degree, up to 12 km high.
            let t = [
                (lat + rng.range(-0.5, 0.5)).clamp(-90.0, 90.0),
                wrap(lon + rng.range(-0.5, 0.5)),
                rng.range(0.0, 12_000.0),
            ];
            general.push((obs, t, None));
        }
        for i in 0..50 {
            // Anywhere on Earth, on the ground or in orbit.
            let [tl, tn] = rng.lat_lon();
            general.push((obs, [tl, tn, [0.0, 400e3, 20_200e3, 35_786e3][i % 4]], None));
        }
        // On and just off due north, east, south and west, and on and just off straight up: where
        // formulas built on asin or acos lose precision. Built from a known pointing, which is
        // then the reference.
        let frame = LocalFrame::new(Geodetic::new(obs[0], obs[1], obs[2]));
        let build = |az: f64, el: f64, range: f64| {
            let p = frame.destination(Aer::new(az, el, range));
            (obs, [p.lat_deg, p.lon_deg, p.height_m], Some([az, el, range]))
        };
        for az in [0.0, 90.0, 180.0, 270.0] {
            for off in [0.0, 1e-7, 1e-5] {
                for range in [10e3, 1_000e3] {
                    level.push(build(az + off, 0.0, range));
                }
            }
        }
        for off in [0.0, 1e-7, 1e-6, 1e-5, 1e-3] {
            for range in [1e3, 400e3, 20_200e3, 35_786e3] {
                zenith.push(build(rng.range(0.0, 360.0), 90.0 - off, range));
            }
        }
    }

    println!("## Local frames and look angles\n");
    println!("Every library is compared with geo3d and with map_3d. For targets built from a known pointing, that pointing is the reference too. Disagreements are in micrometres. For look angles, each answer (azimuth, elevation, range) is turned back into a point and the distance between the two points is reported, so a pointing error counts by how far it moves the target.\n");
    let titles = [
        format!(
            "{} observer–target pairs: nearby aircraft, and targets anywhere up to geostationary height",
            general.len()
        ),
        format!(
            "{} pairs: on the horizon, on and just off (1e-7°, 1e-5°) due north, east, south and west, at 10 km and 1 000 km",
            level.len()
        ),
        format!(
            "{} pairs: straight up and just off it (1e-7° to 1e-3°), 1 km to geostationary range",
            zenith.len()
        ),
    ];
    for (title, set) in titles.iter().zip([&general, &level, &zenith]) {
        let worst = |f: &dyn Fn(&Case) -> Option<f64>| {
            um(set.iter().map(|c| f(c).unwrap_or(f64::INFINITY)).fold(0.0, f64::max))
        };
        let mut rows = Vec::new();
        for (name, enu, aer) in FRAMES {
            let mut row = vec![name.to_string()];
            for (_, ref_enu, _) in [FRAMES[0], FRAMES[1]] {
                row.push(worst(&|&(o, t, _)| {
                    Some(dist(finite(guarded(|| enu(o, t)))?, ref_enu(o, t)?))
                }));
            }
            let Some(aer) = aer else {
                row.extend(["—".to_string(), "—".into(), "—".into(), "—".into()]);
                rows.push(row);
                continue;
            };
            let ours = |o, t| finite(guarded(|| aer(o, t))).map(aer_to_enu);
            for (_, _, ref_aer) in [FRAMES[0], FRAMES[1]] {
                let ref_aer = ref_aer.unwrap();
                row.push(worst(&|&(o, t, _)| Some(dist(ours(o, t)?, aer_to_enu(ref_aer(o, t)?)))));
            }
            if set[0].2.is_some() {
                row.push(worst(&|&(o, t, built)| Some(dist(ours(o, t)?, aer_to_enu(built?)))));
                // The angle between the two pointings: the chord between unit directions.
                let unit = |a: AerV| aer_to_enu([a[0], a[1], 1.0]);
                let angle = set
                    .iter()
                    .map(|&(o, t, built)| match finite(guarded(|| aer(o, t))) {
                        Some(a) => dist(unit(a), unit(built.unwrap())).to_degrees(),
                        None => f64::INFINITY,
                    })
                    .fold(0.0, f64::max);
                row.push(sci(angle));
            } else {
                row.extend(["—".to_string(), "—".into()]);
            }
            rows.push(row);
        }
        table(
            title,
            &[
                "Library",
                "ENU vs geo3d (µm)",
                "ENU vs map_3d (µm)",
                "Look angles vs geo3d (µm)",
                "Look angles vs map_3d (µm)",
                "Look angles vs the built pointing (µm)",
                "Pointing vs the built pointing (°)",
            ],
            &rows,
        );
    }
}

// ── 5. Robustness: edge cases ───────────────────────────────────────────────────────────────────

/// What came back, in a word: ok, NaN, rejected, panic, or wrong (with the value).
fn verdict<const N: usize>(
    r: Result<Option<[f64; N]>, ()>,
    ok: impl Fn(&[f64; N]) -> bool,
    show: impl Fn(&[f64; N]) -> String,
) -> String {
    match r {
        Err(()) => "panic".into(),
        Ok(None) => "rejected (`None`)".into(),
        Ok(Some(v)) if v.iter().all(|x| x.is_nan()) => "NaN".into(),
        Ok(Some(v)) if v.iter().any(|x| x.is_nan()) => format!("partly NaN ({})", show(&v)),
        Ok(Some(v)) if ok(&v) => "ok".into(),
        Ok(Some(v)) => format!("wrong ({})", show(&v)),
    }
}

fn lla(v: &[f64; 3]) -> String {
    format!("{:.6}°, {:.6}°, {:.4} m", v[0], v[1], v[2])
}

fn robustness() {
    let close_lla = |want: Lla| {
        move |v: &Lla| (v[0] - want[0]).abs() < 1e-8 && dang(v[1], want[1]) < 1e-8 && (v[2] - want[2]).abs() < 1e-3
    };
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut row = |case: &str, cells: Vec<String>| rows.push([vec![case.to_string()], cells].concat());

    // Conversions and frames: the first five columns; the geodesic-only crates get "—".
    let dash = || vec!["—".to_string(); 2];

    let centre: Vec<String> = CONVERSIONS
        .iter()
        .map(|(_, _, to_geo)| {
            verdict(
                guarded(|| to_geo([0.0, 0.0, 0.0])),
                |v| (v[0].abs() - 90.0).abs() < 1e-8 && (v[2] + WGS84_B).abs() < 1e-3,
                lla,
            )
        })
        .collect();
    row(
        "ECEF → geodetic at the Earth's centre (0, 0, 0); want ±90°, h = −b",
        [centre, dash()].concat(),
    );

    let pole: Vec<String> = CONVERSIONS
        .iter()
        .map(|(_, _, to_geo)| {
            verdict(
                guarded(|| to_geo([0.0, 0.0, WGS84_B + 100.0])),
                |v| (v[0] - 90.0).abs() < 1e-8 && (v[2] - 100.0).abs() < 1e-3,
                lla,
            )
        })
        .collect();
    row(
        "ECEF → geodetic exactly above the North Pole (0, 0, b + 100 m)",
        [pole, dash()].concat(),
    );

    let south: Vec<String> = CONVERSIONS
        .iter()
        .map(|(_, to_ecef, to_geo)| {
            verdict(
                guarded(|| to_geo(to_ecef([-90.0, 0.0, 100.0])?)),
                |v| (v[0] + 90.0).abs() < 1e-8 && (v[2] - 100.0).abs() < 1e-3,
                lla,
            )
        })
        .collect();
    row(
        "Round trip exactly at the South Pole (−90°, 0°, 100 m)",
        [south, dash()].concat(),
    );

    // 0.0002° of longitude along the equator, across ±180°.
    let dl = 2e-4f64.to_radians();
    let want = [6_378_137.0 * dl.sin(), 0.0, 6_378_137.0 * (dl.cos() - 1.0)];
    let anti: Vec<String> = FRAMES
        .iter()
        .map(|(_, enu, _)| {
            verdict(
                guarded(|| enu([0.0, 179.9999, 0.0], [0.0, -179.9999, 0.0])),
                |v| dist(*v, want) < 1e-3,
                |v| format!("E {:.3} m, N {:.3} m, U {:.3} m", v[0], v[1], v[2]),
            )
        })
        .collect();
    row(
        "ENU across the antimeridian, 179.9999° → −179.9999° (want 22.264 m east)",
        [anti, dash()].concat(),
    );

    let far: Vec<String> = CONVERSIONS
        .iter()
        .map(|(_, to_ecef, to_geo)| {
            verdict(
                guarded(|| to_geo(to_ecef([45.0, 45.0, 1e8])?)),
                close_lla([45.0, 45.0, 1e8]),
                lla,
            )
        })
        .collect();
    row(
        "Round trip at 1e8 m altitude (45°, 45°, 100 000 km)",
        [far, dash()].concat(),
    );

    let nan_fwd: Vec<String> = CONVERSIONS
        .iter()
        .map(|(_, to_ecef, _)| verdict(guarded(|| to_ecef([f64::NAN, 0.0, 0.0])), |_| false, lla))
        .collect();
    row("Geodetic → ECEF with a NaN latitude", [nan_fwd, dash()].concat());

    let nan_inv: Vec<String> = CONVERSIONS
        .iter()
        .map(|(_, _, to_geo)| verdict(guarded(|| to_geo([f64::NAN, 0.0, 0.0])), |_| false, lla))
        .collect();
    row("ECEF → geodetic with a NaN x", [nan_inv, dash()].concat());

    // Geodesics: geo3d, geographiclib-rs and geo only.
    let dist_m = |v: &[f64; 3]| format!("{:.4} m", v[0]);
    let (a, b) = ([0.0, 0.0], [0.0, 180.0]);
    let half_meridian = 20_003_931.458_6;
    let ok = |v: &[f64; 3]| (v[0] - half_meridian).abs() < 1e-3;
    let cells = vec![
        verdict(guarded(|| with_geo3d::inverse(a, b)), ok, dist_m),
        "—".into(),
        "—".into(),
        "—".into(),
        "—".into(),
        verdict(guarded(|| with_geographiclib::inverse(a, b)), ok, dist_m),
        verdict(
            guarded(|| with_geo::geodesic_distance_m(a, b).map(|d| [d, 0.0, 0.0])),
            ok,
            dist_m,
        ),
    ];
    row(
        "Geodesic between exact antipodes (0°, 0°) → (0°, 180°); want 20 003 931.459 m",
        cells,
    );

    let (a, b) = ([f64::NAN, 0.0], [10.0, 10.0]);
    let cells = vec![
        verdict(guarded(|| with_geo3d::inverse(a, b)), |_| false, dist_m),
        "—".into(),
        "—".into(),
        "—".into(),
        "—".into(),
        verdict(guarded(|| with_geographiclib::inverse(a, b)), |_| false, dist_m),
        verdict(
            guarded(|| with_geo::geodesic_distance_m(a, b).map(|d| [d])),
            |_| false,
            |v| format!("{:.4} m", v[0]),
        ),
    ];
    row("Geodesic with a NaN latitude", cells);

    println!("## Robustness: edge cases\n");
    println!("\"ok\" means within 1 mm and 1e-8°. For NaN input, NaN out or a rejection (`None`) are both fine; a panic is not.\n");
    table(
        "What each library returns",
        &[
            "Case",
            "geo3d",
            "map_3d",
            "nav-types",
            "sguaba",
            "jord",
            "geographiclib-rs",
            "geo",
        ],
        &rows,
    );
}
