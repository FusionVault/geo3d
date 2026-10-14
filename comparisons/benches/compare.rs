//! Speed of geo3d and the other crates on the same inputs. `cargo bench` (from `comparisons/`).
//!
//! Every benchmark runs one library over the same 1 000 inputs, so each reported time is for
//! 1 000 calls: divide by 1 000 for the time per call. Inputs and results go through `black_box`,
//! and each adapter includes the unit and type conversions its users have to write (see
//! `src/lib.rs`).

use std::hint::black_box;
use std::time::Duration;

use criterion::measurement::WallTime;
use criterion::{criterion_group, criterion_main, BenchmarkGroup, Criterion, Throughput};
use geo3d_comparisons::*;

const N: usize = 1_000;

/// Time `f` over every input, once per iteration.
fn run<I: Copy, O>(g: &mut BenchmarkGroup<WallTime>, name: &str, inputs: &[I], f: impl Fn(I) -> Option<O>) {
    g.bench_function(name, |b| {
        b.iter(|| {
            for &x in inputs {
                black_box(f(black_box(x)));
            }
        })
    });
}

fn group<'a>(c: &'a mut Criterion, name: &str) -> BenchmarkGroup<'a, WallTime> {
    let mut g = c.benchmark_group(name);
    g.throughput(Throughput::Elements(N as u64));
    g
}

fn bench(c: &mut Criterion) {
    let mut rng = Rng::new(7);
    // Positions anywhere on Earth, from 500 m below the ellipsoid to airliner height.
    let points: Vec<Lla> = (0..N)
        .map(|_| {
            let [lat, lon] = rng.lat_lon();
            [lat, lon, rng.range(-500.0, 12_000.0)]
        })
        .collect();
    let ecefs: Vec<Xyz> = points.iter().map(|&p| with_geo3d::to_ecef(p).unwrap()).collect();
    // Satellites: 400 km to geostationary height.
    let orbits: Vec<Xyz> = (0..N)
        .map(|_| {
            let [lat, lon] = rng.lat_lon();
            with_geo3d::to_ecef([lat, lon, rng.range(400e3, 35_786e3)]).unwrap()
        })
        .collect();
    // One radar, and aircraft within a degree of it.
    let obs: Lla = [-33.94, 151.18, 21.0];
    let targets: Vec<Lla> = (0..N)
        .map(|_| {
            [
                obs[0] + rng.range(-1.0, 1.0),
                obs[1] + rng.range(-1.0, 1.0),
                rng.range(0.0, 12_000.0),
            ]
        })
        .collect();
    // Point pairs anywhere on Earth, and lines from a point on a bearing for up to 20 000 km.
    let pairs: Vec<(LatLon, LatLon)> = (0..N).map(|_| (rng.lat_lon(), rng.lat_lon())).collect();
    let lines: Vec<(LatLon, f64, f64)> = (0..N)
        .map(|_| (rng.lat_lon(), rng.range(0.0, 360.0), rng.range(1e3, 20_000e3)))
        .collect();

    let mut g = group(c, "geodetic_to_ecef");
    run(&mut g, "geo3d", &points, with_geo3d::to_ecef);
    run(&mut g, "map_3d", &points, with_map_3d::to_ecef);
    run(&mut g, "nav-types", &points, with_nav_types::to_ecef);
    run(&mut g, "sguaba", &points, with_sguaba::to_ecef);
    run(&mut g, "jord", &points, with_jord::to_ecef);
    g.finish();

    for (name, set) in [("ecef_to_geodetic", &ecefs), ("ecef_to_geodetic_orbit", &orbits)] {
        let mut g = group(c, name);
        run(&mut g, "geo3d", set, with_geo3d::to_geodetic);
        run(&mut g, "map_3d", set, with_map_3d::to_geodetic);
        run(&mut g, "nav-types", set, with_nav_types::to_geodetic);
        run(&mut g, "sguaba", set, with_sguaba::to_geodetic);
        run(&mut g, "jord", set, with_jord::to_geodetic);
        g.finish();
    }

    // One observer, many targets: each library's frame (or whatever it lets you keep) is built
    // once, outside the timed loop.
    let (g3, m3, nt, sg, jo) = (
        with_geo3d::Frame::new(obs).unwrap(),
        with_map_3d::Frame::new(obs).unwrap(),
        with_nav_types::Frame::new(obs).unwrap(),
        with_sguaba::Frame::new(obs).unwrap(),
        with_jord::Frame::new(obs).unwrap(),
    );
    let mut g = group(c, "enu_one_observer");
    run(&mut g, "geo3d", &targets, |t| g3.enu(t));
    run(&mut g, "map_3d", &targets, |t| m3.enu(t));
    run(&mut g, "nav-types", &targets, |t| nt.enu(t));
    run(&mut g, "sguaba", &targets, |t| sg.enu(t));
    run(&mut g, "jord", &targets, |t| jo.enu(t));
    g.finish();

    let mut g = group(c, "aer_one_observer");
    run(&mut g, "geo3d", &targets, |t| g3.aer(t));
    run(&mut g, "map_3d", &targets, |t| m3.aer(t));
    run(&mut g, "sguaba", &targets, |t| sg.aer(t));
    run(&mut g, "jord", &targets, |t| jo.aer(t));
    g.finish();

    // A single call: the frame is built on every call (the observer goes through `black_box`, so
    // it cannot be hoisted out of the loop).
    let mut g = group(c, "enu_single_call");
    run(&mut g, "geo3d", &targets, |t| {
        with_geo3d::Frame::new(black_box(obs))?.enu(t)
    });
    run(&mut g, "map_3d", &targets, |t| {
        with_map_3d::Frame::new(black_box(obs))?.enu(t)
    });
    run(&mut g, "nav-types", &targets, |t| {
        with_nav_types::Frame::new(black_box(obs))?.enu(t)
    });
    run(&mut g, "sguaba", &targets, |t| {
        with_sguaba::Frame::new(black_box(obs))?.enu(t)
    });
    run(&mut g, "jord", &targets, |t| {
        with_jord::Frame::new(black_box(obs))?.enu(t)
    });
    g.finish();

    let mut g = group(c, "aer_single_call");
    run(&mut g, "geo3d", &targets, |t| {
        with_geo3d::Frame::new(black_box(obs))?.aer(t)
    });
    run(&mut g, "map_3d", &targets, |t| {
        with_map_3d::Frame::new(black_box(obs))?.aer(t)
    });
    run(&mut g, "sguaba", &targets, |t| {
        with_sguaba::Frame::new(black_box(obs))?.aer(t)
    });
    run(&mut g, "jord", &targets, |t| {
        with_jord::Frame::new(black_box(obs))?.aer(t)
    });
    g.finish();

    // geo3d and geographiclib-rs return the distance and both azimuths; geo returns the distance only.
    let mut g = group(c, "geodesic_inverse");
    run(&mut g, "geo3d", &pairs, |(a, b)| with_geo3d::inverse(a, b));
    run(&mut g, "geographiclib-rs", &pairs, |(a, b)| {
        with_geographiclib::inverse(a, b)
    });
    run(&mut g, "geo", &pairs, |(a, b)| with_geo::geodesic_distance_m(a, b));
    g.finish();

    // geo3d and geographiclib-rs return the position and final azimuth; geo the position only.
    let mut g = group(c, "geodesic_direct");
    run(&mut g, "geo3d", &lines, |(a, az, s)| with_geo3d::direct(a, az, s));
    run(&mut g, "geographiclib-rs", &lines, |(a, az, s)| {
        with_geographiclib::direct(a, az, s)
    });
    run(&mut g, "geo", &lines, |(a, az, s)| {
        with_geo::geodesic_destination(a, az, s)
    });
    g.finish();

    let mut g = group(c, "great_circle_distance");
    run(&mut g, "geo3d", &pairs, |(a, b)| with_geo3d::great_circle_m(a, b));
    run(&mut g, "map_3d", &pairs, |(a, b)| with_map_3d::great_circle_m(a, b));
    run(&mut g, "nav-types", &pairs, |(a, b)| {
        with_nav_types::great_circle_m(a, b)
    });
    run(&mut g, "sguaba", &pairs, |(a, b)| with_sguaba::great_circle_m(a, b));
    run(&mut g, "jord", &pairs, |(a, b)| with_jord::great_circle_m(a, b));
    run(&mut g, "geo", &pairs, |(a, b)| with_geo::great_circle_m(a, b));
    g.finish();

    // Polygons of 10 vertices, 1 km to 100 km across, anywhere. Only the geodesic (ellipsoidal)
    // areas: the spherical ones are faster but off by about 0.7 % (see the accuracy run).
    let mut rng = Rng::new(11);
    let polygons: Vec<Vec<LatLon>> = (0..N)
        .map(|_| {
            let centre = [rng.range(-70.0, 70.0), rng.range(-175.0, 175.0)];
            let size = 10f64.powf(rng.range(-2.0, 0.0));
            star(&mut rng, centre, size, 10)
        })
        .collect();
    let polygons: Vec<&[LatLon]> = polygons.iter().map(Vec::as_slice).collect();
    let mut g = group(c, "polygon_area");
    run(&mut g, "geo3d", &polygons, with_geo3d::polygon_area_m2);
    run(
        &mut g,
        "geographiclib-rs",
        &polygons,
        with_geographiclib::polygon_area_m2,
    );
    run(&mut g, "geo", &polygons, with_geo::geodesic_area_m2);
    g.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(1))
        .sample_size(50);
    targets = bench
}
criterion_main!(benches);
