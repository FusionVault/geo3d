//! Micro-benchmarks for the hot paths. Run with `cargo bench`; for the FMA-enabled numbers add
//! `RUSTFLAGS="-C target-cpu=native"`.
//!
//! These exist to make speed measurable and to catch regressions — not to chase a headline number.
//! The per-call cost of most of these is dominated by hardware `sin`/`cos`/`atan2`/`sqrt`.

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};
use geo3d::{Geodetic, Helmert7, LocalFrame, Vec3, WGS84};

fn bench(c: &mut Criterion) {
    let g = Geodetic::new(-33.8688, 151.2093, 100.0);
    let e = WGS84.to_ecef(g);
    let frame = LocalFrame::new(g);
    let target = Geodetic::new(-33.80, 151.30, 2_500.0);
    let (a, b) = (Vec3::new(1.0, 2.0, 3.0), Vec3::new(-4.0, 5.0, -6.0));
    let helmert = Helmert7::position_vector(-446.448, 125.157, -542.060, -0.1502, -0.2470, -0.8421, 20.4894);
    let lax = Geodetic::new(33.9425, -118.4081, 0.0);
    let jfk = Geodetic::new(40.6413, -73.7781, 0.0);

    c.bench_function("geodetic_to_ecef", |bn| bn.iter(|| WGS84.to_ecef(black_box(g))));
    c.bench_function("ecef_to_geodetic", |bn| bn.iter(|| WGS84.to_geodetic(black_box(e))));
    c.bench_function("look_angles", |bn| bn.iter(|| frame.look_angles(black_box(target))));
    c.bench_function("geodesic_inverse", |bn| {
        bn.iter(|| WGS84.geodesic_inverse(black_box(lax), black_box(jfk)))
    });
    c.bench_function("geodesic_direct", |bn| {
        bn.iter(|| WGS84.geodesic_direct(black_box(lax), black_box(65.0), black_box(4e6)))
    });
    c.bench_function("helmert_apply", |bn| bn.iter(|| helmert.apply(black_box(e))));
    c.bench_function("ray_intersect", |bn| {
        bn.iter(|| WGS84.ray_intersect(black_box(e), black_box(Vec3::new(-1.0, 0.0, 0.0))))
    });
    c.bench_function("great_circle_distance", |bn| {
        bn.iter(|| WGS84.great_circle_distance_m(black_box(lax), black_box(jfk)))
    });
    c.bench_function("vec3_dot", |bn| bn.iter(|| black_box(a).dot(black_box(b))));
    c.bench_function("vec3_cross", |bn| bn.iter(|| black_box(a).cross(black_box(b))));
    c.bench_function("vec3_norm", |bn| bn.iter(|| black_box(a).norm()));
}

criterion_group! {
    name = benches;
    // Keep the run quick: short warm-up and measurement windows are plenty for ops this cheap.
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_millis(700))
        .sample_size(60);
    targets = bench
}
criterion_main!(benches);
