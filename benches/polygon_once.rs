//! A polygon used once, as each row of a column of polygons is: `polygon_containment`, which has no
//! index to build, against `geo`'s `Contains` and against preparing a `Footprint` for the one test.
//! Each iteration tests 4,096 points, so the time per element is the time per test. Run with
//! `cargo bench --bench polygon_once --features footprint,geo-types`.

use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use geo::{Contains, Coord, LineString, Point, Polygon};
use geo3d::{polygon_containment, Containment, Footprint, Geodetic};

/// A wavy ring of `n` vertices around (151, −33), radius 1°, closed.
fn wavy(n: usize) -> Vec<Coord<f64>> {
    let mut v: Vec<Coord<f64>> = (0..n)
        .map(|k| {
            let a = std::f64::consts::TAU * k as f64 / n as f64;
            let r = 1.0 + 0.05 * (a * 40.0).sin();
            Coord {
                x: 151.0 + r * a.cos(),
                y: -33.0 + r * a.sin(),
            }
        })
        .collect();
    v.push(v[0]);
    v
}

/// 4,096 points uniform in the 2° box around the ring.
fn points() -> Vec<Coord<f64>> {
    let mut s = 42u64;
    let mut next = || {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..4096)
        .map(|_| Coord {
            x: 150.0 + 2.0 * next(),
            y: -34.0 + 2.0 * next(),
        })
        .collect()
}

fn bench(c: &mut Criterion) {
    let pts = points();
    for n in [4, 50, 5_000] {
        let ring = wavy(n);
        let polygon = Polygon::new(LineString::new(ring.clone()), vec![]);
        let at: Vec<Geodetic> = ring.iter().map(|p| Geodetic::from_lon_lat(p.x, p.y, 0.0)).collect();
        let point = |p: &Coord<f64>| Geodetic::from_lon_lat(p.x, p.y, 0.0);
        // every method agrees on these points before any is timed
        let count = |f: &dyn Fn(&Coord<f64>) -> bool| pts.iter().filter(|p| f(p)).count();
        let expected = count(&|p| polygon.contains(&Point::from(*p)));
        assert!(expected > 0);
        assert_eq!(
            count(&|p| polygon_containment([&at[..]], point(p)) == Containment::Inside),
            expected
        );
        assert_eq!(
            count(&|p| Footprint::polygon(&at, &[]).containment(point(p)) == Containment::Inside),
            expected
        );

        let mut g = c.benchmark_group(format!("used_once_{n}_vertices"));
        g.throughput(Throughput::Elements(pts.len() as u64));
        g.bench_function("polygon_containment", |bn| {
            bn.iter(|| {
                pts.iter()
                    .filter(|p| polygon_containment([&at[..]], point(p)) == Containment::Inside)
                    .count()
            })
        });
        g.bench_function("geo_contains", |bn| {
            bn.iter(|| pts.iter().filter(|p| polygon.contains(&Point::from(**p))).count())
        });
        g.bench_function("footprint_prepared_per_test", |bn| {
            bn.iter(|| {
                pts.iter()
                    .filter(|p| Footprint::polygon(&at, &[]).containment(point(p)) == Containment::Inside)
                    .count()
            })
        });
        g.finish();
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(2))
        .sample_size(30);
    targets = bench
}
criterion_main!(benches);
