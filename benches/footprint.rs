//! A prepared footprint against the usual alternatives: `geo`'s `Contains`, and a plain scan of
//! every edge. Each iteration tests the same 4,096 points, so the time per element is the time per
//! point. Run with `cargo bench --bench footprint --features footprint,geo-types`.

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use geo::{Contains, Coord, LineString, Point, Polygon};
use geo3d::{Containment, Footprint, Geodetic};

/// A wavy ring of `n` vertices around (151, −33), radius 1°, like a coastline.
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

/// The even-odd rule over every edge: what a footprint test costs without an index.
fn scan(ring: &[Coord<f64>], c: Coord<f64>) -> bool {
    let mut inside = false;
    for w in ring.windows(2) {
        let (a, b) = (w[0], w[1]);
        if (a.y > c.y) != (b.y > c.y) && c.x < (b.x - a.x) * (c.y - a.y) / (b.y - a.y) + a.x {
            inside = !inside;
        }
    }
    inside
}

fn bench(c: &mut Criterion) {
    let pts = points();
    for n in [50, 5_000] {
        let ring = wavy(n);
        let polygon = Polygon::new(LineString::new(ring.clone()), vec![]);
        let footprint = Footprint::from(&polygon);
        // every method agrees on these points before any is timed
        let count = |f: &dyn Fn(Coord<f64>) -> bool| pts.iter().filter(|p| f(**p)).count();
        let expected = count(&|p| polygon.contains(&Point::from(p)));
        assert_eq!(
            count(&|p| footprint.containment(Geodetic::from_lon_lat(p.x, p.y, 0.0)) == Containment::Inside),
            expected
        );
        assert_eq!(count(&|p| scan(&ring, p)), expected);

        let mut g = c.benchmark_group(format!("footprint_{n}_vertices"));
        g.throughput(Throughput::Elements(pts.len() as u64));
        g.bench_function("prepared", |bn| {
            bn.iter(|| {
                pts.iter()
                    .filter(|p| footprint.containment(Geodetic::from_lon_lat(p.x, p.y, 0.0)) == Containment::Inside)
                    .count()
            })
        });
        g.bench_function("geo_contains", |bn| {
            bn.iter(|| pts.iter().filter(|p| polygon.contains(&Point::from(**p))).count())
        });
        g.bench_function("edge_scan", |bn| {
            bn.iter(|| pts.iter().filter(|p| scan(black_box(&ring), **p)).count())
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
