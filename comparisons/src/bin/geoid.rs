//! Geoid heights of geo3d, egm96 and egm2008 against the gravity models, as a Markdown table.
//!
//! Needs GeographicLib's geoid data: run `./fetch-geoids.sh`, then
//! `cargo run --release --features geoid --bin geoid > results/geoid.md`.

use geo3d_comparisons::report::{sci, table, wrap};
use geo3d_comparisons::*;

fn main() {
    geoid();
}
fn geoid() {
    // GeographicLib's test set: 500 000 random points with the geoid height of each model from NGA's
    // own harmonic synthesis, to 1 µm (latitude, longitude, EGM84, EGM96, EGM2008).
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/data/GeoidHeights.dat");
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e} (run fetch-geoids.sh)"));
    let truth: Vec<[f64; 5]> = text
        .lines()
        .map(|l| {
            let mut v = l.split_whitespace().map(|t| t.parse().unwrap());
            [(); 5].map(|_| v.next().unwrap())
        })
        .collect();
    type F = fn(f64, f64) -> Option<f64>;
    // Library, model and method, the function, the model's column, and every how-many-th point.
    let libs: [(&str, &str, F, usize, usize); 7] = [
        (
            "geo3d",
            "EGM96, egm96-5 grid (19 MB), cubic",
            geo3d_geoids::egm96_m,
            3,
            1,
        ),
        (
            "geo3d",
            "EGM96, egm96-5 grid, bilinear",
            geo3d_geoids::egm96_bilinear_m,
            3,
            1,
        ),
        (
            "geo3d",
            "EGM96, egm96-15 grid (2 MB), cubic",
            geo3d_geoids::egm96_15_m,
            3,
            1,
        ),
        (
            "geo3d",
            "EGM2008, egm2008-1 grid (470 MB), cubic",
            geo3d_geoids::egm2008_m,
            4,
            1,
        ),
        (
            "egm96",
            "EGM96, its 5′ raster (the same egm96-5 grid), bilinear: its default",
            with_egm96::geoid_m,
            3,
            1,
        ),
        (
            "egm96",
            "EGM96, spherical harmonics to degree 360",
            with_egm96::harmonics_m,
            3,
            50,
        ),
        ("egm2008", "EGM2008, its built-in 3° grid", with_egm2008::geoid_m, 4, 1),
    ];
    let mut rows = Vec::new();
    let mut harmonics_us = 0.0;
    for (name, what, f, column, every) in libs {
        let (mut count, mut failures, mut max, mut sum2) = (0, 0, 0f64, 0.0);
        let start = std::time::Instant::now();
        for t in truth.iter().step_by(every) {
            count += 1;
            // Longitudes from −180° to 180°: the test set's run to 360°, which egm2008 refuses.
            match f(t[0], wrap(t[1])).filter(|n| n.is_finite()) {
                Some(n) => {
                    max = max.max((n - t[column]).abs());
                    sum2 += (n - t[column]).powi(2);
                }
                None => failures += 1,
            }
        }
        if what.contains("harmonics") {
            harmonics_us = start.elapsed().as_secs_f64() * 1e6 / count as f64;
        }
        rows.push(vec![
            name.to_string(),
            what.to_string(),
            count.to_string(),
            failures.to_string(),
            sci(max),
            sci((sum2 / (count - failures) as f64).sqrt()),
        ]);
    }
    println!("## Geoid heights\n");
    println!("The geoid's height above the ellipsoid, against the gravity model each crate implements, at GeographicLib's 500 000 random test points (NGA's harmonic synthesis, to 1 µm). egm96's harmonic synthesis runs on every 50th point: it took {harmonics_us:.0} µs per point here. Longitudes are given from −180° to 180°; egm2008 refuses the same places written as 180° to 360°, and geo3d takes either.\n");
    table(
        "Geoid height against the model",
        &[
            "Library",
            "Model and method",
            "Points",
            "Failures",
            "Max error (m)",
            "RMS error (m)",
        ],
        &rows,
    );
}
