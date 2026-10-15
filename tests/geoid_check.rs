//! Checks [`GeoidGrid`] on GeographicLib's real geoid grids, which are too big to ship.
//!
//! Download them from <https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst>, with
//! the test set `GeoidHeights.dat` (gunzipped) from the same page, and run
//!
//! ```sh
//! GEOGRAPHICLIB_GEOID_PATH=/path/to/geoids cargo test --release --test geoid_check -- --ignored --nocapture
//! ```
//!
//! The path defaults to `comparisons/data/geoids`, then GeographicLib's own default,
//! `/usr/local/share/GeographicLib/geoids`. `GeoidHeights.dat` goes in the same folder or the one
//! above it.

use std::path::PathBuf;

use geo3d::{Geoid, GeoidGrid, GeoidInterpolation};

fn dir() -> PathBuf {
    let candidates = [
        std::env::var_os("GEOGRAPHICLIB_GEOID_PATH").map(PathBuf::from),
        Some(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("comparisons/data/geoids")),
        Some(PathBuf::from("/usr/local/share/GeographicLib/geoids")),
    ];
    candidates
        .into_iter()
        .flatten()
        .find(|d| d.join("egm96-5.pgm").exists())
        .expect("no geoid grids found: see the top of tests/geoid_check.rs")
}

fn grid(name: &str) -> GeoidGrid<Vec<u8>> {
    let path = dir().join(format!("{name}.pgm"));
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    GeoidGrid::from_pgm(bytes).unwrap()
}

/// GeographicLib's online GeoidEval (version 2.5.2, cubic, egm2008-1, egm96-5 and egm84-15),
/// queried 2026-10-15: latitude, longitude, then EGM2008, EGM96 and EGM84 heights in metres,
/// rounded to 0.1 mm.
#[rustfmt::skip]
const GEOIDEVAL: [[f64; 5]; 26] = [
    [-33.87, 151.21, 22.3526, 22.4595, 22.1025],
    [16.776, -3.009, 28.7880, 28.7067, 31.2979],
    [90.0, 0.0, 14.8980, 13.6050, 13.0980],
    [-90.0, 0.0, -30.1500, -29.5350, -29.7120],
    [89.999, 45.0, 14.8991, 13.6049, 13.0991],
    [-89.999, -135.0, -30.1505, -29.5358, -29.7145],
    [0.0, 0.0, 17.2260, 17.1630, 18.3296],
    [0.0, 180.0, 21.2813, 21.1534, 21.7089],
    [0.0, -180.0, 21.2813, 21.1534, 21.7089],
    [27.9881, 86.925, -28.4939, -28.7408, -30.4587],
    [-0.5, 359.99, 17.4580, 17.3567, 18.4278],
    [5.5, 80.5, -103.3321, -103.1341, -103.3920],
    [-6.0, 145.0, 81.0990, 80.5192, 78.0223],
    [40.7128, -74.006, -32.7222, -32.7612, -33.4307],
    [51.4778, -0.0015, 45.8978, 45.7992, 47.2480],
    [64.1466, -21.9426, 66.4794, 66.4283, 67.2118],
    [-77.85, 166.67, -52.4732, -53.3846, -56.3638],
    [35.36, 138.73, 42.2707, 41.2849, 40.1702],
    [19.82, -155.47, 26.2207, 24.8180, 18.0423],
    [-54.8, -68.3, 13.2421, 13.0647, 12.4768],
    [71.29, -156.79, -0.8529, -1.2566, -0.1674],
    [1.3521, 103.8198, 8.0297, 7.6559, 7.3943],
    [36.5785, -118.2923, -25.1726, -26.5626, -25.9880],
    [11.35, 142.2, 34.3742, 34.7933, 40.7344],
    [60.0, 179.99, 0.9294, 0.4319, 2.0334],
    [-45.0, -0.001, 20.9880, 21.2432, 22.3297],
];

#[test]
#[ignore = "needs GeographicLib's geoid grids; see the top of this file"]
fn matches_geoideval() {
    for (column, name) in [(2, "egm2008-1"), (3, "egm96-5"), (4, "egm84-15")] {
        let g = grid(name);
        for row in GEOIDEVAL {
            let n = g.undulation_m(row[0], row[1]);
            // Equal once rounded to 0.1 mm, as GeoidEval prints it.
            assert!(
                (n - row[column]).abs() <= 0.5e-4 + 1e-9,
                "{name} at {}, {}: {n} vs {}",
                row[0],
                row[1],
                row[column]
            );
        }
    }
}

/// Against the gravity models themselves (NGA's harmonic synthesis, in GeographicLib's 500 000
/// point test set), every grid stays within the largest error its file states, to the 1 mm it's
/// stated to.
#[test]
#[ignore = "needs GeographicLib's geoid grids; see the top of this file"]
fn within_stated_errors_of_the_models() {
    let dir = dir();
    let path = [dir.join("GeoidHeights.dat"), dir.join("../GeoidHeights.dat")]
        .into_iter()
        .find(|p| p.exists())
        .expect("GeoidHeights.dat (gunzipped) not found");
    let text = std::fs::read_to_string(path).unwrap();
    let rows: Vec<[f64; 5]> = text
        .lines()
        .map(|l| {
            let mut v = l.split_whitespace().map(|t| t.parse::<f64>().unwrap());
            [(); 5].map(|_| v.next().unwrap())
        })
        .collect();
    assert_eq!(rows.len(), 500_000);
    // Column 2 is EGM84, 3 EGM96 and 4 EGM2008.
    for (column, name) in [
        (2, "egm84-15"),
        (3, "egm96-15"),
        (3, "egm96-5"),
        (4, "egm2008-5"),
        (4, "egm2008-2_5"),
        (4, "egm2008-1"),
    ] {
        let Ok(bytes) = std::fs::read(dir.join(format!("{name}.pgm"))) else {
            println!("{name}: not downloaded, skipped");
            continue;
        };
        for interpolation in [GeoidInterpolation::Cubic, GeoidInterpolation::Bilinear] {
            let g = GeoidGrid::from_pgm(bytes.as_slice())
                .unwrap()
                .with_interpolation(interpolation);
            let (mut max, mut sum2) = (0f64, 0.0);
            for r in &rows {
                let e = g.undulation_m(r[0], r[1]) - r[column];
                max = max.max(e.abs());
                sum2 += e * e;
            }
            let rms = (sum2 / rows.len() as f64).sqrt();
            let stated = g.max_error_m().unwrap();
            println!(
                "{name:12} {interpolation:?}: max {:.2} mm, RMS {:.3} mm (file states max {} mm)",
                max * 1e3,
                rms * 1e3,
                stated * 1e3
            );
            // The files state errors to the nearest 1 mm.
            assert!(max <= stated + 0.5e-3, "{name} {interpolation:?}: {max} > {stated}");
        }
    }
}
