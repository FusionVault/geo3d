//! Turn GPS heights into heights above sea level, with a geoid grid from GeographicLib.
//!
//! Download a grid, such as egm96-5 (19 MB), from
//! <https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst>, unpack it, then
//!
//! `cargo run --example sea_level -- path/to/egm96-5.pgm`

use geo3d::{Geodetic, Geoid, GeoidGrid};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let Some(path) = std::env::args().nth(1) else {
        println!("Usage: cargo run --example sea_level -- path/to/egm96-5.pgm");
        println!("Grids: https://geographiclib.sourceforge.io/C++/doc/geoid.html#geoidinst");
        return Ok(());
    };
    let geoid = GeoidGrid::from_pgm(std::fs::read(path)?)?;
    println!("Using {}", geoid.description().unwrap_or("a geoid grid"));

    // GPS readings: heights above the WGS84 ellipsoid.
    let readings = [
        ("Sydney Opera House steps", Geodetic::new(-33.8568, 151.2153, 25.0)),
        ("The summit of Everest", Geodetic::new(27.9881, 86.9250, 8_820.1)),
        ("The shore of the Dead Sea", Geodetic::new(31.5590, 35.4732, -412.0)),
    ];
    for (place, gps) in readings {
        let above_sea = geoid.sea_level_height_m(gps);
        println!(
            "{place}: GPS height {:.1} m, height above sea level {above_sea:.1} m",
            gps.height_m
        );
    }
    Ok(())
}
