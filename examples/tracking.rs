//! A ground station tracking an aircraft: slant range, look angles, the ground point under its
//! line of sight, and a proximity/conflict check against a second aircraft.
//!
//! `cargo run --example tracking`

use geo3d::{closest_approach, distance3d_m, Cone, Geodetic, LocalFrame, Sphere, Vec3, Volume, WGS84};

fn main() {
    let station = Geodetic::new(-33.94, 151.18, 21.0); // ellipsoidal height, metres
    let aircraft = Geodetic::new(-33.80, 151.30, 2_500.0);

    println!("slant range: {:.1} m", distance3d_m(station, aircraft));

    let frame = LocalFrame::new(station);
    let aer = frame.look_angles(aircraft);
    println!(
        "look angles: az {:.2}°, el {:.2}°, range {:.1} m",
        aer.azimuth_deg, aer.elevation_deg, aer.range_m
    );
    println!("destination of that pointing: {:?}", frame.destination(aer));

    match frame.ground_strike(aer.azimuth_deg, -2.0) {
        Some(g) => println!(
            "a 2° down-look along that bearing meets the surface at {:.5}, {:.5}",
            g.lat_deg, g.lon_deg
        ),
        None => println!("that bearing never meets the surface"),
    }

    let approach_cone = Cone::new(station, aer.azimuth_deg, 5.0, 15.0, 40_000.0);
    println!("in the 15° approach cone: {}", approach_cone.contains(aircraft));
    println!(
        "within 30 km bubble: {}",
        Sphere::new(station, 30_000.0).contains(aircraft)
    );

    // Conflict check: the aircraft heading south-west at 120 m/s vs another heading north at 100 m/s.
    let other = Geodetic::new(-33.90, 151.25, 2_400.0);
    let (a, b) = (WGS84.to_ecef(aircraft), WGS84.to_ecef(other));
    let va = frame.dir_to_ecef(geo3d::Aer::new(225.0, 0.0, 120.0).to_enu());
    let vb = frame.dir_to_ecef(geo3d::Aer::new(0.0, 0.0, 100.0).to_enu());
    let cpa = closest_approach(a, va, b, vb);
    println!("closest approach in {:.0} s at {:.0} m", cpa.time_s, cpa.distance_m);
    let _ = Vec3::ZERO;
}
