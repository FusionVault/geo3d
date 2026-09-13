//! A ground station tracking an aircraft: slant range, look angles, the ground point under its
//! line of sight, and a proximity/conflict check against a second aircraft.
//!
//! `cargo run --example tracking`

use geo3d::{distance3d_m, eci, Cone, Course, Geodetic, LocalFrame, Sphere, Track, Volume, WGS84};

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

    // Conflict check: the aircraft heading south-west at 120 m/s vs another heading north at 100 m/s,
    // each expressed as a course over ground at its own position and turned into an ECEF velocity.
    let other = Geodetic::new(-33.90, 151.25, 2_400.0);
    let a = Track::new(
        WGS84.to_ecef(aircraft),
        LocalFrame::new(aircraft).velocity_of_course(Course::new(225.0, 120.0, 0.0)),
    );
    let b = Track::new(
        WGS84.to_ecef(other),
        LocalFrame::new(other).velocity_of_course(Course::new(0.0, 100.0, 0.0)),
    );
    let cpa = a.closest_approach(b);
    println!(
        "closest approach in {:.0} s at {:.0} m; conflict inside 1 km: {}",
        cpa.time_s,
        cpa.distance_m,
        a.conflict(b, 1_000.0).is_some()
    );
    println!("aircraft in 60 s (dead reckoning): {:?}", WGS84.to_geodetic(a.at(60.0)));

    // Line of sight to a satellite and its inertial position at a given time.
    let sat_ecef = WGS84.to_ecef(Geodetic::new(-20.0, 140.0, 500_000.0));
    println!(
        "station sees the satellite: {}",
        WGS84.line_of_sight(WGS84.to_ecef(station), sat_ecef)
    );
    let gmst = eci::gmst_from_unix(1_800_000_000.0);
    println!("satellite in ECI at that time: {:?}", sat_ecef.to_eci(gmst));
    println!(
        "great-circle to the other aircraft: {:.0} m, bearing {:.1}°",
        WGS84.great_circle_distance_m(aircraft, other),
        WGS84.initial_bearing_deg(aircraft, other)
    );
}
