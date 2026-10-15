//! Plan a long flight: the shortest route from Sydney to Los Angeles, waypoints along it, and how
//! far a plane has strayed from it.
//!
//! `cargo run --example flight_route`

use geo3d::{Geodetic, WGS84};

fn main() {
    let sydney = Geodetic::new(-33.9461, 151.1772, 0.0);
    let los_angeles = Geodetic::new(33.9416, -118.4085, 0.0);

    // The shortest path over the Earth's surface: a geodesic.
    let route = WGS84.geodesic_inverse(sydney, los_angeles).unwrap();
    println!("Distance: {:.0} km", route.distance_m / 1_000.0);
    println!("Set off heading {:.1}°", route.initial_bearing_deg);

    // A waypoint every quarter of the way.
    for quarter in 1..=3 {
        let distance_m = route.distance_m * quarter as f64 / 4.0;
        let (waypoint, _heading) = WGS84.geodesic_direct(sydney, route.initial_bearing_deg, distance_m);
        println!("Waypoint {quarter}: {waypoint:#.2}");
    }

    // A plane reports its position. How far off the route is it, and how far along?
    let plane = Geodetic::new(-3.0, -169.0, 11_000.0);
    let off_m = WGS84.cross_track_distance_m(sydney, los_angeles, plane);
    let along_m = WGS84.along_track_distance_m(sydney, los_angeles, plane);
    let side = if off_m > 0.0 { "right" } else { "left" };
    println!(
        "The plane is {:.0} km {side} of the route, {:.0} km along it",
        off_m.abs() / 1_000.0,
        along_m / 1_000.0
    );
}
