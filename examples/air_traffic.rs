//! Air traffic control: which flights are inside the airport's controlled airspace, and will two of
//! them come too close?
//!
//! `cargo run --example air_traffic`

use geo3d::{AltitudeBand, Course, Cylinder, Geodetic, LocalFrame, Track, Volume, WGS84};

fn main() {
    // Controlled airspace: 15 km around Sydney Airport, from the ground up to 3 000 m.
    let airport = Geodetic::new(-33.9461, 151.1772, 21.0);
    let zone = Cylinder::new(airport, 15_000.0, AltitudeBand::below(3_000.0));

    let qf1 = Geodetic::new(-33.80, 151.30, 2_500.0); // coming in to land
    let va2 = Geodetic::new(-33.90, 151.25, 2_400.0); // climbing out
    println!("QF1 in controlled airspace: {}", zone.contains(qf1));
    println!("VA2 in controlled airspace: {}", zone.contains(va2));

    // Each flight as a track: where it is, and its heading, speed (m/s) and climb (m/s).
    let qf1 = flying(qf1, Course::new(225.0, 120.0, -4.0));
    let va2 = flying(va2, Course::new(0.0, 100.0, 6.0));

    // When will they be closest, and how close? Near airports, aircraft keep 3 nautical miles
    // (5 556 m) apart.
    let closest = qf1.closest_approach(va2);
    println!("Closest in {:.0} s, {:.0} m apart", closest.time_s, closest.distance_m);
    if qf1.conflict(va2, 5_556.0).is_some() {
        println!("Conflict: they'll come within 3 nautical miles");
    }

    // Where QF1 will be in a minute, if it holds its course.
    println!("QF1 in 60 s: {:.4}", WGS84.to_geodetic(qf1.at(60.0)));
}

/// A track from where an aircraft is and the course it's flying.
fn flying(at: Geodetic, course: Course) -> Track {
    Track::new(WGS84.to_ecef(at), LocalFrame::new(at).velocity_of_course(course))
}
