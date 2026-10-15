//! Search and rescue at sea: a yacht has lost power and is drifting. How soon can the rescue boat
//! reach it, which way should it steer, and from how far will the crew see the yacht's mast?
//!
//! `cargo run --example search_and_rescue`

use geo3d::{Geodetic, LocalFrame, Track, WGS84};

fn main() {
    // Two position reports from the yacht, 10 minutes apart, give its drift.
    let at_0900 = Geodetic::new(-33.9550, 151.4450, 0.0);
    let at_0910 = Geodetic::new(-33.9493, 151.4519, 0.0);
    let yacht = Track::from_fixes(WGS84.to_ecef(at_0900), WGS84.to_ecef(at_0910), 600.0);

    // The rescue boat leaves Watsons Bay, in Sydney Harbour, at 15 m/s (about 30 knots).
    let base = LocalFrame::new(Geodetic::new(-33.8418, 151.2810, 0.0));
    let rescue = yacht.intercept_from(base.origin_ecef(), 15.0).unwrap();

    // Which way to steer, and where they meet.
    let steer = base.course_of_velocity(rescue.velocity);
    let meet = WGS84.to_geodetic(rescue.position);
    println!("The yacht drifts at {:.1} m/s", yacht.velocity.norm());
    println!(
        "Reach it in {:.0} minutes, steering {:.0}°",
        rescue.time_s / 60.0,
        steer.course_deg
    );
    println!("Meet at {:.4}°, {:.4}°", meet.lat_deg, meet.lon_deg);

    // The crew look out from 3 m up, and the yacht's mast is 15 m tall. Each sees over the
    // horizon from its own height, so the mast shows above it from this far away:
    let sighting_m = WGS84.horizon_distance_m(3.0) + WGS84.horizon_distance_m(15.0);
    println!(
        "Its mast shows over the horizon from {:.1} km away",
        sighting_m / 1_000.0
    );
}
