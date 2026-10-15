//! Measure a paddock from the GPS positions of its corners: its area, the length and heading of
//! each fence, and grid references for a map.
//!
//! `cargo run --example field_survey`

use geo3d::{Geodetic, WGS84};

fn main() {
    // The corners, walked anticlockwise with a GPS, near Dubbo in New South Wales.
    let corners = [
        Geodetic::new(-32.2430, 148.6010, 270.0),
        Geodetic::new(-32.2465, 148.6052, 268.0),
        Geodetic::new(-32.2431, 148.6101, 271.0),
        Geodetic::new(-32.2392, 148.6060, 275.0),
    ];

    let paddock = WGS84.polygon_area(&corners).unwrap();
    println!("Area:  {:.2} hectares", paddock.area_m2 / 10_000.0);
    println!("Fence: {:.0} m all the way round", paddock.perimeter_m);

    // Each fence line, and the compass heading to walk it.
    for i in 0..corners.len() {
        let (from, to) = (corners[i], corners[(i + 1) % corners.len()]);
        let fence = WGS84.geodesic_inverse(from, to).unwrap();
        println!(
            "Fence {}: {:.1} m, heading {:.1}°",
            i + 1,
            fence.distance_m,
            fence.initial_bearing_deg
        );
    }

    // Grid references for a topographic map (UTM zone, easting and northing), and the middle.
    for corner in corners {
        println!("Corner: {:.0}", WGS84.to_utm_ups(corner).unwrap());
    }
    println!("Middle: {:.5}", WGS84.mean_position(&corners).unwrap());
}
