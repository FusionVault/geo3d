//! Where is a drone's camera looking? From the drone's GPS position and its yaw, pitch and roll (as
//! the autopilot reports them), find the camera's line of sight, and where the thing its laser
//! rangefinder hits is.
//!
//! `cargo run --example drone_camera`

use geo3d::{Attitude, Frd, Geodetic, LocalFrame};

fn main() {
    // The drone, 120 m above a paddock, flying east with its nose 5° down and rolled 2° right.
    // (Autopilots like PX4 send a quaternion: use `Attitude::from_quaternion` for that.)
    let drone = LocalFrame::new(Geodetic::new(-32.2430, 148.6010, 390.0));
    let attitude = Attitude::new(90.0, -5.0, 2.0);

    // The camera points 30° below the drone's nose: forward and down, in the drone's own frame.
    let (down, forward) = 30f64.to_radians().sin_cos();
    let camera = Frd::new(forward, 0.0, down);

    // That direction as a compass bearing and an angle below the horizon.
    let look = attitude.body_to_ned(camera).to_aer();
    println!(
        "Camera points at {:.1}°, {:.1}° below the horizon",
        look.azimuth_deg, -look.elevation_deg
    );

    // The laser rangefinder reads 209 m along that line. Where is the spot it hit?
    let target = drone.geodetic_of_body(attitude, camera * 209.0);
    println!("Target: {target:.6}");
}
