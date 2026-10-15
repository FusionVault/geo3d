//! An airport radar tracking a plane: which way to point, how fast the plane is closing (what the
//! radar's Doppler measures), and how sure we can be of where it is.
//!
//! `cargo run --example radar`

use geo3d::{Course, Cov3, Geodetic, LocalFrame, Track, WGS84};

fn main() {
    // The radar at Sydney Airport, and a plane coming in from the north-east.
    let radar = LocalFrame::new(Geodetic::new(-33.9461, 151.1772, 21.0));
    let plane = Geodetic::new(-33.80, 151.30, 2_500.0);

    // Which way to point the dish.
    let look = radar.look_angles(plane);
    println!("Point at:  {look:.1}");

    // The plane broadcasts its heading, speed and climb (ADS-B). Turn them into a 3D velocity.
    let velocity = LocalFrame::new(plane).velocity_of_course(Course::new(225.0, 120.0, -4.0));
    let track = Track::new(WGS84.to_ecef(plane), velocity);

    // How fast it's closing, and the Doppler shift that gives on a 2.8 GHz radar.
    let range_rate = radar.range_rate_mps(track); // negative while it comes closer
    let doppler_hz = -2.0 * range_rate * 2.8e9 / 299_792_458.0;
    println!(
        "Closing:   {:.0} m/s, a Doppler shift of {doppler_hz:+.0} Hz",
        -range_rate
    );

    // The radar measures angles to ±0.1° and range to ±50 m. How far off could the plane be?
    let noise = Cov3::from_std_devs([0.1, 0.1, 50.0]); // azimuth °, elevation °, range m
    let [east, north, up] = look.cov_to_enu(noise).std_devs();
    println!("Accuracy:  ±{east:.0} m east, ±{north:.0} m north, ±{up:.0} m up");
}
