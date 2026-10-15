//! Watch the International Space Station pass over Brisbane: where to point the antenna, and the
//! Doppler shift to tune the radio for.
//!
//! The `sgp4` crate predicts the orbit from a two-line element set (TLE). It gives positions and
//! velocities in the TEME frame, which is geo3d's `Eci`.
//!
//! `cargo run --example satellite`

use geo3d::{eci, Eci, Geodetic, LocalFrame, Track, Vec3};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // The ISS on 12 July 2020.
    let tle = sgp4::Elements::from_tle(
        Some("ISS (ZARYA)".to_owned()),
        b"1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992",
        b"2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008",
    )?;
    let orbit = sgp4::Constants::from_elements(&tle)?;
    let epoch = tle.datetime.and_utc().timestamp_micros() as f64 / 1e6; // Unix seconds

    let station = LocalFrame::new(Geodetic::new(-27.47, 153.03, 30.0)); // Brisbane
    let radio_hz = 145.8e6; // the ISS's voice downlink

    for minute in 0..=95 {
        // SGP4 gives kilometres and kilometres per second.
        let state = orbit.propagate(sgp4::MinutesSinceEpoch(minute as f64))?;
        let [x, y, z] = state.position.map(|km| km * 1e3);
        let [vx, vy, vz] = state.velocity.map(|km_s| km_s * 1e3);

        // How far the Earth has turned then. UT1 was 0.2 s behind UTC that day (IERS Bulletin A).
        let turned = eci::gmst_from_unix_with_dut1(epoch + 60.0 * minute as f64, -0.2);
        let iss = Track::from_eci(Eci::new(x, y, z), Vec3::new(vx, vy, vz), turned);

        // Only while it's above the horizon.
        let look = station.look_angles_ecef(iss.position);
        if look.elevation_deg > 0.0 {
            let doppler_hz = -station.range_rate_mps(iss) * radio_hz / 299_792_458.0;
            println!("{minute} min: point at {look:.1}, Doppler {doppler_hz:+.0} Hz");
        }
    }
    Ok(())
}
