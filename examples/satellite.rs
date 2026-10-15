//! Track the International Space Station from a ground station.
//!
//! SGP4, the standard model for satellites published as two-line elements (TLEs), gives positions
//! and velocities in the TEME frame. That is geo3d's `Eci`, so they plug straight in:
//! `Track::from_eci` turns them Earth-fixed with the sidereal angle. Then a `LocalFrame` gives the
//! way to point, and the range rate gives the Doppler shift to tune for.
//!
//! `cargo run --example satellite`

use geo3d::{ecef_to_geodetic, eci, Eci, Geodetic, LocalFrame, Track, Vec3};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // The ISS on 12 July 2020, as a two-line element set.
    let elements = sgp4::Elements::from_tle(
        Some("ISS (ZARYA)".to_owned()),
        b"1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992",
        b"2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008",
    )?;
    let orbit = sgp4::Constants::from_elements(&elements)?;
    let epoch = elements.datetime.and_utc();
    let epoch_unix = epoch.timestamp() as f64 + f64::from(epoch.timestamp_subsec_nanos()) * 1e-9;

    // DUT1 = UT1 − UTC for that date, from IERS Bulletin A (about −0.2 s in mid-2020). Leaving it
    // at 0 would put the satellite up to ~150 m out along its track.
    let dut1_s = -0.2;

    // Where the ISS is, `minutes` after the epoch, as an Earth-fixed track.
    let iss_at = |minutes: f64| -> Result<Track, sgp4::Error> {
        let state = orbit.propagate(sgp4::MinutesSinceEpoch(minutes))?; // km and km/s
        let [x, y, z] = state.position.map(|km| km * 1e3);
        let [vx, vy, vz] = state.velocity.map(|km_s| km_s * 1e3);
        let gmst = eci::gmst_from_unix_with_dut1(epoch_unix + minutes * 60.0, dut1_s);
        Ok(Track::from_eci(Eci::new(x, y, z), Vec3::new(vx, vy, vz), gmst))
    };

    let station = LocalFrame::new(Geodetic::new(-27.47, 153.03, 30.0)); // Brisbane
    println!("min   beneath the ISS          height    azimuth  elevation       range");
    for minutes in (0..=90).step_by(5) {
        let iss = iss_at(minutes as f64)?;
        let below = ecef_to_geodetic(iss.position);
        let look = station.look_angles_ecef(iss.position);
        println!(
            "{minutes:3}  {:7.2}°, {:8.2}°  {:6.1} km  {:7.1}°  {:8.1}°  {:8.0} km{}",
            below.lat_deg,
            below.lon_deg,
            below.height_m / 1e3,
            look.azimuth_deg,
            look.elevation_deg,
            look.range_m / 1e3,
            if look.elevation_deg > 10.0 { "  ← visible" } else { "" }
        );
    }

    // The pass over Brisbane, minute by minute, with the Doppler shift on the ISS's 145.8 MHz
    // voice downlink: it's one way, so the shift is −range rate × frequency / c.
    let (downlink_hz, c) = (145.8e6, 299_792_458.0);
    println!("\nThe pass:\nmin   elevation       range   range rate   Doppler at 145.8 MHz");
    for minutes in 46..=56 {
        let iss = iss_at(minutes as f64)?;
        let look = station.look_angles_ecef(iss.position);
        let rate = station.range_rate_mps(iss);
        // A check: it's the slope of the range, from SGP4's positions half a second either side.
        // (SGP4's velocity itself differs from the slope of its positions by about 1.5 cm/s.)
        let range = |m: f64| -> Result<f64, sgp4::Error> { Ok(station.look_angles_ecef(iss_at(m)?.position).range_m) };
        let slope = range(minutes as f64 + 0.5 / 60.0)? - range(minutes as f64 - 0.5 / 60.0)?; // over 1 s
        assert!((rate - slope).abs() < 0.05, "{rate} vs {slope}");
        println!(
            "{minutes:3}  {:8.1}°  {:8.0} km  {:7.2} km/s  {:+9.2} kHz",
            look.elevation_deg,
            look.range_m / 1e3,
            rate / 1e3,
            -rate * downlink_hz / c / 1e3
        );
    }
    Ok(())
}
