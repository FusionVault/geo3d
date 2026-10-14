//! Track the International Space Station from a ground station.
//!
//! SGP4, the standard model for satellites published as two-line elements (TLEs), gives positions
//! in the TEME frame. That is geo3d's `Eci`, so they plug straight in: rotate them to Earth-fixed
//! with the sidereal angle, then point at them from a `LocalFrame`.
//!
//! `cargo run --example satellite`

use geo3d::{ecef_to_geodetic, eci, Eci, Geodetic, LocalFrame};

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

    let station = LocalFrame::new(Geodetic::new(-27.47, 153.03, 30.0)); // Brisbane
    println!("min   beneath the ISS          height    azimuth  elevation       range");
    for minutes in (0..=90).step_by(5) {
        let teme = orbit.propagate(sgp4::MinutesSinceEpoch(minutes as f64))?.position; // km
        let unix = epoch_unix + minutes as f64 * 60.0;

        let in_space = Eci::new(teme[0] * 1e3, teme[1] * 1e3, teme[2] * 1e3);
        let fixed = in_space.to_ecef(eci::gmst_from_unix_with_dut1(unix, dut1_s));
        let below = ecef_to_geodetic(fixed);
        let look = station.look_angles_ecef(fixed);

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
    Ok(())
}
