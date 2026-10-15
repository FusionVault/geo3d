//! Moving objects: [`Course`] (heading, speed, climb), [`Track`] (position and velocity),
//! [`closest_approach`] for near-miss checks, range rate for Doppler, and [`Intercept`]s. Motion is
//! assumed constant: straight lines at steady speed.

use crate::local::{wrap_360, Enu};
#[cfg(all(feature = "libm", not(feature = "std"), not(test)))]
use crate::math::Float as _;
use crate::position::Ecef;
use crate::vector::{Dd, Vec3};

/// Motion over the ground at a position: course in degrees clockwise from north, ground speed and
/// climb rate in m/s — what ADS-B, AIS and GNSS receivers report. Turn it into an ECEF velocity at
/// a position with [`LocalFrame::velocity_of_course`](crate::LocalFrame::velocity_of_course).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Course {
    /// Track over ground, degrees clockwise from north in `[0, 360)`.
    pub course_deg: f64,
    /// Horizontal speed, m/s.
    pub ground_speed_mps: f64,
    /// Vertical rate, m/s, positive up.
    pub climb_mps: f64,
}

impl Course {
    /// Build from a course in degrees and speeds in m/s.
    #[inline]
    pub const fn new(course_deg: f64, ground_speed_mps: f64, climb_mps: f64) -> Self {
        Course {
            course_deg,
            ground_speed_mps,
            climb_mps,
        }
    }

    /// From a velocity in a local ENU frame (m/s). Zero horizontal speed gives course 0.
    pub fn from_enu_velocity(v: Enu) -> Course {
        let ground = (v.east * v.east + v.north * v.north).sqrt();
        let course_deg = if ground == 0.0 {
            0.0
        } else {
            wrap_360(v.east.atan2(v.north).to_degrees())
        };
        Course {
            course_deg,
            ground_speed_mps: ground,
            climb_mps: v.up,
        }
    }

    /// As a velocity in a local ENU frame (m/s).
    #[inline]
    pub fn to_enu_velocity(self) -> Enu {
        let (s, c) = self.course_deg.to_radians().sin_cos();
        Enu::new(s * self.ground_speed_mps, c * self.ground_speed_mps, self.climb_mps)
    }

    /// Total speed (horizontal and vertical together), m/s.
    #[inline]
    pub fn speed_mps(self) -> f64 {
        (self.ground_speed_mps * self.ground_speed_mps + self.climb_mps * self.climb_mps).sqrt()
    }
}

/// A moving object: ECEF position (metres) and velocity (m/s), assumed constant.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Track {
    /// Position now, metres.
    pub position: Ecef,
    /// Velocity, m/s.
    pub velocity: Vec3,
}

impl Track {
    /// A track from a position and velocity.
    #[inline]
    pub const fn new(position: Ecef, velocity: Vec3) -> Self {
        Track { position, velocity }
    }

    /// A track from two position fixes `dt_s` seconds apart (the second fix is "now"). A
    /// non-positive `dt_s` gives zero velocity.
    pub fn from_fixes(earlier: Ecef, now: Ecef, dt_s: f64) -> Self {
        let velocity = if dt_s > 0.0 {
            (now.vec() - earlier.vec()) / dt_s
        } else {
            Vec3::ZERO
        };
        Track {
            position: now,
            velocity,
        }
    }

    /// Where the object will be `dt_s` seconds from now (negative for the past) — dead reckoning.
    #[inline]
    pub fn at(self, dt_s: f64) -> Ecef {
        self.position.offset(self.velocity * dt_s)
    }

    /// When and how close this track gets to another.
    #[inline]
    pub fn closest_approach(self, other: Track) -> Cpa {
        closest_approach(self.position, self.velocity, other.position, other.velocity)
    }

    /// A near-miss check: `Some(cpa)` if the two tracks come within `threshold_m` of each other now
    /// or in the future, `None` otherwise.
    pub fn conflict(self, other: Track, threshold_m: f64) -> Option<Cpa> {
        let cpa = self.closest_approach(other);
        (cpa.time_s >= 0.0 && cpa.distance_m <= threshold_m).then_some(cpa)
    }

    /// How fast the distance to `other` is changing now, m/s: positive while they move apart,
    /// negative while they close. It's what a radar's Doppler measures.
    ///
    /// The Doppler shift at frequency `f` follows from it: `−2 · range_rate · f / c` for a radar
    /// (there and back), and `−range_rate · f / c` for a one-way signal, such as a satellite's
    /// downlink. With `c` = 299 792 458 m/s.
    ///
    /// At zero range it gives the relative speed: from then on, they move apart.
    ///
    /// ```
    /// use geo3d::{Ecef, Track, Vec3};
    /// let radar = Track::new(Ecef::new(0.0, 0.0, 0.0), Vec3::ZERO);
    /// // 10 km away, flying at 250 m/s: 200 m/s towards the radar, 150 m/s across.
    /// let plane = Track::new(Ecef::new(10_000.0, 0.0, 0.0), Vec3::new(-200.0, 150.0, 0.0));
    /// let rate = radar.range_rate_mps(plane);
    /// assert!((rate + 200.0).abs() < 1e-9); // closing at 200 m/s
    /// let doppler_hz = -2.0 * rate * 3e9 / 299_792_458.0; // an S-band radar hears +4 kHz
    /// assert!((doppler_hz - 4002.8).abs() < 0.1);
    /// ```
    #[doc(alias = "doppler")]
    #[doc(alias = "closing_speed")]
    pub fn range_rate_mps(self, other: Track) -> f64 {
        let r = other.position.vec() - self.position.vec();
        let v = other.velocity - self.velocity;
        let range = r.norm();
        if range == 0.0 {
            v.norm()
        } else {
            r.dot(v) / range
        }
    }

    /// The soonest a pursuer at `from` can reach this track, flying straight at `speed_mps`: when,
    /// where, and which way to go. `None` if it never can: when this track is at least as fast and
    /// not coming towards it, or slower than [`slowest_intercept_from`](Self::slowest_intercept_from)
    /// needs.
    ///
    /// The answer is exact to within a rounding or two, even when the two speeds are equal, where
    /// the textbook formula divides by zero, or nearly so.
    ///
    /// ```
    /// use geo3d::{Ecef, Track, Vec3};
    /// // A target 3 km east, heading north at 300 m/s; a pursuer at 500 m/s.
    /// let target = Track::new(Ecef::new(3_000.0, 0.0, 0.0), Vec3::new(0.0, 300.0, 0.0));
    /// let hit = target.intercept_from(Ecef::new(0.0, 0.0, 0.0), 500.0).unwrap();
    /// assert!((hit.time_s - 7.5).abs() < 1e-12); // meets it 2 250 m north, in 7.5 s
    /// assert!((hit.velocity.x - 400.0).abs() < 1e-9 && (hit.velocity.y - 300.0).abs() < 1e-9);
    /// ```
    pub fn intercept_from(self, from: Ecef, speed_mps: f64) -> Option<Intercept> {
        // Solve |d + v·t| = s·t, that is (v·v − s²) t² + 2 (d·v) t + d·d = 0, for the first t ≥ 0.
        let d = self.position.vec() - from.vec();
        let v = self.velocity;
        let c = Dd::dot(d, d);
        if c.0 == 0.0 {
            return Some(Intercept::now(self.position));
        }
        if speed_mps.is_nan() || speed_mps < 0.0 {
            return None;
        }
        // Every term to 32 digits: near equal speeds, or a grazing intercept, they cancel.
        let b = Dd::dot(d, v);
        let a = Dd::dot(v, v) - Dd::prod(speed_mps, speed_mps);
        let disc = (b * b - a * c).f64();
        if disc.is_nan() || disc < 0.0 {
            return None;
        }
        let (a, b, c) = (a.f64(), b.f64(), c.f64());
        // Each root in the form that doesn't cancel.
        let t = if b <= 0.0 {
            c / (disc.sqrt() - b)
        } else if a < 0.0 {
            -(b + disc.sqrt()) / a
        } else {
            return None; // moving away, at least as fast
        };
        Intercept::at(self, d, t)
    }

    /// The slowest a pursuer at `from` can fly straight and still reach this track, with when and
    /// where it gets there. The speed is [`Intercept::speed_mps`]: the part of this track's velocity
    /// across the line between them.
    ///
    /// `None` if this track isn't coming towards `from`: then any speed faster than its own
    /// catches it in the end, and nothing slower does.
    ///
    /// ```
    /// use geo3d::{Ecef, Track, Vec3};
    /// // A ship 10 km north, heading south-east at 10 m/s.
    /// let ship = Track::new(Ecef::new(0.0, 10_000.0, 0.0), Vec3::new(8.0, -6.0, 0.0));
    /// let slowest = ship.slowest_intercept_from(Ecef::new(0.0, 0.0, 0.0)).unwrap();
    /// assert!((slowest.speed_mps() - 8.0).abs() < 1e-12); // 8 m/s, heading due east
    /// assert!((slowest.time_s - 10_000.0 / 6.0).abs() < 1e-9);
    /// ```
    pub fn slowest_intercept_from(self, from: Ecef) -> Option<Intercept> {
        let d = self.position.vec() - from.vec();
        let (c, b) = (Dd::dot(d, d).f64(), Dd::dot(d, self.velocity).f64());
        if c == 0.0 {
            return Some(Intercept::now(self.position));
        }
        if b.is_nan() || b >= 0.0 {
            return None;
        }
        // The slowest speed needs the longest time: when the discriminant is zero, t = −d·d / d·v.
        Intercept::at(self, d, -c / b)
    }
}

/// How a pursuer flying straight at a steady speed meets a moving target: from
/// [`Track::intercept_from`] or [`Track::slowest_intercept_from`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Intercept {
    /// Seconds from now until they meet.
    pub time_s: f64,
    /// Where they meet: ECEF, metres.
    pub position: Ecef,
    /// The pursuer's velocity, m/s: the way to fly. Zero if it's already there. Turn it into a
    /// course over the ground with [`LocalFrame::course_of_velocity`](crate::LocalFrame::course_of_velocity).
    pub velocity: Vec3,
}

impl Intercept {
    /// The pursuer's speed, m/s.
    #[inline]
    pub fn speed_mps(&self) -> f64 {
        self.velocity.norm()
    }

    fn now(position: Ecef) -> Intercept {
        Intercept {
            time_s: 0.0,
            position,
            velocity: Vec3::ZERO,
        }
    }

    /// Meeting `target` at time `t`, starting `d` from it; `None` unless `t` is a time to come.
    fn at(target: Track, d: Vec3, t: f64) -> Option<Intercept> {
        (t.is_finite() && t > 0.0).then(|| Intercept {
            time_s: t,
            position: target.at(t),
            velocity: (d + target.velocity * t) / t,
        })
    }
}

/// The closest point of approach of two moving objects: when, and how close.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Cpa {
    /// Seconds from now until the closest approach; negative when it is already past, zero when
    /// the two move in parallel.
    pub time_s: f64,
    /// The separation at that moment, metres.
    pub distance_m: f64,
}

/// When and how close two objects get, from their ECEF positions (metres) and velocities (m/s),
/// assuming both keep going in a straight line at steady speed. Compare `distance_m` with a
/// separation threshold for aircraft, vessels or rendezvous — or use [`Track::conflict`].
///
/// ```
/// use geo3d::{closest_approach, Ecef, Vec3};
/// // Two objects converging along x, offset 5 m in y: closest in 5 s, 5 m apart.
/// let cpa = closest_approach(
///     Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0),
///     Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0),
/// );
/// assert!((cpa.time_s - 5.0).abs() < 1e-9 && (cpa.distance_m - 5.0).abs() < 1e-9);
/// ```
#[doc(alias = "cpa")]
pub fn closest_approach(p1: Ecef, v1: Vec3, p2: Ecef, v2: Vec3) -> Cpa {
    let prel = p1.vec() - p2.vec();
    let vrel = v1 - v2;
    let vv = vrel.dot(vrel);
    let time_s = if vv == 0.0 { 0.0 } else { -prel.dot(vrel) / vv };
    Cpa {
        time_s,
        distance_m: (prel + vrel * time_s).norm(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_util::close;

    #[test]
    fn parallel_non_closing_is_now() {
        let c = closest_approach(
            Ecef::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Ecef::new(0.0, 3.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        );
        assert_eq!(c.time_s, 0.0);
        assert!(close(c.distance_m, 3.0, 1e-9));
    }

    #[test]
    fn diverging_reports_a_past_time() {
        let c = closest_approach(
            Ecef::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Ecef::new(-10.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
        );
        assert!(c.time_s < 0.0 && c.distance_m < 1e-9);
    }

    #[test]
    fn tracks_dead_reckon_and_detect_conflicts() {
        let a = Track::new(Ecef::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(a.at(10.0), Ecef::new(10.0, 0.0, 0.0));
        assert_eq!(a.at(-2.0), Ecef::new(-2.0, 0.0, 0.0));
        let b = Track::new(Ecef::new(10.0, 5.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        let cpa = a.closest_approach(b);
        assert!(close(cpa.time_s, 5.0, 1e-9) && close(cpa.distance_m, 5.0, 1e-9));
        assert!(a.conflict(b, 6.0).is_some());
        assert!(a.conflict(b, 4.0).is_none());
        // A closest approach in the past is not a conflict.
        let past = Track::new(Ecef::new(-10.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
        assert!(a.conflict(past, 100.0).is_none());
    }

    #[test]
    fn from_fixes_infers_velocity() {
        let t = Track::from_fixes(Ecef::new(0.0, 0.0, 0.0), Ecef::new(20.0, 0.0, 0.0), 10.0);
        assert_eq!(t.velocity, Vec3::new(2.0, 0.0, 0.0));
        assert_eq!(t.position, Ecef::new(20.0, 0.0, 0.0));
        assert_eq!(
            Track::from_fixes(Ecef::new(0.0, 0.0, 0.0), Ecef::new(20.0, 0.0, 0.0), 0.0).velocity,
            Vec3::ZERO
        );
    }

    /// xorshift64*, in `[0, 1)`.
    fn rnd(s: &mut u64) -> f64 {
        *s ^= *s >> 12;
        *s ^= *s << 25;
        *s ^= *s >> 27;
        (s.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A direction spread evenly over the sphere.
    fn unit(s: &mut u64) -> Vec3 {
        let z = 2.0 * rnd(s) - 1.0;
        let a = core::f64::consts::TAU * rnd(s);
        let r = (1.0 - z * z).sqrt();
        Vec3::new(r * a.cos(), r * a.sin(), z)
    }

    #[test]
    fn range_rate_is_how_fast_the_distance_changes() {
        let mut seed = 0x9e37_79b9_7f4a_7c15;
        for _ in 0..20_000 {
            let a = Track::new(
                Ecef::from(unit(&mut seed) * 6.4e6),
                unit(&mut seed) * 10f64.powf(4.0 * rnd(&mut seed)),
            );
            let b = Track::new(
                a.position
                    .offset(unit(&mut seed) * 10f64.powf(1.0 + 6.0 * rnd(&mut seed))),
                unit(&mut seed) * 10f64.powf(4.0 * rnd(&mut seed)),
            );
            // A central difference of the separation, kept relative so nothing large cancels.
            let (r, v) = (b.position.vec() - a.position.vec(), b.velocity - a.velocity);
            let h = 1e-4 * r.norm() / v.norm().max(1.0);
            let gap = |t: f64| (r + v * t).norm();
            let numeric = (gap(h) - gap(-h)) / (2.0 * h);
            let rate = a.range_rate_mps(b);
            assert!(close(rate, numeric, 1e-6 * v.norm()), "{rate} vs {numeric}");
            assert_eq!(rate.to_bits(), (-(-rate)).to_bits());
            // It's the same seen from either end.
            assert!(close(b.range_rate_mps(a), rate, 1e-12 * v.norm()));
        }
    }

    #[test]
    fn range_rate_signs_and_edges() {
        let here = Track::new(Ecef::new(1.0, 2.0, 3.0), Vec3::ZERO);
        let at = |x: f64, v: Vec3| Track::new(Ecef::new(1.0 + x, 2.0, 3.0), v);
        assert_eq!(here.range_rate_mps(at(100.0, Vec3::new(-5.0, 0.0, 0.0))), -5.0); // closing
        assert_eq!(here.range_rate_mps(at(100.0, Vec3::new(5.0, 0.0, 0.0))), 5.0); // opening
        assert_eq!(here.range_rate_mps(at(100.0, Vec3::new(0.0, 7.0, 0.0))), 0.0); // across
                                                                                   // At zero range they can only move apart.
        assert_eq!(here.range_rate_mps(at(0.0, Vec3::new(0.0, 3.0, 4.0))), 5.0);
        assert_eq!(here.range_rate_mps(here), 0.0);
        assert!(here.range_rate_mps(at(f64::NAN, Vec3::ZERO)).is_nan());
        assert!(here
            .range_rate_mps(at(1.0, Vec3::new(f64::INFINITY, 0.0, 0.0)))
            .is_infinite());
    }

    /// Build each case backwards from where and when they meet, so the answer is known, then
    /// check the solver finds it from the start.
    #[test]
    fn intercepts_built_backwards_from_the_meeting_point() {
        let mut seed = 42;
        for i in 0..20_000 {
            let meet = Ecef::from(unit(&mut seed) * 6.4e6);
            let when = 10f64.powf(-1.0 + 5.0 * rnd(&mut seed));
            let v = unit(&mut seed) * 10f64.powf(4.0 * rnd(&mut seed));
            // Pursuer speeds from a fifth of the target's to five times it, and exactly equal.
            let ratio = if i % 5 == 0 {
                1.0
            } else {
                10f64.powf(1.4 * rnd(&mut seed) - 0.7)
            };
            let s = v.norm() * ratio;
            let heading = unit(&mut seed);
            let target = Track::new(meet.offset(-v * when), v);
            let from = meet.offset(heading * (-s * when));
            let hit = target.intercept_from(from, s).expect("an intercept exists");
            // The gap left at the time it gives: how far the pursuer would miss by.
            let d = target.position.vec() - from.vec();
            let gap = |t: f64| (d + v * t).norm() - s * t;
            // Positions near 6 400 km are rounded to about a nanometre.
            let slack = 1e-8 + 1e-12 * s * when;
            assert!(gap(hit.time_s).abs() <= slack, "case {i}: miss {}", gap(hit.time_s));
            assert!(close(hit.speed_mps(), s, 1e-12 * s + 1e-9 / hit.time_s));
            assert!(close(hit.position.distance_to(target.at(hit.time_s)), 0.0, 1e-9));
            // Nothing sooner: the gap stays open until then.
            for k in 0..10 {
                assert!(gap(hit.time_s * k as f64 / 10.0) > -slack, "case {i}: an earlier root");
            }
            // The rounding moves the meeting time by about slack / (how fast the gap closes there).
            let shift = slack / (s - heading.dot(v)).abs() + 1e-12 * when;
            if ratio > 1.2 {
                // A faster pursuer has only one way to meet the target: the time we built.
                assert!(close(hit.time_s, when, shift), "case {i}: {} vs {when}", hit.time_s);
            } else {
                // A slower one may get there sooner, but never later.
                assert!(hit.time_s <= when + shift, "case {i}: {} vs {when}", hit.time_s);
            }
        }
    }

    #[test]
    fn intercept_at_the_same_speed_and_the_edge_cases() {
        let from = Ecef::new(0.0, 0.0, 0.0);
        // Coming straight at us at our own speed: we meet half way. (The textbook formula divides
        // by zero here.)
        let head_on = Track::new(Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(-50.0, 0.0, 0.0));
        let hit = head_on.intercept_from(from, 50.0).unwrap();
        assert_eq!(hit.time_s, 10.0);
        assert_eq!(hit.position, Ecef::new(500.0, 0.0, 0.0));
        assert_eq!(hit.velocity, Vec3::new(50.0, 0.0, 0.0));
        // Crossing at our own speed, or moving away no slower: never.
        let across = Track::new(Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(0.0, 50.0, 0.0));
        assert_eq!(across.intercept_from(from, 50.0), None);
        let away = Track::new(Ecef::new(1_000.0, 0.0, 0.0), Vec3::new(30.0, 40.0, 0.0));
        assert_eq!(away.intercept_from(from, 50.0), None);
        assert_eq!(away.intercept_from(from, 49.0), None);
        assert!(away.intercept_from(from, 51.0).is_some());
        // A still target: straight there.
        let still = Track::new(Ecef::new(0.0, 300.0, 400.0), Vec3::ZERO);
        assert_eq!(still.intercept_from(from, 10.0).unwrap().time_s, 50.0);
        assert_eq!(still.intercept_from(from, 0.0), None);
        // Already there.
        let there = Track::new(from, Vec3::new(1.0, 2.0, 3.0));
        let now = there.intercept_from(from, 0.0).unwrap();
        assert_eq!((now.time_s, now.velocity), (0.0, Vec3::ZERO));
        assert_eq!(there.slowest_intercept_from(from).unwrap().time_s, 0.0);
        // A target that passes straight through us can be met standing still.
        let through = Track::new(Ecef::new(0.0, 0.0, 600.0), Vec3::new(0.0, 0.0, -20.0));
        assert_eq!(through.intercept_from(from, 0.0).unwrap().time_s, 30.0);
        // Bad speeds.
        for s in [-1.0, f64::NAN] {
            assert_eq!(head_on.intercept_from(from, s), None);
        }
        assert!(head_on
            .intercept_from(from, f64::INFINITY)
            .is_none_or(|h| h.time_s == 0.0));
        let nan = Track::new(Ecef::new(f64::NAN, 0.0, 0.0), Vec3::ZERO);
        assert_eq!(nan.intercept_from(from, 10.0), None);
        assert_eq!(nan.slowest_intercept_from(from), None);
    }

    #[test]
    fn the_slowest_intercept_is_the_slowest() {
        let mut seed = 7;
        for _ in 0..20_000 {
            let from = Ecef::from(unit(&mut seed) * 6.4e6);
            let target = Track::new(
                from.offset(unit(&mut seed) * 10f64.powf(1.0 + 5.0 * rnd(&mut seed))),
                unit(&mut seed) * 10f64.powf(4.0 * rnd(&mut seed)),
            );
            let d = target.position.vec() - from.vec();
            let Some(slowest) = target.slowest_intercept_from(from) else {
                // Not coming towards us: only speeds above its own work.
                assert!(d.dot(target.velocity) >= 0.0);
                let v = target.velocity.norm();
                assert!(target.intercept_from(from, v * 0.999).is_none());
                assert!(target.intercept_from(from, v * 1.001).is_some());
                continue;
            };
            let s = slowest.speed_mps();
            // It's the part of the target's velocity across the line between them.
            let across = d.cross(target.velocity).norm() / d.norm();
            assert!(close(s, across, 1e-12 * target.velocity.norm()), "{s} vs {across}");
            let miss = ((d + target.velocity * slowest.time_s).norm() - s * slowest.time_s).abs();
            assert!(miss <= 1e-8 + 1e-11 * s * slowest.time_s, "miss {miss}");
            // A little slower can't make it; a little faster can, and sooner.
            assert!(target.intercept_from(from, s * (1.0 - 1e-6)).is_none());
            let faster = target.intercept_from(from, s * (1.0 + 1e-6)).unwrap();
            assert!(faster.time_s <= slowest.time_s * (1.0 + 1e-12));
        }
    }

    #[test]
    fn course_conventions() {
        // Hovering: course 0 by convention.
        assert_eq!(
            Course::from_enu_velocity(Enu::new(0.0, 0.0, 2.0)),
            Course::new(0.0, 0.0, 2.0)
        );
        assert!(close(Course::new(270.0, 3.0, 4.0).speed_mps(), 5.0, 1e-12));
        let w = Course::new(270.0, 10.0, 0.0).to_enu_velocity();
        assert!(close(w.east, -10.0, 1e-9) && close(w.north, 0.0, 1e-9));
        let c = Course::from_enu_velocity(Course::new(123.0, 50.0, -3.0).to_enu_velocity());
        assert!(close(c.course_deg, 123.0, 1e-9) && close(c.ground_speed_mps, 50.0, 1e-9) && c.climb_mps == -3.0);
    }
}
