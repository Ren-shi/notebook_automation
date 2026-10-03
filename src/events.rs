//! Events: moments when a scalar function of the state crosses zero, such as a
//! periapsis passage (radial velocity changes sign) or a particle reaching the ground.
//!
//! [`World::run_with_events`](crate::World::run_with_events) checks every event after each
//! step. When one changes sign, the crossing is located by re-stepping the integrator from
//! the start of the step with a fraction of `dt` (Illinois root finding), so the event
//! state is the integrator's own solution at the event time, as accurate as any step.

use crate::error::{invalid, Result};
use crate::vec3::Vec3;

/// A scalar function of the state whose zeros are events.
pub trait EventFunction: Send + Sync {
    fn value(&self, t: f64, pos: &[Vec3], vel: &[Vec3], mass: &[f64]) -> Result<f64>;
    fn name(&self) -> String;
}

/// Which sign changes count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Negative to positive (or reaching zero from below).
    Rising,
    /// Positive to negative (or reaching zero from above).
    Falling,
    Either,
}

impl Direction {
    /// `+1` rising, `-1` falling, `0` either.
    pub fn from_sign(sign: i32) -> Result<Direction> {
        match sign {
            1 => Ok(Direction::Rising),
            -1 => Ok(Direction::Falling),
            0 => Ok(Direction::Either),
            _ => invalid(format!("direction must be -1, 0 or +1, got {sign}")),
        }
    }
}

pub struct Event {
    pub function: Box<dyn EventFunction>,
    pub direction: Direction,
    /// Stop the run at this event.
    pub terminal: bool,
}

impl Event {
    pub fn new(
        function: impl EventFunction + 'static,
        direction: Direction,
        terminal: bool,
    ) -> Self {
        Self {
            function: Box::new(function),
            direction,
            terminal,
        }
    }

    /// Whether going from `g0` to `g1` over a step is a crossing of this event.
    /// Starting exactly at zero is not a crossing, so an event is never reported twice.
    pub fn triggers(&self, g0: f64, g1: f64) -> bool {
        let rising = g0 < 0.0 && g1 >= 0.0;
        let falling = g0 > 0.0 && g1 <= 0.0;
        match self.direction {
            Direction::Rising => rising,
            Direction::Falling => falling,
            Direction::Either => rising || falling,
        }
    }
}

/// One detected event: which event fired, when, and the state at that moment.
#[derive(Debug, Clone)]
pub struct EventHit {
    /// Index into the events passed to the run.
    pub event: usize,
    pub t: f64,
    pub pos: Vec<Vec3>,
    pub vel: Vec<Vec3>,
}

fn check_index(i: usize, n: usize, what: &str) -> Result<()> {
    if i >= n {
        return invalid(format!(
            "{what}: particle index {i} out of range (have {n} particles)"
        ));
    }
    Ok(())
}

/// `(r_i - r_j) · (v_i - v_j)`, proportional to the rate of change of the separation.
/// Rising crossings are periapses (closest approach), falling ones apoapses.
/// With `j = None` the reference point is the origin, at rest.
#[derive(Debug, Clone)]
pub struct RadialVelocity {
    pub i: usize,
    pub j: Option<usize>,
}

impl EventFunction for RadialVelocity {
    fn value(&self, _t: f64, pos: &[Vec3], vel: &[Vec3], _mass: &[f64]) -> Result<f64> {
        check_index(self.i, pos.len(), "RadialVelocity")?;
        let (r, v) = match self.j {
            Some(j) => {
                check_index(j, pos.len(), "RadialVelocity")?;
                (pos[self.i] - pos[j], vel[self.i] - vel[j])
            }
            None => (pos[self.i], vel[self.i]),
        };
        Ok(r.dot(v))
    }

    fn name(&self) -> String {
        match self.j {
            Some(j) => format!("RadialVelocity({}, {j})", self.i),
            None => format!("RadialVelocity({})", self.i),
        }
    }
}

/// `pos_i[axis] - value`: a particle crossing a plane (Poincaré sections, hitting the ground).
#[derive(Debug, Clone)]
pub struct CoordinateCrossing {
    pub i: usize,
    pub axis: usize,
    pub value: f64,
}

impl EventFunction for CoordinateCrossing {
    fn value(&self, _t: f64, pos: &[Vec3], _vel: &[Vec3], _mass: &[f64]) -> Result<f64> {
        check_index(self.i, pos.len(), "CoordinateCrossing")?;
        if self.axis > 2 {
            return invalid(format!(
                "CoordinateCrossing: axis must be 0, 1 or 2, got {}",
                self.axis
            ));
        }
        Ok(pos[self.i][self.axis] - self.value)
    }

    fn name(&self) -> String {
        format!(
            "CoordinateCrossing({}, {}={})",
            self.i,
            ["x", "y", "z"][self.axis.min(2)],
            self.value
        )
    }
}

/// `|r_i - r_j| - distance`: two particles reaching a given separation
/// (falling: approach to contact; rising: escape).
#[derive(Debug, Clone)]
pub struct Separation {
    pub i: usize,
    pub j: usize,
    pub distance: f64,
}

impl EventFunction for Separation {
    fn value(&self, _t: f64, pos: &[Vec3], _vel: &[Vec3], _mass: &[f64]) -> Result<f64> {
        check_index(self.i, pos.len(), "Separation")?;
        check_index(self.j, pos.len(), "Separation")?;
        Ok((pos[self.i] - pos[self.j]).norm() - self.distance)
    }

    fn name(&self) -> String {
        format!("Separation({}, {})", self.i, self.j)
    }
}
