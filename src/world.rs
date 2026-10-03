use std::fmt;

use crate::error::{invalid, Result, SimError};
use crate::events::{Event, EventHit};
use crate::forces::{Force, ForceId, ForceSet, Param};
use crate::integrators::{self, Integrator};
use crate::state::State;
use crate::vec3::Vec3;

/// A simulation: particle state + forces + time integrator.
pub struct World {
    pub state: State,
    pub forces: ForceSet,
    integrator: Box<dyn Integrator>,
    /// State before the current step, restored if the step fails.
    backup: State,
}

impl World {
    pub fn new(integrator: Box<dyn Integrator>) -> Self {
        Self {
            state: State::new(),
            forces: ForceSet::new(),
            integrator,
            backup: State::new(),
        }
    }

    /// Creates a world using an integrator from [`integrators::NAMES`].
    pub fn with_integrator(name: &str) -> Result<Self> {
        Ok(Self::new(integrators::by_name(name)?))
    }

    pub fn integrator(&self) -> &dyn Integrator {
        self.integrator.as_ref()
    }

    pub fn set_integrator(&mut self, integrator: Box<dyn Integrator>) {
        self.integrator = integrator;
    }

    /// Adds a particle and returns its index. A mass of zero makes a test particle:
    /// it feels gravity but does not source it.
    pub fn add_particle(&mut self, pos: Vec3, vel: Vec3, mass: f64) -> Result<usize> {
        check_mass(mass)?;
        Ok(self.state.add_particle(pos, vel, mass))
    }

    /// Removes particle `i`; particles above it shift down by one index, and forces that
    /// refer to particles by index are renumbered. Fails if a force refers to particle `i`.
    pub fn remove_particle(&mut self, i: usize) -> Result<()> {
        self.check_particle(i)?;
        let users = self.forces.referencing(i);
        if !users.is_empty() {
            return invalid(format!(
                "cannot remove particle {i}: still used by {}; remove those forces first",
                users.join(", ")
            ));
        }
        self.forces.particle_removed(i);
        self.state.remove_particle(i);
        Ok(())
    }

    pub fn set_masses(&mut self, masses: Vec<f64>) -> Result<()> {
        if masses.len() != self.state.len() {
            return invalid(format!(
                "expected {} masses, got {}",
                self.state.len(),
                masses.len()
            ));
        }
        masses.iter().try_for_each(|&m| check_mass(m))?;
        self.state.mass = masses;
        Ok(())
    }

    /// Sets velocities; pinned particles must be given zero velocity.
    pub fn set_velocities(&mut self, vel: Vec<Vec3>) -> Result<()> {
        if vel.len() != self.state.len() {
            return invalid(format!(
                "expected {} velocities, got {}",
                self.state.len(),
                vel.len()
            ));
        }
        if let Some(i) = (0..vel.len()).find(|&i| self.state.pinned[i] && vel[i] != Vec3::ZERO) {
            return invalid(format!("particle {i} is pinned; its velocity must be zero"));
        }
        self.state.vel = vel;
        Ok(())
    }

    /// Pins particle `i` in place (velocity set to zero) or releases it.
    pub fn pin(&mut self, i: usize, pinned: bool) -> Result<()> {
        self.check_particle(i)?;
        self.state.pinned[i] = pinned;
        if pinned {
            self.state.vel[i] = Vec3::ZERO;
        }
        Ok(())
    }

    fn check_particle(&self, i: usize) -> Result<()> {
        if i >= self.state.len() {
            return invalid(format!(
                "particle index {i} out of range (have {} particles)",
                self.state.len()
            ));
        }
        Ok(())
    }

    pub fn add_force(&mut self, force: impl Force + 'static) -> ForceId {
        self.forces.add(Box::new(force))
    }

    pub fn remove_force(&mut self, id: ForceId) -> Result<Box<dyn Force>> {
        self.forces.remove(id)
    }

    /// Sets parameters of one force, e.g. `[("G", Param::Scalar(2.0))]`. All-or-nothing.
    pub fn set_force_params(&mut self, id: ForceId, values: &[(String, Param)]) -> Result<()> {
        self.forces.set_params(id, values)
    }

    /// Advances one step. On failure the state is left as it was before the step.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        if !(dt.is_finite() && dt != 0.0) {
            return invalid(format!("dt must be finite and non-zero, got {dt}"));
        }
        // Keep pinned particles exactly at rest even if `state.vel` was edited directly.
        for (v, &p) in self.state.vel.iter_mut().zip(&self.state.pinned) {
            if p {
                *v = Vec3::ZERO;
            }
        }
        self.backup.clone_from(&self.state);
        let result = self.integrator.step(&mut self.state, &self.forces, dt).and_then(|()| {
            let finite = self.state.pos.iter().chain(&self.state.vel).all(|v| {
                v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
            });
            if finite {
                Ok(())
            } else {
                invalid(format!(
                    "state became non-finite while stepping from t = {} (time step too large or singular force?)",
                    self.backup.t
                ))
            }
        });
        if result.is_err() {
            std::mem::swap(&mut self.state, &mut self.backup);
        }
        result
    }

    /// Takes `steps` steps of size `dt`, recording the initial state, every
    /// `record_every`-th step, and the final state.
    ///
    /// If a step fails, the world is left at the last completed step and the error
    /// comes back together with the frames recorded so far, ending at that step.
    pub fn run(
        &mut self,
        dt: f64,
        steps: usize,
        record_every: usize,
    ) -> std::result::Result<Trajectory, RunFailure> {
        self.run_with_events(dt, steps, record_every, &[])
    }

    /// Like [`World::run`], also detecting `events`. Every crossing is located to
    /// ~1e-12 of a step and stored in [`Trajectory::events`] with the state at that time.
    /// A terminal event stops the run there: the world is left at the event state, which
    /// is recorded as the final frame, and [`Trajectory::terminated_by`] names the event.
    /// Two crossings of the same event within one step are not detected (use a smaller `dt`).
    pub fn run_with_events(
        &mut self,
        dt: f64,
        steps: usize,
        record_every: usize,
        events: &[Event],
    ) -> std::result::Result<Trajectory, RunFailure> {
        let mut traj = Trajectory::new(self.state.len());
        let fail = |error, trajectory| RunFailure {
            error,
            trajectory: Box::new(trajectory),
        };
        if record_every == 0 {
            return Err(fail(
                SimError::Invalid("record_every must be at least 1".into()),
                traj,
            ));
        }
        if let Err(e) = traj.record(self) {
            return Err(fail(e, traj));
        }
        let mut g = match self.event_values(events) {
            Ok(g) => g,
            Err(e) => return Err(fail(e, traj)),
        };
        let mut start = State::new();
        for k in 1..=steps {
            if !events.is_empty() {
                start.clone_from(&self.state);
            }
            let stepped = self.step(dt).and_then(|()| {
                if events.is_empty() {
                    return Ok(false);
                }
                let g_new = self.event_values(events)?;
                let terminated =
                    self.handle_crossings(events, &start, &g, &g_new, dt, &mut traj)?;
                g = g_new;
                Ok(terminated)
            });
            match stepped {
                Ok(true) => {
                    if traj.t.last() != Some(&self.state.t) {
                        if let Err(e) = traj.record(self) {
                            return Err(fail(e, traj));
                        }
                    }
                    return Ok(traj);
                }
                Ok(false) => {}
                Err(e) => {
                    if !events.is_empty() && self.state.t != start.t {
                        // The step itself succeeded but event handling failed: undo it too.
                        self.state.clone_from(&start);
                    }
                    if traj.t.last() != Some(&self.state.t) {
                        let _ = traj.record(self);
                    }
                    return Err(fail(e, traj));
                }
            }
            if k % record_every == 0 || k == steps {
                if let Err(e) = traj.record(self) {
                    return Err(fail(e, traj));
                }
            }
        }
        Ok(traj)
    }

    fn event_values(&self, events: &[Event]) -> Result<Vec<f64>> {
        let s = &self.state;
        events
            .iter()
            .map(|e| e.function.value(s.t, &s.pos, &s.vel, &s.mass))
            .collect()
    }

    /// Locates the events that fired during the step from `start` (values `g0`) to the
    /// current state (values `g1`) and records them in time order. Returns `true` if a
    /// terminal event fired, in which case the world is left at that event.
    fn handle_crossings(
        &mut self,
        events: &[Event],
        start: &State,
        g0: &[f64],
        g1: &[f64],
        dt: f64,
        traj: &mut Trajectory,
    ) -> Result<bool> {
        let fired: Vec<usize> = (0..events.len())
            .filter(|&e| events[e].triggers(g0[e], g1[e]))
            .collect();
        if fired.is_empty() {
            return Ok(false);
        }
        let mut hits = Vec::with_capacity(fired.len());
        for e in fired {
            let (theta, state) = self.locate(&events[e], start, g0[e], g1[e], dt)?;
            hits.push((theta, e, state));
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let stop = hits.iter().position(|(_, e, _)| events[*e].terminal);
        let keep = stop.map_or(hits.len(), |s| s + 1);
        for (_, e, state) in hits.drain(..).take(keep) {
            traj.events.push(EventHit {
                event: e,
                t: state.t,
                pos: state.pos.clone(),
                vel: state.vel.clone(),
            });
            if events[e].terminal {
                traj.terminated_by = Some(e);
                self.state = state;
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Finds the fraction `theta` of the step at which `event` crosses zero, by stepping
    /// from `start` with `theta * dt` (Illinois variant of regula falsi). Returns the end of
    /// the final bracket that lies past the crossing, and the state there.
    fn locate(
        &mut self,
        event: &Event,
        start: &State,
        g0: f64,
        g1: f64,
        dt: f64,
    ) -> Result<(f64, State)> {
        const TOLERANCE: f64 = 1e-12;
        let (mut lo, mut f_lo, mut hi, mut f_hi) = (0.0, g0, 1.0, g1);
        let mut last_side = 0;
        for _ in 0..100 {
            if hi - lo <= TOLERANCE || f_hi == 0.0 {
                break;
            }
            let mut theta = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
            if !(theta > lo && theta < hi) {
                theta = 0.5 * (lo + hi);
                if !(theta > lo && theta < hi) {
                    break; // bracket cannot shrink further in floating point
                }
            }
            let s = self.advance(start, theta * dt)?;
            let f = event.function.value(s.t, &s.pos, &s.vel, &s.mass)?;
            if f * f_lo > 0.0 {
                (lo, f_lo) = (theta, f);
                if last_side == -1 {
                    f_hi *= 0.5;
                }
                last_side = -1;
            } else {
                (hi, f_hi) = (theta, f);
                if last_side == 1 {
                    f_lo *= 0.5;
                }
                last_side = 1;
            }
        }
        Ok((hi, self.advance(start, hi * dt)?))
    }

    /// The integrator's solution after a step of `h` from `start` (the world is unchanged).
    fn advance(&mut self, start: &State, h: f64) -> Result<State> {
        let mut s = start.clone();
        self.integrator.step(&mut s, &self.forces, h)?;
        Ok(s)
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.state.kinetic_energy()
    }

    /// Potential energy of all conservative forces.
    pub fn potential_energy(&self) -> Result<f64> {
        self.forces
            .potential(self.state.t, &self.state.pos, &self.state.mass)
    }

    pub fn total_energy(&self) -> Result<f64> {
        Ok(self.kinetic_energy() + self.potential_energy()?)
    }

    /// Current total acceleration of every particle.
    pub fn accelerations(&self) -> Result<Vec<Vec3>> {
        let s = &self.state;
        let mut acc = Vec::new();
        self.forces
            .accelerations(s.t, &s.pos, &s.vel, &s.mass, &s.pinned, &mut acc)?;
        Ok(acc)
    }
}

fn check_mass(mass: f64) -> Result<()> {
    if !(mass.is_finite() && mass >= 0.0) {
        return invalid(format!("mass must be non-negative and finite, got {mass}"));
    }
    Ok(())
}

/// A failed [`World::run`]: the error plus the frames recorded before it.
#[derive(Debug)]
pub struct RunFailure {
    pub error: SimError,
    pub trajectory: Box<Trajectory>,
}

impl fmt::Display for RunFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} (after {} recorded frames)",
            self.error,
            self.trajectory.n_frames()
        )
    }
}

impl std::error::Error for RunFailure {}

impl From<RunFailure> for SimError {
    fn from(f: RunFailure) -> SimError {
        f.error
    }
}

/// Recorded frames of a run. Positions and velocities are stored frame-major:
/// particle `p` of frame `f` is at index `f * n_particles + p`.
#[derive(Debug, Clone, Default)]
pub struct Trajectory {
    pub n_particles: usize,
    pub t: Vec<f64>,
    pub pos: Vec<Vec3>,
    pub vel: Vec<Vec3>,
    pub kinetic: Vec<f64>,
    pub potential: Vec<f64>,
    /// Events detected during the run, in time order.
    pub events: Vec<EventHit>,
    /// The terminal event that stopped the run, if any.
    pub terminated_by: Option<usize>,
}

impl Trajectory {
    pub fn new(n_particles: usize) -> Self {
        Self {
            n_particles,
            ..Default::default()
        }
    }

    pub fn n_frames(&self) -> usize {
        self.t.len()
    }

    /// Appends the current state. Nothing is appended if computing the energies fails.
    fn record(&mut self, world: &World) -> Result<()> {
        let potential = world.potential_energy()?;
        self.t.push(world.state.t);
        self.pos.extend_from_slice(&world.state.pos);
        self.vel.extend_from_slice(&world.state.vel);
        self.kinetic.push(world.kinetic_energy());
        self.potential.push(potential);
        Ok(())
    }

    pub fn total_energy(&self) -> Vec<f64> {
        self.kinetic
            .iter()
            .zip(&self.potential)
            .map(|(k, u)| k + u)
            .collect()
    }

    /// Positions of frame `f`.
    pub fn frame_pos(&self, f: usize) -> &[Vec3] {
        &self.pos[f * self.n_particles..(f + 1) * self.n_particles]
    }
}
