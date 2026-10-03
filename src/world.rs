use std::fmt;

use crate::constraints::{self, Anchor, ConstraintId, Constraints, Rod};
use crate::error::{invalid, Result, SimError};
use crate::events::{Event, EventHit};
use crate::forces::{Force, ForceId, ForceSet, Param};
use crate::integrators::{self, Integrator};
use crate::state::State;
use crate::vec3::Vec3;

/// A simulation: particle state + forces + constraints + time integrator.
pub struct World {
    pub state: State,
    pub forces: ForceSet,
    /// Rigid rods, enforced with RATTLE (needs the `verlet` or `yoshida4` integrator).
    pub constraints: Constraints,
    integrator: Box<dyn Integrator>,
    /// Tension in each constraint at the end of the last step.
    pub(crate) tension: Vec<f64>,
    /// State before the current step, restored if the step fails.
    backup: State,
    backup_tension: Vec<f64>,
}

impl World {
    pub fn new(integrator: Box<dyn Integrator>) -> Self {
        Self {
            state: State::new(),
            forces: ForceSet::new(),
            constraints: Constraints::new(),
            integrator,
            tension: Vec::new(),
            backup: State::new(),
            backup_tension: Vec::new(),
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
        let mut users = self.forces.referencing(i);
        users.extend(self.constraints.referencing(i));
        if !users.is_empty() {
            return invalid(format!(
                "cannot remove particle {i}: still used by {}; remove those first",
                users.join(", ")
            ));
        }
        self.forces.particle_removed(i);
        self.constraints.particle_removed(i);
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
        let old = std::mem::replace(&mut self.state.mass, masses);
        if let Err(e) = self.constraints.validate(&self.state) {
            self.state.mass = old;
            return Err(e);
        }
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

    /// Adds a rigid rod keeping particle `i` at a fixed distance from `anchor` (another
    /// particle or a fixed point) and returns its id. `length` defaults to the current
    /// distance; if given, it must match the current distance to the constraint tolerance.
    ///
    /// Velocity components along the rod are removed in the first step.
    pub fn add_rod(
        &mut self,
        i: usize,
        anchor: Anchor,
        length: Option<f64>,
    ) -> Result<ConstraintId> {
        let mut rod = Rod {
            i,
            anchor,
            length: 1.0,
        };
        // Index and mass checks first, so `separation` cannot index out of range.
        constraints::check_rod(&rod, &self.state)?;
        let current = rod.separation(&self.state.pos).norm();
        rod.length = length.unwrap_or(current);
        constraints::check_rod(&rod, &self.state)?;
        if (current / rod.length - 1.0).abs() > self.constraints.tolerance.max(1e-12) {
            return invalid(format!(
                "{}: particles are {current} apart but the length is {}; place them first",
                rod.name(),
                rod.length
            ));
        }
        let id = self.constraints.add(rod);
        self.tension = vec![0.0; self.constraints.len()];
        Ok(id)
    }

    pub fn remove_constraint(&mut self, id: ConstraintId) -> Result<Rod> {
        let rod = self.constraints.remove(id)?;
        self.tension = vec![0.0; self.constraints.len()];
        Ok(rod)
    }

    pub fn clear_constraints(&mut self) {
        self.constraints.clear();
        self.tension.clear();
    }

    /// Tension in each constraint (in [`Constraints::iter`] order) at the end of the last
    /// step: the force pulling the rod's ends together, negative when the rod pushes them
    /// apart. Zero before the first step.
    pub fn constraint_tensions(&self) -> &[f64] {
        &self.tension
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
        self.backup_tension.clone_from(&self.tension);
        let result = Self::integrate(
            self.integrator.as_mut(),
            &self.forces,
            &self.constraints,
            &mut self.state,
            &mut self.tension,
            dt,
        )
        .and_then(|()| {
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
            std::mem::swap(&mut self.tension, &mut self.backup_tension);
        }
        result
    }

    /// One integrator step, constrained when there are constraints.
    fn integrate(
        integrator: &mut dyn Integrator,
        forces: &ForceSet,
        constraints: &Constraints,
        state: &mut State,
        tension: &mut Vec<f64>,
        dt: f64,
    ) -> Result<()> {
        if constraints.is_empty() {
            tension.clear();
            integrator.step(state, forces, dt)
        } else {
            constraints.validate(state)?;
            integrator.step_constrained(state, forces, constraints, tension, dt)
        }
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
        self.run_with_options(
            dt,
            steps,
            &RunOptions {
                record_every,
                events,
                energies: true,
            },
        )
    }

    /// [`World::run_with_events`] with all options, collecting frames in memory.
    pub fn run_with_options(
        &mut self,
        dt: f64,
        steps: usize,
        options: &RunOptions<'_>,
    ) -> std::result::Result<Trajectory, RunFailure> {
        let mut traj = Trajectory::new(self.state.len(), options.energies);
        if let Some(frames) = steps.checked_div(options.record_every) {
            traj.reserve(frames + 2);
        }
        match self.run_into(dt, steps, options, &mut traj) {
            Ok(terminated_by) => {
                traj.terminated_by = terminated_by;
                Ok(traj)
            }
            Err(error) => Err(RunFailure {
                error,
                trajectory: Box::new(traj),
            }),
        }
    }

    /// The run loop behind [`World::run`]: hands every recorded frame and event to
    /// `recorder` instead of keeping them, so memory stays flat however long the run.
    /// Returns the index of the terminal event that stopped the run, if any.
    ///
    /// On error the world is left at the last completed step, and that step has been
    /// recorded (unless recording itself failed).
    pub fn run_into(
        &mut self,
        dt: f64,
        steps: usize,
        options: &RunOptions<'_>,
        recorder: &mut dyn Recorder,
    ) -> Result<Option<usize>> {
        let RunOptions {
            record_every,
            events,
            energies,
        } = *options;
        if record_every == 0 {
            return invalid("record_every must be at least 1");
        }
        self.emit(recorder, energies)?;
        let mut last_recorded = self.state.t;
        let mut g = self.event_values(events)?;
        let mut start = State::new();
        let mut start_tension = Vec::new();
        for k in 1..=steps {
            if !events.is_empty() {
                start.clone_from(&self.state);
                start_tension.clone_from(&self.tension);
            }
            let stepped = self.step(dt).and_then(|()| {
                if events.is_empty() {
                    return Ok(None);
                }
                let g_new = self.event_values(events)?;
                let terminated = self.handle_crossings(events, &start, &g, &g_new, dt, recorder)?;
                g = g_new;
                Ok(terminated)
            });
            match stepped {
                Ok(Some(e)) => {
                    if last_recorded != self.state.t {
                        self.emit(recorder, energies)?;
                    }
                    return Ok(Some(e));
                }
                Ok(None) => {}
                Err(e) => {
                    if !events.is_empty() && self.state.t != start.t {
                        // The step itself succeeded but event handling failed: undo it too.
                        self.state.clone_from(&start);
                        self.tension.clone_from(&start_tension);
                    }
                    if last_recorded != self.state.t {
                        let _ = self.emit(recorder, energies);
                    }
                    return Err(e);
                }
            }
            if k % record_every == 0 || k == steps {
                self.emit(recorder, energies)?;
                last_recorded = self.state.t;
            }
        }
        Ok(None)
    }

    /// Passes the current state to `recorder`. Nothing is passed if computing the energies fails.
    fn emit(&self, recorder: &mut dyn Recorder, energies: bool) -> Result<()> {
        let (kinetic, potential) = if energies {
            (Some(self.kinetic_energy()), Some(self.potential_energy()?))
        } else {
            (None, None)
        };
        recorder.frame(&Frame {
            t: self.state.t,
            pos: &self.state.pos,
            vel: &self.state.vel,
            kinetic,
            potential,
            tension: &self.tension,
        })
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
    /// terminal event fired, its index; the world is then left at that event.
    fn handle_crossings(
        &mut self,
        events: &[Event],
        start: &State,
        g0: &[f64],
        g1: &[f64],
        dt: f64,
        recorder: &mut dyn Recorder,
    ) -> Result<Option<usize>> {
        let fired: Vec<usize> = (0..events.len())
            .filter(|&e| events[e].triggers(g0[e], g1[e]))
            .collect();
        if fired.is_empty() {
            return Ok(None);
        }
        let mut hits = Vec::with_capacity(fired.len());
        for e in fired {
            let (theta, state, tension) = self.locate(&events[e], start, g0[e], g1[e], dt)?;
            hits.push((theta, e, state, tension));
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let stop = hits.iter().position(|(_, e, _, _)| events[*e].terminal);
        let keep = stop.map_or(hits.len(), |s| s + 1);
        for (_, e, state, tension) in hits.drain(..).take(keep) {
            recorder.event(EventHit {
                event: e,
                t: state.t,
                pos: state.pos.clone(),
                vel: state.vel.clone(),
            })?;
            if events[e].terminal {
                self.state = state;
                self.tension = tension;
                return Ok(Some(e));
            }
        }
        Ok(None)
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
    ) -> Result<(f64, State, Vec<f64>)> {
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
            let (s, _) = self.advance(start, theta * dt)?;
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
        let (state, tension) = self.advance(start, hi * dt)?;
        Ok((hi, state, tension))
    }

    /// The integrator's solution after a step of `h` from `start`, and the constraint
    /// tensions there (the world is unchanged).
    fn advance(&mut self, start: &State, h: f64) -> Result<(State, Vec<f64>)> {
        let mut s = start.clone();
        let mut tension = Vec::new();
        Self::integrate(
            self.integrator.as_mut(),
            &self.forces,
            &self.constraints,
            &mut s,
            &mut tension,
            h,
        )?;
        Ok((s, tension))
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

/// Options for [`World::run_with_options`] and [`World::run_into`].
#[derive(Clone, Copy)]
pub struct RunOptions<'a> {
    /// Record every `record_every`-th step (the initial and final states are always recorded).
    pub record_every: usize,
    /// Events to detect; see [`World::run_with_events`].
    pub events: &'a [Event],
    /// Compute kinetic and potential energy for each frame. The potential is a full pass over
    /// the forces (O(N²) for gravity), so turning this off makes frequent recording cheap.
    pub energies: bool,
}

impl Default for RunOptions<'_> {
    fn default() -> Self {
        Self {
            record_every: 1,
            events: &[],
            energies: true,
        }
    }
}

/// One recorded frame. The slices borrow the world's state, so copy what you keep.
#[derive(Debug, Clone, Copy)]
pub struct Frame<'a> {
    pub t: f64,
    pub pos: &'a [Vec3],
    pub vel: &'a [Vec3],
    /// `None` when the run does not compute energies.
    pub kinetic: Option<f64>,
    pub potential: Option<f64>,
    /// Tension in each constraint (see [`World::constraint_tensions`]).
    pub tension: &'a [f64],
}

/// Receives the output of [`World::run_into`] as it is produced: write it to disk, reduce it
/// on the fly, or collect it (as [`Trajectory`] does). An error stops the run.
pub trait Recorder {
    fn frame(&mut self, frame: &Frame<'_>) -> Result<()>;

    /// An event detected during the run, in time order. Ignored by default.
    fn event(&mut self, _hit: EventHit) -> Result<()> {
        Ok(())
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
    /// Whether energies were recorded; if not, `kinetic` and `potential` stay empty.
    pub energies: bool,
    pub kinetic: Vec<f64>,
    pub potential: Vec<f64>,
    /// Constraint tensions, frame-major: `n_constraints` values per frame.
    pub n_constraints: usize,
    pub tension: Vec<f64>,
    /// Events detected during the run, in time order.
    pub events: Vec<EventHit>,
    /// The terminal event that stopped the run, if any.
    pub terminated_by: Option<usize>,
}

impl Trajectory {
    pub fn new(n_particles: usize, energies: bool) -> Self {
        Self {
            n_particles,
            energies,
            ..Default::default()
        }
    }

    pub fn n_frames(&self) -> usize {
        self.t.len()
    }

    /// Kinetic plus potential energy per frame (empty if energies were not recorded).
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

    /// Makes room for `frames` more frames, so recording does not repeatedly reallocate and
    /// copy. Capped at 1 GiB of positions and velocities: a long run that stops early on a
    /// terminal event should not claim memory it never uses.
    pub fn reserve(&mut self, frames: usize) {
        const CAP_BYTES: usize = 1 << 30;
        let per_frame = 2 * self.n_particles.max(1) * std::mem::size_of::<Vec3>();
        let frames = frames.min(CAP_BYTES / per_frame);
        self.t.reserve(frames);
        self.pos.reserve(frames * self.n_particles);
        self.vel.reserve(frames * self.n_particles);
        if self.energies {
            self.kinetic.reserve(frames);
            self.potential.reserve(frames);
        }
    }

    /// Removes all frames and events, keeping `n_particles` and `energies`.
    pub fn clear(&mut self) {
        self.t.clear();
        self.pos.clear();
        self.vel.clear();
        self.kinetic.clear();
        self.potential.clear();
        self.tension.clear();
        self.events.clear();
        self.terminated_by = None;
    }
}

impl Recorder for Trajectory {
    fn frame(&mut self, f: &Frame<'_>) -> Result<()> {
        if f.pos.len() != self.n_particles || f.vel.len() != self.n_particles {
            return invalid(format!(
                "frame has {} particles, trajectory expects {}",
                f.pos.len(),
                self.n_particles
            ));
        }
        if self.t.is_empty() {
            self.n_constraints = f.tension.len();
        } else if f.tension.len() != self.n_constraints {
            return invalid(format!(
                "frame has {} constraints, trajectory expects {}",
                f.tension.len(),
                self.n_constraints
            ));
        }
        if self.energies {
            match (f.kinetic, f.potential) {
                (Some(k), Some(u)) => {
                    self.kinetic.push(k);
                    self.potential.push(u);
                }
                _ => return invalid("trajectory records energies but the frame has none"),
            }
        }
        self.tension.extend_from_slice(f.tension);
        self.t.push(f.t);
        self.pos.extend_from_slice(f.pos);
        self.vel.extend_from_slice(f.vel);
        Ok(())
    }

    fn event(&mut self, hit: EventHit) -> Result<()> {
        self.events.push(hit);
        Ok(())
    }
}
