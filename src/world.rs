use std::fmt;

use crate::error::{invalid, Result, SimError};
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
        for k in 1..=steps {
            if let Err(e) = self.step(dt) {
                if traj.t.last() != Some(&self.state.t) {
                    let _ = traj.record(self);
                }
                return Err(fail(e, traj));
            }
            if k % record_every == 0 || k == steps {
                if let Err(e) = traj.record(self) {
                    return Err(fail(e, traj));
                }
            }
        }
        Ok(traj)
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
