use crate::error::{invalid, Result};
use crate::forces::{Force, ForceSet};
use crate::integrators::{self, Integrator};
use crate::state::State;
use crate::vec3::Vec3;

/// A simulation: particle state + forces + time integrator.
pub struct World {
    pub state: State,
    pub forces: ForceSet,
    integrator: Box<dyn Integrator>,
}

impl World {
    pub fn new(integrator: Box<dyn Integrator>) -> Self {
        Self {
            state: State::new(),
            forces: ForceSet::new(),
            integrator,
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

    pub fn add_particle(&mut self, pos: Vec3, vel: Vec3, mass: f64) -> Result<usize> {
        if !(mass.is_finite() && mass > 0.0) {
            return invalid(format!("mass must be positive and finite, got {mass}"));
        }
        Ok(self.state.add_particle(pos, vel, mass))
    }

    pub fn add_force(&mut self, force: impl Force + 'static) {
        self.forces.add(Box::new(force));
    }

    pub fn step(&mut self, dt: f64) -> Result<()> {
        if !(dt.is_finite() && dt != 0.0) {
            return invalid(format!("dt must be finite and non-zero, got {dt}"));
        }
        self.integrator.step(&mut self.state, &self.forces, dt)?;
        let finite = self
            .state
            .pos
            .iter()
            .chain(&self.state.vel)
            .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite());
        if !finite {
            return invalid(format!(
                "state became non-finite at t = {} (time step too large or singular force?)",
                self.state.t
            ));
        }
        Ok(())
    }

    /// Takes `steps` steps of size `dt`, recording the initial state, every
    /// `record_every`-th step, and the final state.
    pub fn run(&mut self, dt: f64, steps: usize, record_every: usize) -> Result<Trajectory> {
        if record_every == 0 {
            return invalid("record_every must be at least 1");
        }
        let mut traj = Trajectory::new(self.state.len());
        traj.record(self)?;
        for k in 1..=steps {
            self.step(dt)?;
            if k % record_every == 0 || k == steps {
                traj.record(self)?;
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

    fn record(&mut self, world: &World) -> Result<()> {
        self.t.push(world.state.t);
        self.pos.extend_from_slice(&world.state.pos);
        self.vel.extend_from_slice(&world.state.vel);
        self.kinetic.push(world.kinetic_energy());
        self.potential.push(world.potential_energy()?);
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
