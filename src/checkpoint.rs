//! Checkpoints: everything needed to continue a run later, bit for bit.
//!
//! A [`Checkpoint`] is plain data (the state, the integrator's name and descriptions of the
//! forces), so it can be written in any format; the Python package saves it as `.npz`.
//! Restoring and continuing gives exactly the same numbers as never having stopped, because
//! the integrators keep no history between steps other than caches of values they would
//! otherwise recompute identically.

use crate::constraints::Constraints;
use crate::error::{invalid, Result};
use crate::forces::{BuiltinForce, Force, ForceId, ForceSet};
use crate::integrators::{self, Scheme};
use crate::state::State;
use crate::world::World;

/// A force as stored in a checkpoint.
#[derive(Debug, Clone)]
pub enum SavedForce {
    /// A built-in force, stored completely.
    Builtin(BuiltinForce),
    /// Any other force (e.g. one defined in Python): only its name is stored, and the
    /// force itself must be supplied again when the checkpoint is restored.
    External { name: String },
}

/// The complete state of a [`World`].
#[derive(Debug, Clone)]
pub struct Checkpoint {
    pub state: State,
    /// Integrator name, as accepted by [`crate::integrators::by_name`].
    pub integrator: String,
    /// Coefficients of a user-defined integrator (then used instead of `integrator`).
    pub integrator_scheme: Option<Scheme>,
    /// Forces in evaluation order, with their ids.
    pub forces: Vec<(ForceId, SavedForce)>,
    /// The id the next added force will get.
    pub next_force_id: ForceId,
    /// Rigid rods, with their ids and the solver settings.
    pub constraints: Constraints,
}

impl World {
    pub fn checkpoint(&self) -> Checkpoint {
        Checkpoint {
            state: self.state.clone(),
            integrator: self.integrator().name().to_string(),
            integrator_scheme: self.integrator().scheme(),
            forces: self
                .forces
                .iter()
                .map(|(id, f)| {
                    let saved = match f.builtin() {
                        Some(b) => SavedForce::Builtin(b),
                        None => SavedForce::External { name: f.name() },
                    };
                    (id, saved)
                })
                .collect(),
            next_force_id: self.forces.next_id(),
            constraints: self.constraints.clone(),
        }
    }

    /// Rebuilds a world from a checkpoint. `supply(id, name)` is called for every
    /// [`SavedForce::External`] force and must return it.
    pub fn from_checkpoint(
        checkpoint: Checkpoint,
        mut supply: impl FnMut(ForceId, &str) -> Result<Box<dyn Force>>,
    ) -> Result<World> {
        let Checkpoint {
            state,
            integrator,
            integrator_scheme,
            forces,
            next_force_id,
            constraints,
        } = checkpoint;
        let n = state.pos.len();
        if [state.vel.len(), state.mass.len(), state.pinned.len()] != [n, n, n] {
            return invalid(format!(
                "inconsistent checkpoint: {n} positions, {} velocities, {} masses, {} pinned flags",
                state.vel.len(),
                state.mass.len(),
                state.pinned.len()
            ));
        }
        let mut world = World::new(match integrator_scheme {
            Some(scheme) => scheme.build()?,
            None => integrators::by_name(&integrator)?,
        });
        let forces = forces
            .into_iter()
            .map(|(id, saved)| {
                let force = match saved {
                    SavedForce::Builtin(b) => b.into_force(),
                    SavedForce::External { name } => supply(id, &name)?,
                };
                Ok((id, force))
            })
            .collect::<Result<Vec<_>>>()?;
        world.forces = ForceSet::from_parts(forces, next_force_id)?;
        // Validate through the public setters so a corrupt file cannot bypass the checks.
        world.state.t = state.t;
        for ((&pos, &vel), &mass) in state.pos.iter().zip(&state.vel).zip(&state.mass) {
            world.add_particle(pos, vel, mass)?;
        }
        for (i, &p) in state.pinned.iter().enumerate() {
            if p {
                world.pin(i, true)?;
            }
        }
        // Pinning zeroes velocities; a valid checkpoint already has them at zero.
        world.set_velocities(state.vel)?;
        constraints.validate(&world.state)?;
        world.tension = vec![0.0; constraints.len()];
        world.constraints = constraints;
        Ok(world)
    }
}
