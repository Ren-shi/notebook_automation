//! Symplectic schemes built from the exact flows of the kinetic and potential parts:
//! compositions of velocity Verlet substeps, and general kick/drift splittings.

use super::{kdk, rattle, AccelCache, Integrator};
use crate::constraints::Constraints;
use crate::error::{invalid, Result};
use crate::forces::ForceSet;
use crate::state::State;
use crate::vec3::Vec3;

/// One flow of a splitting scheme, for a fraction of the step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Op {
    /// `x += c h v` (time advances by `c h`).
    Drift(f64),
    /// `v += c h a(t, x)`.
    Kick(f64),
}

/// The coefficients of a user-defined scheme: what checkpoints store to rebuild it.
#[derive(Debug, Clone, PartialEq)]
pub enum Scheme {
    /// Velocity Verlet substeps of `weights[k] * dt` in turn.
    Composition {
        name: String,
        order: u32,
        weights: Vec<f64>,
    },
    /// Kicks and drifts in the given order.
    Splitting {
        name: String,
        order: u32,
        ops: Vec<Op>,
    },
}

impl Scheme {
    /// Builds the integrator. Names of built-in integrators are reserved.
    pub fn build(&self) -> Result<Box<dyn Integrator>> {
        let name = match self {
            Scheme::Composition { name, .. } | Scheme::Splitting { name, .. } => name,
        };
        if super::by_name(name).is_ok() {
            return invalid(format!(
                "{name:?} is the name of a built-in integrator; choose another"
            ));
        }
        Ok(match self {
            Scheme::Composition {
                name,
                order,
                weights,
            } => Box::new(Composition::new(name, *order, weights.clone())?),
            Scheme::Splitting { name, order, ops } => {
                Box::new(Splitting::new(name, *order, ops.clone())?)
            }
        })
    }
}

fn check_sum(what: &str, values: impl Iterator<Item = f64>) -> Result<()> {
    let mut sum = 0.0;
    for v in values {
        if !v.is_finite() {
            return invalid(format!("{what} must be finite, got {v}"));
        }
        sum += v;
    }
    if (sum - 1.0).abs() > 1e-12 {
        return invalid(format!("{what} must sum to 1 (consistency), got {sum}"));
    }
    Ok(())
}

/// A composition of velocity Verlet substeps with the given weights (which must sum to 1).
/// Symmetric weights give a time-reversible symplectic scheme; Yoshida's construction
/// raises the order by two with each level. Supports constraints (RATTLE substeps).
pub struct Composition {
    name: String,
    order: u32,
    weights: Vec<f64>,
    cache: AccelCache,
    old: Vec<Vec3>,
}

impl Composition {
    pub fn new(name: impl Into<String>, order: u32, weights: Vec<f64>) -> Result<Self> {
        if weights.is_empty() {
            return invalid("a composition needs at least one weight");
        }
        check_sum("composition weights", weights.iter().copied())?;
        Ok(Self {
            name: name.into(),
            order,
            weights,
            cache: AccelCache::default(),
            old: Vec::new(),
        })
    }

    /// Yoshida (1990) sixth order, solution A: 7 substeps.
    pub fn yoshida6() -> Self {
        let w = [-1.17767998417887, 0.235573213359357, 0.784513610477560];
        Self::symmetric("yoshida6", 6, &w)
    }

    /// Yoshida (1990) eighth order, solution D: 15 substeps.
    pub fn yoshida8() -> Self {
        let w = [
            0.102799849391985,
            -1.96061023297549,
            1.93813913762276,
            -0.158240635368243,
            -1.44485223686048,
            0.253693336566229,
            0.914844246229740,
        ];
        Self::symmetric("yoshida8", 8, &w)
    }

    /// `[w_m .. w_1, w_0, w_1 .. w_m]` with `w_0 = 1 - 2 Σ w_k`, from `[w_1 .. w_m]`.
    fn symmetric(name: &str, order: u32, w: &[f64]) -> Self {
        let w0 = 1.0 - 2.0 * w.iter().sum::<f64>();
        let mut weights: Vec<f64> = w.iter().rev().copied().collect();
        weights.push(w0);
        weights.extend(w.iter().copied());
        Self::new(name, order, weights).expect("built-in weights are consistent")
    }
}

impl Integrator for Composition {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        for &w in &self.weights {
            kdk(s, forces, &mut self.cache, w * dt)?;
        }
        Ok(())
    }
    fn step_constrained(
        &mut self,
        s: &mut State,
        forces: &ForceSet,
        constraints: &Constraints,
        tension: &mut Vec<f64>,
        dt: f64,
    ) -> Result<()> {
        for &w in &self.weights {
            rattle(
                s,
                forces,
                constraints,
                &mut self.cache,
                &mut self.old,
                tension,
                w * dt,
            )?;
        }
        Ok(())
    }
    fn name(&self) -> &str {
        &self.name
    }
    fn order(&self) -> u32 {
        self.order
    }
    fn symplectic(&self) -> bool {
        true
    }
    fn scheme(&self) -> Option<Scheme> {
        match self.name.as_str() {
            "yoshida6" | "yoshida8" => None,
            _ => Some(Scheme::Composition {
                name: self.name.clone(),
                order: self.order,
                weights: self.weights.clone(),
            }),
        }
    }
}

/// A general splitting: kicks and drifts applied in order, each for a fraction of the step
/// (drift fractions and kick fractions must each sum to 1). Symplectic for
/// velocity-independent forces. A kick at the same positions as the previous one (e.g. the
/// last kick of one step and the first of the next) reuses its force evaluation.
pub struct Splitting {
    name: String,
    order: u32,
    ops: Vec<Op>,
    cache: AccelCache,
}

impl Splitting {
    pub fn new(name: impl Into<String>, order: u32, ops: Vec<Op>) -> Result<Self> {
        let drifts = ops.iter().filter_map(|o| match o {
            Op::Drift(c) => Some(*c),
            Op::Kick(_) => None,
        });
        check_sum("drift coefficients", drifts)?;
        let kicks = ops.iter().filter_map(|o| match o {
            Op::Kick(c) => Some(*c),
            Op::Drift(_) => None,
        });
        check_sum("kick coefficients", kicks)?;
        Ok(Self {
            name: name.into(),
            order,
            ops,
            cache: AccelCache::default(),
        })
    }

    /// Omelyan, Mryglod & Folk (2002) position-extended Forest-Ruth-like scheme (PEFRL):
    /// fourth order, 4 force evaluations, error constant ~100 times smaller than Yoshida's.
    pub fn pefrl() -> Self {
        let xi = 0.1786178958448091;
        let lambda = -0.2123418310626054;
        let chi = -0.066_264_582_669_818_5;
        let ops = vec![
            Op::Drift(xi),
            Op::Kick(0.5 * (1.0 - 2.0 * lambda)),
            Op::Drift(chi),
            Op::Kick(lambda),
            Op::Drift(1.0 - 2.0 * (chi + xi)),
            Op::Kick(lambda),
            Op::Drift(chi),
            Op::Kick(0.5 * (1.0 - 2.0 * lambda)),
            Op::Drift(xi),
        ];
        Self::new("pefrl", 4, ops).expect("built-in coefficients are consistent")
    }

    /// Blanes & Moan (2002) optimised fourth-order Runge-Kutta-Nyström splitting
    /// (SRKN₆ᵇ): 6 force evaluations, error constant far below Yoshida's.
    pub fn blanes_moan4() -> Self {
        let b1 = 0.0792036964311957;
        let b2 = 0.353172906049774;
        let b3 = -0.0420650803577195;
        let b4 = 1.0 - 2.0 * (b1 + b2 + b3);
        let a1 = 0.209515106613362;
        let a2 = -0.143851773179818;
        let a3 = 0.5 - (a1 + a2);
        let ops = vec![
            Op::Kick(b1),
            Op::Drift(a1),
            Op::Kick(b2),
            Op::Drift(a2),
            Op::Kick(b3),
            Op::Drift(a3),
            Op::Kick(b4),
            Op::Drift(a3),
            Op::Kick(b3),
            Op::Drift(a2),
            Op::Kick(b2),
            Op::Drift(a1),
            Op::Kick(b1),
        ];
        Self::new("blanes_moan4", 4, ops).expect("built-in coefficients are consistent")
    }
}

impl Integrator for Splitting {
    fn step(&mut self, s: &mut State, forces: &ForceSet, dt: f64) -> Result<()> {
        let t0 = s.t;
        let mut elapsed = 0.0;
        for op in &self.ops {
            match *op {
                Op::Drift(c) => {
                    for (x, v) in s.pos.iter_mut().zip(&s.vel) {
                        *x += *v * (c * dt);
                    }
                    elapsed += c;
                    s.t = t0 + elapsed * dt;
                }
                Op::Kick(c) => {
                    let acc = self.cache.get(s, forces)?;
                    for (v, a) in s.vel.iter_mut().zip(acc) {
                        *v += *a * (c * dt);
                    }
                }
            }
        }
        s.t = t0 + dt;
        Ok(())
    }
    fn name(&self) -> &str {
        &self.name
    }
    fn order(&self) -> u32 {
        self.order
    }
    fn symplectic(&self) -> bool {
        true
    }
    fn scheme(&self) -> Option<Scheme> {
        match self.name.as_str() {
            "pefrl" | "blanes_moan4" => None,
            _ => Some(Scheme::Splitting {
                name: self.name.clone(),
                order: self.order,
                ops: self.ops.clone(),
            }),
        }
    }
}
