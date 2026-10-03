//! Adaptive time stepping: Dormand-Prince 5(4) with error control and dense output.
//!
//! [`World::run_adaptive`] integrates to a final time, choosing each step so that the
//! estimated local error stays below `atol + rtol * |y|` in every position and velocity
//! component. Output comes either at every accepted step or at requested times (from the
//! method's fourth-order continuous extension, so no extra force evaluations), and events
//! are located on the same continuous extension.
//!
//! Adaptive Runge-Kutta is not symplectic: on a conservative problem the energy error
//! grows roughly linearly with time instead of staying bounded. It wins when the required
//! step varies a lot (eccentric orbits, close encounters, transients); for very long
//! integrations of smooth conservative motion a fixed-step symplectic scheme is usually better.

use crate::error::{invalid, Result};
use crate::events::{Event, EventHit};
use crate::integrators::Dopri5;
use crate::vec3::Vec3;
use crate::world::{Recorder, RunFailure, Trajectory, World};

/// When [`World::run_adaptive`] records a frame.
#[derive(Debug, Clone, Copy)]
pub enum Output<'a> {
    /// The initial state and every accepted step (the step sizes are then visible in `t`).
    Steps,
    /// Exactly these times, which must be ordered in the direction of integration and lie
    /// between the current time and the final time (inclusive). Frames between steps come
    /// from the continuous extension.
    Times(&'a [f64]),
}

/// Settings for [`World::run_adaptive`].
#[derive(Clone, Copy)]
pub struct AdaptiveOptions<'a> {
    /// Relative tolerance per component of position and velocity.
    pub rtol: f64,
    /// Absolute tolerance per component (same units for positions and velocities, so scale
    /// the problem or rely on `rtol` when the two differ by orders of magnitude).
    pub atol: f64,
    /// First step to try; chosen automatically if `None`.
    pub first_step: Option<f64>,
    /// Largest step allowed (positive).
    pub max_step: f64,
    /// Fail after this many accepted steps.
    pub max_steps: usize,
    pub output: Output<'a>,
    pub events: &'a [Event],
    /// Record kinetic and potential energy with each frame.
    pub energies: bool,
}

impl Default for AdaptiveOptions<'_> {
    fn default() -> Self {
        Self {
            rtol: 1e-9,
            atol: 1e-12,
            first_step: None,
            max_step: f64::INFINITY,
            max_steps: 10_000_000,
            output: Output::Steps,
            events: &[],
            energies: true,
        }
    }
}

/// Work done by an adaptive run.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AdaptiveStats {
    pub accepted: usize,
    pub rejected: usize,
    /// Force evaluations, including the one used to choose the first step.
    pub evaluations: u64,
}

/// Result of [`World::run_adaptive_into`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AdaptiveOutcome {
    /// Index of the terminal event that stopped the run, if any.
    pub terminated_by: Option<usize>,
    pub stats: AdaptiveStats,
}

/// Result of [`World::run_adaptive`].
pub struct AdaptiveRun {
    pub trajectory: Trajectory,
    pub stats: AdaptiveStats,
}

// Step-size controller constants (Hairer's DOPRI5 defaults, with PI stabilisation).
const SAFETY: f64 = 0.9;
const BETA: f64 = 0.04;
const EXPONENT: f64 = 0.2 - BETA * 0.75;
const MIN_FACTOR: f64 = 0.2;
const MAX_FACTOR: f64 = 10.0;

impl World {
    /// Integrates from the current time to `t_end` with adaptive Dormand-Prince 5(4) steps,
    /// collecting frames in memory. The world's own integrator is not used. Events are
    /// located on the continuous extension to ~1e-12 of a step; a terminal event stops the
    /// run there, leaves the world at the event state and, with [`Output::Steps`], records
    /// it as the final frame. Integrating backwards (`t_end` before the current time) works.
    ///
    /// Constraints are not supported. On error the world is left at the last accepted step
    /// and the frames recorded so far come back with the error.
    pub fn run_adaptive(
        &mut self,
        t_end: f64,
        options: &AdaptiveOptions<'_>,
    ) -> std::result::Result<AdaptiveRun, RunFailure> {
        let mut traj = Trajectory::new(self.state.len(), options.energies);
        if let Output::Times(times) = options.output {
            traj.reserve(times.len());
        }
        match self.run_adaptive_into(t_end, options, &mut traj) {
            Ok(outcome) => {
                traj.terminated_by = outcome.terminated_by;
                Ok(AdaptiveRun {
                    trajectory: traj,
                    stats: outcome.stats,
                })
            }
            Err(error) => Err(RunFailure {
                error,
                trajectory: Box::new(traj),
            }),
        }
    }

    /// [`World::run_adaptive`], handing frames and events to `recorder` instead of keeping them.
    pub fn run_adaptive_into(
        &mut self,
        t_end: f64,
        options: &AdaptiveOptions<'_>,
        recorder: &mut dyn Recorder,
    ) -> Result<AdaptiveOutcome> {
        let o = *options;
        self.check_adaptive(t_end, &o)?;
        self.tension.clear();
        for (v, &p) in self.state.vel.iter_mut().zip(&self.state.pinned) {
            if p {
                *v = Vec3::ZERO;
            }
        }
        let t0 = self.state.t;
        let direction = if t_end >= t0 { 1.0 } else { -1.0 };
        let times = match o.output {
            Output::Times(times) => times,
            Output::Steps => &[],
        };
        let mut next_time = 0;
        match o.output {
            Output::Steps => self.emit(recorder, o.energies)?,
            Output::Times(_) => {
                while next_time < times.len() && times[next_time] == t0 {
                    self.emit(recorder, o.energies)?;
                    next_time += 1;
                }
            }
        }

        let mut dp = Dopri5::default();
        let mut stats = AdaptiveStats::default();
        if t_end == t0 {
            return Ok(AdaptiveOutcome {
                terminated_by: None,
                stats,
            });
        }
        let mut g = self.event_values(o.events)?;
        dp.begin(&self.state, &self.forces)?;
        let mut h = match o.first_step {
            Some(h) => h,
            None => dp.initial_step(
                &self.forces,
                &self.state.mass,
                &self.state.pinned,
                direction,
                o.rtol,
                o.atol,
            )?,
        };
        let mut previous_error: f64 = 1e-4;
        let mut rejected_last = false;
        let (mut pos, mut vel) = (Vec::new(), Vec::new());

        let result = loop {
            let t = self.state.t;
            if (t_end - t) * direction <= 0.0 {
                break Ok(None);
            }
            if stats.accepted >= o.max_steps {
                break invalid(format!(
                    "reached max_steps = {} at t = {t} before t_end = {t_end}",
                    o.max_steps
                ));
            }
            h = h.min(o.max_step);
            // Stretch the step to land exactly on t_end rather than leave a sliver.
            let last = (t + 1.01 * h * direction - t_end) * direction >= 0.0;
            if last {
                h = (t_end - t).abs();
            }
            if h <= 10.0 * f64::EPSILON * t.abs().max(f64::MIN_POSITIVE) {
                break invalid(format!(
                    "step size underflow at t = {t} (h = {h:e}): the problem is stiff or singular there, or the tolerances are too tight"
                ));
            }
            if let Err(e) = dp.attempt(
                &self.forces,
                &self.state.mass,
                &self.state.pinned,
                h * direction,
            ) {
                break Err(e);
            }
            let err = dp.error_norm(o.rtol, o.atol);
            let grow = err.powf(EXPONENT);
            if err > 1.0 {
                stats.rejected += 1;
                rejected_last = true;
                h /= (grow / SAFETY).min(1.0 / MIN_FACTOR);
                continue;
            }

            // Accepted.
            stats.accepted += 1;
            let factor = (grow / previous_error.powf(BETA) / SAFETY)
                .clamp(1.0 / MAX_FACTOR, 1.0 / MIN_FACTOR);
            let mut h_next = h / factor;
            if rejected_last {
                h_next = h_next.min(h);
            }
            previous_error = err.max(1e-4);
            rejected_last = false;
            let start_t = t;
            let end_t = if last { t_end } else { t + h * direction };
            if !o.events.is_empty() || !times.is_empty() {
                dp.prepare_dense();
            }
            dp.commit(&mut self.state, end_t);

            // Events, located on the continuous extension.
            let mut stop = None;
            if !o.events.is_empty() {
                let g_new = match self.event_values(o.events) {
                    Ok(g) => g,
                    Err(e) => break Err(e),
                };
                match self.dense_crossings(&dp, o.events, start_t, end_t, &g, &g_new) {
                    Ok(hits) => {
                        for (_, e, hit) in hits {
                            let terminal = o.events[e].terminal;
                            let (t_hit, p_hit, v_hit) = (hit.t, hit.pos.clone(), hit.vel.clone());
                            if let Err(err) = recorder.event(hit) {
                                stop = Some(Err(err));
                                break;
                            }
                            if terminal {
                                stop = Some(Ok((e, t_hit, p_hit, v_hit)));
                                break;
                            }
                        }
                    }
                    Err(e) => break Err(e),
                }
                g = g_new;
            }
            let stop = match stop {
                Some(Err(e)) => break Err(e),
                Some(Ok(s)) => Some(s),
                None => None,
            };

            // Requested output times inside this step (up to a terminal event).
            let horizon = stop.as_ref().map_or(end_t, |s| s.1);
            let mut failed = None;
            while next_time < times.len() && (times[next_time] - horizon) * direction <= 0.0 {
                let tk = times[next_time];
                let r = if tk == end_t {
                    self.emit(recorder, o.energies)
                } else {
                    dp.interpolate((tk - start_t) / (end_t - start_t), &mut pos, &mut vel);
                    self.emit_at(recorder, o.energies, tk, &pos, &vel)
                };
                if let Err(e) = r {
                    failed = Some(e);
                    break;
                }
                next_time += 1;
            }
            if let Some(e) = failed {
                break Err(e);
            }

            if let Some((e, t_hit, p_hit, v_hit)) = stop {
                self.state.t = t_hit;
                self.state.pos = p_hit;
                self.state.vel = v_hit;
                if matches!(o.output, Output::Steps) {
                    if let Err(err) = self.emit(recorder, o.energies) {
                        break Err(err);
                    }
                }
                break Ok(Some(e));
            }
            if matches!(o.output, Output::Steps) {
                if let Err(e) = self.emit(recorder, o.energies) {
                    break Err(e);
                }
            }
            if let Err(e) = dp.begin(&self.state, &self.forces) {
                break Err(e);
            }
            h = h_next;
        };
        stats.evaluations = dp.evaluations();
        result.map(|terminated_by| AdaptiveOutcome {
            terminated_by,
            stats,
        })
    }

    fn check_adaptive(&self, t_end: f64, o: &AdaptiveOptions<'_>) -> Result<()> {
        if !self.constraints.is_empty() {
            return invalid(
                "adaptive stepping does not support constraints; use run() with verlet or yoshida4",
            );
        }
        if !t_end.is_finite() {
            return invalid(format!("t_end must be finite, got {t_end}"));
        }
        let tol_ok = |x: f64| x.is_finite() && x >= 0.0;
        if !tol_ok(o.rtol) || !tol_ok(o.atol) || (o.rtol == 0.0 && o.atol == 0.0) {
            return invalid(format!(
                "rtol and atol must be finite, non-negative and not both zero (got {}, {})",
                o.rtol, o.atol
            ));
        }
        if o.rtol > 0.0 && o.rtol < 10.0 * f64::EPSILON {
            return invalid(format!(
                "rtol = {:e} is below what double precision can deliver (use at least 1e-15)",
                o.rtol
            ));
        }
        if o.max_step.is_nan() || o.max_step <= 0.0 {
            return invalid(format!("max_step must be positive, got {}", o.max_step));
        }
        if let Some(h) = o.first_step {
            if !(h > 0.0 && h.is_finite()) {
                return invalid(format!(
                    "first_step must be positive and finite, got {h} (the direction comes from t_end)"
                ));
            }
        }
        if let Output::Times(times) = o.output {
            let t0 = self.state.t;
            let (lo, hi) = (t0.min(t_end), t0.max(t_end));
            if let Some(t) = times.iter().find(|t| !(**t >= lo && **t <= hi)) {
                return invalid(format!(
                    "output time {t} is outside the run [{t0}, {t_end}]"
                ));
            }
            let sign = if t_end >= t0 { 1.0 } else { -1.0 };
            if times.windows(2).any(|w| (w[1] - w[0]) * sign < 0.0) {
                return invalid("output times must be ordered in the direction of integration");
            }
        }
        Ok(())
    }

    /// Events that fired during the accepted step from `t0` to `t1` (values `g0`, `g1`),
    /// located by Illinois root finding on the continuous extension. Sorted by time, each
    /// with its step fraction and the state there.
    fn dense_crossings(
        &self,
        dp: &Dopri5,
        events: &[Event],
        t0: f64,
        t1: f64,
        g0: &[f64],
        g1: &[f64],
    ) -> Result<Vec<(f64, usize, EventHit)>> {
        const TOLERANCE: f64 = 1e-12;
        let mut hits = Vec::new();
        let (mut pos, mut vel) = (Vec::new(), Vec::new());
        let mass = &self.state.mass;
        for (e, event) in events.iter().enumerate() {
            if !event.triggers(g0[e], g1[e]) {
                continue;
            }
            let time = |theta: f64| {
                if theta == 1.0 {
                    t1
                } else {
                    t0 + theta * (t1 - t0)
                }
            };
            let (mut lo, mut f_lo, mut hi, mut f_hi) = (0.0, g0[e], 1.0, g1[e]);
            let mut last_side = 0;
            for _ in 0..100 {
                if hi - lo <= TOLERANCE || f_hi == 0.0 {
                    break;
                }
                let mut theta = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
                if !(theta > lo && theta < hi) {
                    theta = 0.5 * (lo + hi);
                    if !(theta > lo && theta < hi) {
                        break;
                    }
                }
                dp.interpolate(theta, &mut pos, &mut vel);
                let f = event.function.value(time(theta), &pos, &vel, mass)?;
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
            let (p, v) = if hi == 1.0 {
                (self.state.pos.clone(), self.state.vel.clone())
            } else {
                dp.interpolate(hi, &mut pos, &mut vel);
                (pos.clone(), vel.clone())
            };
            hits.push((
                hi,
                e,
                EventHit {
                    event: e,
                    t: time(hi),
                    pos: p,
                    vel: v,
                },
            ));
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        Ok(hits)
    }
}
