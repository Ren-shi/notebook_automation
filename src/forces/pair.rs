//! Short-range pair potentials with a cutoff, optionally in a periodic box: the building
//! block of molecular dynamics.

use std::sync::{Arc, Mutex};

use super::{BuiltinForce, Force, Param};
use crate::error::{invalid, Result};
use crate::parallel;
use crate::vec3::Vec3;

/// The interaction law `V(r)` of a [`PairPotential`].
#[derive(Debug, Clone, PartialEq)]
pub enum PairKind {
    /// `4ε [(σ/r)¹² − (σ/r)⁶]`.
    LennardJones { epsilon: f64, sigma: f64 },
    /// `D [(1 − e^{−a (r − r0)})² − 1]`.
    Morse { depth: f64, a: f64, r0: f64 },
    /// `V` tabulated at equally spaced `r = r_min + k Δr`, with the derivative `dV/dr` at the
    /// same points; interpolated by cubic Hermite polynomials (so the force is continuous),
    /// and continued linearly below `r_min`.
    Table {
        r_min: f64,
        dr: f64,
        v: Vec<f64>,
        dv: Vec<f64>,
    },
}

impl PairKind {
    /// `(V(r), −V'(r) / r)`: the potential and the force factor (force on j is `f · d`
    /// with `d = x_j − x_i`).
    fn eval(&self, r2: f64) -> (f64, f64) {
        match *self {
            PairKind::LennardJones { epsilon, sigma } => {
                let s2 = sigma * sigma / r2;
                let s6 = s2 * s2 * s2;
                let s12 = s6 * s6;
                (
                    4.0 * epsilon * (s12 - s6),
                    24.0 * epsilon * (2.0 * s12 - s6) / r2,
                )
            }
            PairKind::Morse { depth, a, r0 } => {
                let r = r2.sqrt();
                let e = (-a * (r - r0)).exp();
                let v = depth * ((1.0 - e) * (1.0 - e) - 1.0);
                let dv = 2.0 * depth * a * e * (1.0 - e);
                (v, -dv / r)
            }
            PairKind::Table {
                r_min,
                dr,
                ref v,
                ref dv,
            } => {
                let r = r2.sqrt();
                if r < r_min {
                    // Below the table: continue linearly from its first point.
                    return (v[0] + dv[0] * (r - r_min), -dv[0] / r);
                }
                let last = v.len() - 1;
                let x = ((r - r_min) / dr).min(last as f64);
                let k = (x.floor() as usize).min(last - 1);
                let t = x - k as f64;
                // Cubic Hermite on [k, k+1].
                let (p0, p1, m0, m1) = (v[k], v[k + 1], dv[k] * dr, dv[k + 1] * dr);
                let (t2, t3) = (t * t, t * t * t);
                let val = (2.0 * t3 - 3.0 * t2 + 1.0) * p0
                    + (t3 - 2.0 * t2 + t) * m0
                    + (-2.0 * t3 + 3.0 * t2) * p1
                    + (t3 - t2) * m1;
                let deriv = ((6.0 * t2 - 6.0 * t) * p0
                    + (3.0 * t2 - 4.0 * t + 1.0) * m0
                    + (-6.0 * t2 + 6.0 * t) * p1
                    + (3.0 * t2 - 2.0 * t) * m1)
                    / dr;
                (val, -deriv / r)
            }
        }
    }

    fn validate(&self) -> Result<()> {
        let ok = match self {
            PairKind::LennardJones { epsilon, sigma } => {
                epsilon.is_finite() && *sigma > 0.0 && sigma.is_finite()
            }
            PairKind::Morse { depth, a, r0 } => {
                depth.is_finite() && *a > 0.0 && a.is_finite() && r0.is_finite()
            }
            PairKind::Table { r_min, dr, v, dv } => {
                *r_min > 0.0
                    && *dr > 0.0
                    && v.len() >= 2
                    && v.len() == dv.len()
                    && v.iter().chain(dv).all(|x| x.is_finite())
            }
        };
        if ok {
            Ok(())
        } else {
            invalid(format!("invalid pair potential parameters: {self:?}"))
        }
    }
}

/// A pair potential acting between every pair of particles closer than `cutoff`, optionally
/// shifted so `V(cutoff) = 0` (the energy is then continuous; the force still jumps), and
/// optionally periodic in the box `[0, L)³` (minimum-image convention; coordinates need not
/// be wrapped). Neighbours come from cell lists, so the cost is O(N) at fixed density; each
/// particle gathers its own force, so the evaluation is parallel and deterministic.
#[derive(Debug)]
pub struct PairPotential {
    pub kind: PairKind,
    pub cutoff: f64,
    pub shift: bool,
    /// Periodic box lengths (`None`: open boundaries).
    pub period: Option<Vec3>,
    /// Verlet neighbour list (a cache: results never depend on it).
    list: Mutex<Option<Arc<NeighbourList>>>,
}

impl Clone for PairPotential {
    fn clone(&self) -> Self {
        Self::new(self.kind.clone(), self.cutoff, self.shift, self.period)
    }
}

/// Below this many particles the force loop runs serially over pairs (half the work).
const SERIAL_BELOW: usize = 8192;

/// Verlet list: every pair closer than `cutoff + skin` when it was built, each particle's
/// neighbours sorted by index. Valid while no particle has moved more than `skin / 2`.
#[derive(Debug)]
struct NeighbourList {
    built_at: Vec<Vec3>,
    skin: f64,
    cutoff: f64,
    period: Option<Vec3>,
    neighbours: Vec<Vec<u32>>,
    /// The same lists restricted to `j > i` (each pair once).
    half: Vec<Vec<u32>>,
}

/// Neighbour search: a cell list over the (wrapped) positions.
struct Cells {
    dims: [usize; 3],
    size: Vec3,
    /// Particle indices sorted by cell; `starts[c]..starts[c + 1]`.
    order: Vec<usize>,
    starts: Vec<usize>,
    wrapped: Vec<Vec3>,
    periodic: bool,
}

impl Cells {
    fn build(pos: &[Vec3], cutoff: f64, period: Option<Vec3>) -> Option<Cells> {
        let (lo, extent, periodic) = match period {
            Some(l) => (Vec3::ZERO, l, true),
            None => {
                let (mut lo, mut hi) = (pos[0], pos[0]);
                for p in pos {
                    lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
                    hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
                }
                (lo, hi - lo + Vec3::new(1e-9, 1e-9, 1e-9), false)
            }
        };
        let mut dims = [0usize; 3];
        for a in 0..3 {
            dims[a] = ((extent[a] / cutoff).floor() as usize).clamp(1, 1 << 20);
            // Periodic cell lists need at least 3 cells per side to avoid double counting.
            if periodic && dims[a] < 3 {
                return None;
            }
        }
        let n_cells = dims[0].checked_mul(dims[1])?.checked_mul(dims[2])?;
        if n_cells > 64 * pos.len() + 64 {
            return None; // sparse system: all pairs is cheaper
        }
        let size = Vec3::new(
            extent.x / dims[0] as f64,
            extent.y / dims[1] as f64,
            extent.z / dims[2] as f64,
        );
        let wrapped: Vec<Vec3> = pos
            .iter()
            .map(|p| {
                if periodic {
                    Vec3::new(
                        p.x.rem_euclid(extent.x),
                        p.y.rem_euclid(extent.y),
                        p.z.rem_euclid(extent.z),
                    )
                } else {
                    *p - lo
                }
            })
            .collect();
        let cell_of = |p: Vec3| {
            let c = |a: usize| ((p[a] / size[a]) as usize).min(dims[a] - 1);
            (c(2) * dims[1] + c(1)) * dims[0] + c(0)
        };
        let mut counts = vec![0usize; n_cells + 1];
        for p in &wrapped {
            counts[cell_of(*p) + 1] += 1;
        }
        for c in 0..n_cells {
            counts[c + 1] += counts[c];
        }
        let starts = counts.clone();
        let mut fill = counts;
        let mut order = vec![0; pos.len()];
        for (i, p) in wrapped.iter().enumerate() {
            let c = cell_of(*p);
            order[fill[c]] = i;
            fill[c] += 1;
        }
        Some(Cells {
            dims,
            size,
            order,
            starts,
            wrapped,
            periodic,
        })
    }

    /// Calls `f(j)` for every particle in the cells around particle `i` (including `i`).
    fn neighbours(&self, i: usize, mut f: impl FnMut(usize)) {
        let p = self.wrapped[i];
        let c = |a: usize| ((p[a] / self.size[a]) as isize).min(self.dims[a] as isize - 1);
        let (cx, cy, cz) = (c(0), c(1), c(2));
        let d = self.dims.map(|x| x as isize);
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let (mut x, mut y, mut z) = (cx + dx, cy + dy, cz + dz);
                    if self.periodic {
                        x = x.rem_euclid(d[0]);
                        y = y.rem_euclid(d[1]);
                        z = z.rem_euclid(d[2]);
                    } else if x < 0 || y < 0 || z < 0 || x >= d[0] || y >= d[1] || z >= d[2] {
                        continue;
                    }
                    let cell = ((z * d[1] + y) * d[0] + x) as usize;
                    for &j in &self.order[self.starts[cell]..self.starts[cell + 1]] {
                        f(j);
                    }
                }
            }
        }
    }
}

impl PairPotential {
    pub fn new(kind: PairKind, cutoff: f64, shift: bool, period: Option<Vec3>) -> Self {
        Self {
            kind,
            cutoff,
            shift,
            period,
            list: Mutex::new(None),
        }
    }

    fn check(&self) -> Result<()> {
        self.kind.validate()?;
        if !(self.cutoff > 0.0 && self.cutoff.is_finite()) {
            return invalid(format!("cutoff must be positive, got {}", self.cutoff));
        }
        if let Some(l) = self.period {
            if !(l.x > 0.0
                && l.y > 0.0
                && l.z > 0.0
                && l.x.is_finite()
                && l.y.is_finite()
                && l.z.is_finite())
            {
                return invalid(format!("periodic box lengths must be positive, got {l:?}"));
            }
            if 2.0 * self.cutoff > l.x.min(l.y).min(l.z) {
                return invalid(format!(
                    "cutoff {} exceeds half the smallest box length {:?} (minimum image)",
                    self.cutoff, l
                ));
            }
        }
        if let PairKind::Table { r_min, dr, v, .. } = &self.kind {
            if r_min + dr * (v.len() - 1) as f64 + 1e-12 < self.cutoff {
                return invalid("the table must extend to the cutoff");
            }
        }
        Ok(())
    }

    /// Separation vector from `i` to `j` under the minimum-image convention.
    fn separation(&self, a: Vec3, b: Vec3) -> Vec3 {
        let d = b - a;
        match self.period {
            Some(l) => {
                let wrap = |x: f64, l: f64| {
                    let half = 0.5 * l;
                    if x > half || x < -half {
                        x - l * (x / l).round()
                    } else {
                        x
                    }
                };
                Vec3::new(wrap(d.x, l.x), wrap(d.y, l.y), wrap(d.z, l.z))
            }
            None => d,
        }
    }

    fn energy_shift(&self) -> f64 {
        if self.shift {
            self.kind.eval(self.cutoff * self.cutoff).0
        } else {
            0.0
        }
    }

    /// Builds the Verlet list for these positions (cell lists, or all pairs if the box is
    /// too small for them).
    fn build_list(&self, pos: &[Vec3]) -> NeighbourList {
        let skin = 0.12 * self.cutoff;
        let reach = self.cutoff + skin;
        let reach2 = reach * reach;
        let n = pos.len();
        let cells = if n > 1 {
            Cells::build(pos, reach, self.period)
        } else {
            None
        };
        let mut neighbours: Vec<Vec<u32>> = vec![Vec::new(); n];
        parallel::for_each_indexed(&mut neighbours, |i, list| {
            let mut consider = |j: usize| {
                if j != i && self.separation(pos[i], pos[j]).norm_squared() < reach2 {
                    list.push(j as u32);
                }
            };
            match &cells {
                Some(c) => c.neighbours(i, consider),
                None => (0..n).for_each(&mut consider),
            }
            list.sort_unstable();
            list.dedup();
        });
        let half = neighbours
            .iter()
            .enumerate()
            .map(|(i, l)| l.iter().copied().filter(|&j| j as usize > i).collect())
            .collect();
        NeighbourList {
            built_at: pos.to_vec(),
            skin,
            cutoff: self.cutoff,
            period: self.period,
            neighbours,
            half,
        }
    }

    /// The current Verlet list, rebuilt if any particle moved more than half the skin.
    fn list(&self, pos: &[Vec3]) -> Arc<NeighbourList> {
        let mut guard = self.list.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(list) = guard.as_ref() {
            let fresh = list.built_at.len() == pos.len()
                && list.cutoff == self.cutoff
                && list.period == self.period
                && {
                    let limit = 0.25 * list.skin * list.skin;
                    list.built_at
                        .iter()
                        .zip(pos)
                        .all(|(a, b)| self.separation(*a, *b).norm_squared() < limit)
                };
            if fresh {
                return list.clone();
            }
        }
        let list = Arc::new(self.build_list(pos));
        *guard = Some(list.clone());
        list
    }

    /// Calls `visit(j, d, r²)` for every neighbour `j` of `i` within the cutoff, in
    /// ascending `j` (so sums do not depend on when the list was built).
    fn for_each_neighbour(
        &self,
        pos: &[Vec3],
        list: &NeighbourList,
        i: usize,
        mut visit: impl FnMut(usize, Vec3, f64),
    ) {
        let rc2 = self.cutoff * self.cutoff;
        for &j in &list.neighbours[i] {
            let j = j as usize;
            let d = self.separation(pos[i], pos[j]);
            let r2 = d.norm_squared();
            if r2 < rc2 && r2 > 0.0 {
                visit(j, d, r2);
            }
        }
    }

    /// Force on every particle (mass times acceleration). Small systems use each pair once
    /// (Newton's third law, serial); large ones gather per particle in parallel. Both orders
    /// depend only on the positions, never on when the list was built.
    fn forces(&self, pos: &[Vec3]) -> Vec<Vec3> {
        let list = self.list(pos);
        let mut force = vec![Vec3::ZERO; pos.len()];
        if pos.len() < SERIAL_BELOW {
            let rc2 = self.cutoff * self.cutoff;
            for i in 0..pos.len() {
                let mut fi = Vec3::ZERO;
                for &j in &list.half[i] {
                    let j = j as usize;
                    let d = self.separation(pos[i], pos[j]);
                    let r2 = d.norm_squared();
                    if r2 < rc2 && r2 > 0.0 {
                        let f = d * self.kind.eval(r2).1;
                        fi -= f;
                        force[j] += f;
                    }
                }
                force[i] += fi;
            }
        } else {
            parallel::for_each_indexed(&mut force, |i, f| {
                self.for_each_neighbour(pos, &list, i, |_, d, r2| {
                    *f -= d * self.kind.eval(r2).1;
                });
            });
        }
        force
    }

    /// Total potential energy and virial.
    fn energy_and_virial(&self, pos: &[Vec3]) -> (f64, f64) {
        let list = self.list(pos);
        let shift = self.energy_shift();
        let mut per: Vec<(f64, f64)> = vec![(0.0, 0.0); pos.len()];
        parallel::for_each_indexed(&mut per, |i, o| {
            self.for_each_neighbour(pos, &list, i, |_, _, r2| {
                let (v, fac) = self.kind.eval(r2);
                o.0 += 0.5 * (v - shift);
                o.1 += 0.5 * fac * r2;
            });
        });
        per.iter().fold((0.0, 0.0), |(u, w), (a, b)| (u + a, w + b))
    }
}

impl Force for PairPotential {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        self.check()?;
        let force = self.forces(pos);
        for (i, (a, f)) in acc.iter_mut().zip(force).enumerate() {
            if f != Vec3::ZERO {
                if mass[i] == 0.0 {
                    return invalid(format!(
                        "PairPotential: particle {i} has neighbours but no mass"
                    ));
                }
                *a += f / mass[i];
            }
        }
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], _mass: &[f64]) -> Result<Option<f64>> {
        self.check()?;
        Ok(Some(self.energy_and_virial(pos).0))
    }

    fn virial(&self, pos: &[Vec3]) -> Result<Option<f64>> {
        self.check()?;
        Ok(Some(self.energy_and_virial(pos).1))
    }

    fn name(&self) -> String {
        match self.kind {
            PairKind::LennardJones { .. } => "LennardJones".into(),
            PairKind::Morse { .. } => "Morse".into(),
            PairKind::Table { .. } => "TabulatedPair".into(),
        }
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::PairPotential(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        let mut p = vec![("cutoff", Param::Scalar(self.cutoff))];
        match self.kind {
            PairKind::LennardJones { epsilon, sigma } => {
                p.push(("epsilon", Param::Scalar(epsilon)));
                p.push(("sigma", Param::Scalar(sigma)));
            }
            PairKind::Morse { depth, a, r0 } => {
                p.push(("depth", Param::Scalar(depth)));
                p.push(("a", Param::Scalar(a)));
                p.push(("r0", Param::Scalar(r0)));
            }
            PairKind::Table { .. } => {}
        }
        if let Some(l) = self.period {
            p.push(("box", Param::Vector(l)));
        }
        p
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        let mut next = self.clone();
        let scalar = |v: Param| match v {
            Param::Scalar(x) if x.is_finite() => Ok(x),
            _ => invalid(format!("{name} must be a finite scalar")),
        };
        match (name, &mut next.kind) {
            ("cutoff", _) => next.cutoff = scalar(value)?,
            ("epsilon", PairKind::LennardJones { epsilon, .. }) => *epsilon = scalar(value)?,
            ("sigma", PairKind::LennardJones { sigma, .. }) => *sigma = scalar(value)?,
            ("depth", PairKind::Morse { depth, .. }) => *depth = scalar(value)?,
            ("a", PairKind::Morse { a, .. }) => *a = scalar(value)?,
            ("r0", PairKind::Morse { r0, .. }) => *r0 = scalar(value)?,
            ("box", _) => match value {
                Param::Vector(l) => next.period = Some(l),
                Param::Scalar(_) => return invalid("box is a 3-vector"),
            },
            _ => return invalid(format!("{} has no parameter {name:?}", self.name())),
        }
        next.check()?;
        *self = next;
        Ok(())
    }
}
