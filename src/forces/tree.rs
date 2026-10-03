//! Barnes-Hut tree gravity: O(N log N) instead of O(N²).
//!
//! Every evaluation builds an octree over the particles. Each node stores its total mass,
//! centre of mass and (optionally) traceless quadrupole moment. For each particle the tree
//! is walked from the root: a node far enough away is replaced by its multipole expansion,
//! otherwise it is opened; leaves are summed directly. A node of size `s` whose centre of
//! mass is at distance `d` is accepted when `d > s/θ + δ`, with `δ` the offset between the
//! centre of mass and the cell centre (which guards against particles inside lopsided cells).
//!
//! The force error grows with the opening angle θ (roughly as θ² with monopoles and faster
//! with quadrupoles). Forces are no longer exactly pairwise antisymmetric, so total momentum is
//! conserved only to that accuracy.

use super::{non_negative, scalar, BuiltinForce, Force, Param};
use crate::error::{invalid, Result};
use crate::parallel;
use crate::vec3::Vec3;

/// Particles per leaf.
const LEAF_SIZE: usize = 8;
/// Particles that share one interaction list (the group walk).
const GROUP_SIZE: usize = 32;
/// Deeper than this, coincident particles stay together in one leaf.
const MAX_DEPTH: usize = 48;

#[derive(Clone, Copy)]
struct Node {
    /// Geometric centre and half-width of the cell.
    center: Vec3,
    half: f64,
    mass: f64,
    com: Vec3,
    /// Traceless quadrupole about the centre of mass: xx, yy, zz, xy, xz, yz.
    quad: [f64; 6],
    /// Acceptance distance `s/θ + δ`, squared.
    open2: f64,
    /// Particle range `order[start..end]`.
    start: usize,
    end: usize,
    /// Index of the first child; children are contiguous, `n_children` of them. 0 for a leaf.
    first_child: usize,
    n_children: usize,
}

struct Tree {
    nodes: Vec<Node>,
    /// Particle indices, grouped so every node covers a contiguous range.
    order: Vec<usize>,
}

impl Tree {
    fn build(pos: &[Vec3], mass: &[f64], theta: f64) -> Tree {
        let n = pos.len();
        let (mut lo, mut hi) = (pos[0], pos[0]);
        for p in pos {
            lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let center = (lo + hi) * 0.5;
        let half = 0.5 * (hi.x - lo.x).max(hi.y - lo.y).max(hi.z - lo.z) * (1.0 + 1e-12) + 1e-300;
        let mut tree = Tree {
            nodes: Vec::with_capacity(2 * n / LEAF_SIZE + 16),
            order: (0..n).collect(),
        };
        tree.nodes.push(Node {
            center,
            half,
            mass: 0.0,
            com: Vec3::ZERO,
            quad: [0.0; 6],
            open2: 0.0,
            start: 0,
            end: n,
            first_child: 0,
            n_children: 0,
        });
        tree.split(0, pos, 0);
        tree.moments(0, pos, mass, theta);
        tree
    }

    /// Recursively splits node `k` into octants.
    fn split(&mut self, k: usize, pos: &[Vec3], depth: usize) {
        let Node {
            center,
            half,
            start,
            end,
            ..
        } = self.nodes[k];
        if end - start <= LEAF_SIZE || depth >= MAX_DEPTH {
            return;
        }
        let octant = |p: Vec3| {
            (usize::from(p.x >= center.x))
                | (usize::from(p.y >= center.y) << 1)
                | (usize::from(p.z >= center.z) << 2)
        };
        // Counting sort of this range by octant.
        let mut counts = [0usize; 8];
        for &i in &self.order[start..end] {
            counts[octant(pos[i])] += 1;
        }
        let mut offsets = [0usize; 8];
        let mut acc = start;
        for o in 0..8 {
            offsets[o] = acc;
            acc += counts[o];
        }
        let mut sorted = vec![0usize; end - start];
        let mut fill = offsets;
        for &i in &self.order[start..end] {
            let o = octant(pos[i]);
            sorted[fill[o] - start] = i;
            fill[o] += 1;
        }
        self.order[start..end].copy_from_slice(&sorted);

        let first = self.nodes.len();
        let mut children = 0;
        for o in 0..8 {
            if counts[o] == 0 {
                continue;
            }
            let q = half * 0.5;
            let c = center
                + Vec3::new(
                    if o & 1 != 0 { q } else { -q },
                    if o & 2 != 0 { q } else { -q },
                    if o & 4 != 0 { q } else { -q },
                );
            self.nodes.push(Node {
                center: c,
                half: q,
                mass: 0.0,
                com: Vec3::ZERO,
                quad: [0.0; 6],
                open2: 0.0,
                start: offsets[o],
                end: offsets[o] + counts[o],
                first_child: 0,
                n_children: 0,
            });
            children += 1;
        }
        self.nodes[k].first_child = first;
        self.nodes[k].n_children = children;
        for c in first..first + children {
            self.split(c, pos, depth + 1);
        }
    }

    /// Mass, centre of mass, quadrupole and acceptance radius of node `k` and its subtree.
    fn moments(&mut self, k: usize, pos: &[Vec3], mass: &[f64], theta: f64) {
        let node = self.nodes[k];
        let mut m = 0.0;
        let mut weighted = Vec3::ZERO;
        if node.n_children == 0 {
            for &i in &self.order[node.start..node.end] {
                m += mass[i];
                weighted += pos[i] * mass[i];
            }
        } else {
            for c in node.first_child..node.first_child + node.n_children {
                self.moments(c, pos, mass, theta);
                m += self.nodes[c].mass;
                weighted += self.nodes[c].com * self.nodes[c].mass;
            }
        }
        let com = if m > 0.0 { weighted / m } else { node.center };
        // Quadrupole about the centre of mass, from the particles directly (exact; O(N log N)).
        let mut quad = [0.0; 6];
        if m > 0.0 {
            for &i in &self.order[node.start..node.end] {
                let r = pos[i] - com;
                let mi = mass[i];
                let r2 = r.norm_squared();
                quad[0] += mi * (3.0 * r.x * r.x - r2);
                quad[1] += mi * (3.0 * r.y * r.y - r2);
                quad[2] += mi * (3.0 * r.z * r.z - r2);
                quad[3] += mi * 3.0 * r.x * r.y;
                quad[4] += mi * 3.0 * r.x * r.z;
                quad[5] += mi * 3.0 * r.y * r.z;
            }
        }
        let open = 2.0 * node.half / theta;
        let n = &mut self.nodes[k];
        n.mass = m;
        n.com = com;
        n.quad = quad;
        n.open2 = open * open;
    }
}

/// One group's interaction list in structure-of-arrays form, so the inner loops vectorise.
#[derive(Default)]
struct Lists {
    /// Point masses: particles of nearby leaves (including the group's own) and, for
    /// monopole-only runs, accepted nodes.
    px: Vec<f64>,
    py: Vec<f64>,
    pz: Vec<f64>,
    pm: Vec<f64>,
    /// Particle index of each point mass (`usize::MAX` for nodes), to skip self-interaction.
    pid: Vec<usize>,
    /// Accepted nodes' centres of mass and quadrupole moments (xx, yy, zz, xy, xz, yz).
    qx: Vec<f64>,
    qy: Vec<f64>,
    qz: Vec<f64>,
    q: [Vec<f64>; 6],
}

impl Lists {
    fn clear(&mut self) {
        for v in [&mut self.px, &mut self.py, &mut self.pz, &mut self.pm] {
            v.clear();
        }
        self.pid.clear();
        for v in [&mut self.qx, &mut self.qy, &mut self.qz] {
            v.clear();
        }
        for v in &mut self.q {
            v.clear();
        }
    }

    fn push_quad(&mut self, at: Vec3, q: &[f64; 6]) {
        self.qx.push(at.x);
        self.qy.push(at.y);
        self.qz.push(at.z);
        for (v, &c) in self.q.iter_mut().zip(q) {
            v.push(c);
        }
    }

    fn push_point(&mut self, x: Vec3, m: f64, id: usize) {
        self.px.push(x.x);
        self.py.push(x.y);
        self.pz.push(x.z);
        self.pm.push(m);
        self.pid.push(id);
    }
}

const LANES: usize = 4;

/// Sums one particle's interaction list (without G). Point masses are summed in `LANES`
/// independent accumulators so the compiler can vectorise; the particle's own entry gets zero
/// weight (and a unit offset in the denominator so it stays finite). Specialised at compile
/// time on whether the potential is wanted.
fn sum<const PHI: bool>(i: usize, x: Vec3, lists: &Lists, eps2: f64) -> (Vec3, f64) {
    let (mut ax, mut ay, mut az, mut ph) = ([0.0; LANES], [0.0; LANES], [0.0; LANES], [0.0; LANES]);
    let (px, py, pz, pm, pid) = (&lists.px, &lists.py, &lists.pz, &lists.pm, &lists.pid);
    let full = px.len() - px.len() % LANES;
    let mut base = 0;
    while base < full {
        let (cx, cy, cz, cm, ci) = (
            &px[base..base + LANES],
            &py[base..base + LANES],
            &pz[base..base + LANES],
            &pm[base..base + LANES],
            &pid[base..base + LANES],
        );
        let mut is_self = [0.0; LANES];
        for l in 0..LANES {
            is_self[l] = if ci[l] == i { 1.0 } else { 0.0 };
        }
        let mut dx = [0.0; LANES];
        let mut dy = [0.0; LANES];
        let mut dz = [0.0; LANES];
        let mut inv = [0.0; LANES];
        let mut w = [0.0; LANES];
        for l in 0..LANES {
            dx[l] = cx[l] - x.x;
            dy[l] = cy[l] - x.y;
            dz[l] = cz[l] - x.z;
            let s2 = dx[l] * dx[l] + dy[l] * dy[l] + dz[l] * dz[l] + eps2 + is_self[l];
            inv[l] = 1.0 / s2.sqrt();
            w[l] = cm[l] * (1.0 - is_self[l]);
        }
        for l in 0..LANES {
            let k = w[l] * inv[l] * inv[l] * inv[l];
            ax[l] += dx[l] * k;
            ay[l] += dy[l] * k;
            az[l] += dz[l] * k;
            if PHI {
                ph[l] -= w[l] * inv[l];
            }
        }
        base += LANES;
    }
    for j in full..px.len() {
        let is_self = if pid[j] == i { 1.0 } else { 0.0 };
        let (dx, dy, dz) = (px[j] - x.x, py[j] - x.y, pz[j] - x.z);
        let inv = 1.0 / (dx * dx + dy * dy + dz * dz + eps2 + is_self).sqrt();
        let m = pm[j] * (1.0 - is_self);
        let k = m * inv * inv * inv;
        ax[0] += dx * k;
        ay[0] += dy * k;
        az[0] += dz * k;
        if PHI {
            ph[0] -= m * inv;
        }
    }
    let mut a = Vec3::new(ax.iter().sum(), ay.iter().sum(), az.iter().sum());
    let mut phi: f64 = ph.iter().sum();
    // Quadrupole terms, also in lanes.
    let (mut qax, mut qay, mut qaz, mut qph) =
        ([0.0; LANES], [0.0; LANES], [0.0; LANES], [0.0; LANES]);
    let nq = lists.qx.len();
    let mut j = 0;
    while j < nq {
        let width = LANES.min(nq - j);
        for l in 0..width {
            let k = j + l;
            // Field point relative to the centre of mass.
            let (rx, ry, rz) = (x.x - lists.qx[k], x.y - lists.qy[k], x.z - lists.qz[k]);
            let inv = 1.0 / (rx * rx + ry * ry + rz * rz + eps2).sqrt();
            let inv2 = inv * inv;
            let inv5 = inv2 * inv2 * inv;
            let q = &lists.q;
            let qrx = q[0][k] * rx + q[3][k] * ry + q[4][k] * rz;
            let qry = q[3][k] * rx + q[1][k] * ry + q[5][k] * rz;
            let qrz = q[4][k] * rx + q[5][k] * ry + q[2][k] * rz;
            let rqr = rx * qrx + ry * qry + rz * qrz;
            let c = 2.5 * rqr * inv5 * inv2;
            qax[l] += qrx * inv5 - rx * c;
            qay[l] += qry * inv5 - ry * c;
            qaz[l] += qrz * inv5 - rz * c;
            if PHI {
                qph[l] -= 0.5 * rqr * inv5;
            }
        }
        j += LANES;
    }
    a += Vec3::new(qax.iter().sum(), qay.iter().sum(), qaz.iter().sum());
    phi += qph.iter().sum::<f64>();
    (a, phi)
}

/// Newtonian gravity by Barnes-Hut tree: the same physics as [`super::NewtonianGravity`]
/// (including Plummer softening and massless test particles) at O(N log N) cost, with an
/// error controlled by the opening angle `theta` in (0, 1] (0.5 is a common default; smaller is
/// more accurate and slower). `quadrupole` adds the quadrupole term of each node.
#[derive(Debug, Clone)]
pub struct TreeGravity {
    pub g: f64,
    pub softening: f64,
    pub theta: f64,
    pub quadrupole: bool,
}

impl TreeGravity {
    fn check(&self) -> Result<()> {
        // Above ~1.15 a cell could accept the particle inside it (self-interaction).
        if !(self.theta > 0.0 && self.theta <= 1.0) {
            return invalid(format!(
                "TreeGravity: theta must be in (0, 1], got {}",
                self.theta
            ));
        }
        Ok(())
    }

    /// Accelerations (and, if `want_phi`, potentials per unit mass) of every particle.
    ///
    /// Group walk (Barnes 1990): for each leaf, one traversal builds an interaction list
    /// against the leaf's bounding box (accepting a node only if it passes for every point of
    /// the box, so at least as strict as the per-particle test), then every particle of the
    /// leaf sums that list in a tight loop. Leaves are processed in parallel blocks; each
    /// particle's sum has a fixed order, so results do not depend on the thread count.
    fn evaluate(&self, pos: &[Vec3], mass: &[f64], want_phi: bool) -> (Vec<Vec3>, Vec<f64>) {
        let n = pos.len();
        let tree = Tree::build(pos, mass, self.theta);
        // Groups: the largest nodes with at most GROUP_SIZE particles (or leaves).
        let mut leaves: Vec<usize> = Vec::new();
        let mut todo = vec![0usize];
        while let Some(k) = todo.pop() {
            let nk = &tree.nodes[k];
            if nk.n_children == 0 || nk.end - nk.start <= GROUP_SIZE {
                leaves.push(k);
            } else {
                todo.extend(nk.first_child..nk.first_child + nk.n_children);
            }
        }
        leaves.sort_unstable();
        const LEAVES_PER_BLOCK: usize = 32;
        let blocks: Vec<(usize, usize)> = (0..leaves.len())
            .step_by(LEAVES_PER_BLOCK)
            .map(|lo| (lo, (lo + LEAVES_PER_BLOCK).min(leaves.len())))
            .collect();
        let eps2 = self.softening * self.softening;
        let results = parallel::map_blocks(&blocks, |lo, hi| {
            let mut out: Vec<(usize, Vec3, f64)> = Vec::new();
            let mut stack: Vec<usize> = Vec::with_capacity(64);
            let mut lists = Lists::default();
            for &leaf in &leaves[lo..hi] {
                let node = &tree.nodes[leaf];
                let members = &tree.order[node.start..node.end];
                // Bounding box of the group's particles.
                let (mut bmin, mut bmax) = (pos[members[0]], pos[members[0]]);
                for &i in members {
                    let p = pos[i];
                    bmin = Vec3::new(bmin.x.min(p.x), bmin.y.min(p.y), bmin.z.min(p.z));
                    bmax = Vec3::new(bmax.x.max(p.x), bmax.y.max(p.y), bmax.z.max(p.z));
                }
                lists.clear();
                stack.clear();
                stack.push(0);
                while let Some(k) = stack.pop() {
                    let nk = &tree.nodes[k];
                    if nk.mass == 0.0 {
                        continue;
                    }
                    let c = nk.com;
                    let gap = |c: f64, lo: f64, hi: f64| (lo - c).max(c - hi).max(0.0);
                    let g = Vec3::new(
                        gap(c.x, bmin.x, bmax.x),
                        gap(c.y, bmin.y, bmax.y),
                        gap(c.z, bmin.z, bmax.z),
                    );
                    // Never accept a cell that overlaps the group (its own particles would be
                    // folded into the multipole).
                    let overlaps = (0..3).all(|ax| {
                        let (c, h) = (nk.center[ax], nk.half);
                        bmax[ax] >= c - h && bmin[ax] <= c + h
                    });
                    if !overlaps && g.norm_squared() > nk.open2 {
                        lists.push_point(nk.com, nk.mass, usize::MAX);
                        if self.quadrupole {
                            lists.push_quad(nk.com, &nk.quad);
                        }
                    } else if nk.n_children == 0 {
                        for &j in &tree.order[nk.start..nk.end] {
                            if mass[j] != 0.0 {
                                lists.push_point(pos[j], mass[j], j);
                            }
                        }
                    } else {
                        stack.extend(nk.first_child..nk.first_child + nk.n_children);
                    }
                }
                for &i in members {
                    let (a, phi) = if want_phi {
                        sum::<true>(i, pos[i], &lists, eps2)
                    } else {
                        sum::<false>(i, pos[i], &lists, eps2)
                    };
                    out.push((i, a * self.g, phi * self.g));
                }
            }
            out
        });
        let mut acc = vec![Vec3::ZERO; n];
        let mut phi = vec![0.0; if want_phi { n } else { 0 }];
        for block in results {
            for (i, a, p) in block {
                acc[i] = a;
                if want_phi {
                    phi[i] = p;
                }
            }
        }
        (acc, phi)
    }
}

impl Force for TreeGravity {
    fn accumulate(
        &self,
        _t: f64,
        pos: &[Vec3],
        _vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> Result<()> {
        self.check()?;
        if pos.len() < 2 {
            return Ok(());
        }
        let (a, _) = self.evaluate(pos, mass, false);
        for (out, a) in acc.iter_mut().zip(a) {
            *out += a;
        }
        Ok(())
    }

    fn potential(&self, _t: f64, pos: &[Vec3], mass: &[f64]) -> Result<Option<f64>> {
        self.check()?;
        if pos.len() < 2 {
            return Ok(Some(0.0));
        }
        let (_, phi) = self.evaluate(pos, mass, true);
        // Each pair is counted from both ends.
        Ok(Some(
            0.5 * phi.iter().zip(mass).map(|(p, m)| p * m).sum::<f64>(),
        ))
    }

    fn name(&self) -> String {
        "TreeGravity".into()
    }

    fn builtin(&self) -> Option<BuiltinForce> {
        Some(BuiltinForce::TreeGravity(self.clone()))
    }

    fn params(&self) -> Vec<(&'static str, Param)> {
        vec![
            ("G", Param::Scalar(self.g)),
            ("softening", Param::Scalar(self.softening)),
            ("theta", Param::Scalar(self.theta)),
        ]
    }

    fn set_param(&mut self, name: &str, value: Param) -> Result<()> {
        match name {
            "G" => self.g = scalar(name, value)?,
            "softening" => self.softening = non_negative(name, value)?,
            "theta" => {
                let t = scalar(name, value)?;
                if !(t > 0.0 && t <= 1.0) {
                    return invalid(format!("theta must be in (0, 1], got {t}"));
                }
                self.theta = t;
            }
            _ => return invalid(format!("TreeGravity has no parameter {name:?}")),
        }
        Ok(())
    }
}
