//! Rigid bodies: orientation as a unit quaternion, angular momentum, principal moments of
//! inertia, and torques from forces applied at body points.
//!
//! A [`RigidSystem`] holds bodies and the [`BodyForce`]s acting on them, and steps them with
//! a symplectic splitting of `H = T_trans + T_rot + V`:
//!
//! - kick: momenta and angular momenta (space frame) change by the force and torque,
//! - drift: positions move with their velocities, and the free rotation is split into the
//!   three exactly solvable rotations about the body axes, `R1(h/2) R2(h/2) R3(h) R2(h/2)
//!   R1(h/2)` (McLachlan; Dullweber, Leimkuhler & McLachlan),
//!
//! composed as kick(h/2) drift(h) kick(h/2): second order, time reversible, symplectic,
//! and exactly conserving the angular momentum of a torque-free body. `order = 4` composes
//! that step by the Yoshida triple jump.
//!
//! A body is either free (its reference point is its centre of mass) or pivoted about a
//! fixed reference point (a top or a compound pendulum), in which case `inertia` is about
//! the pivot and `com` is the centre of mass in body coordinates.

use crate::error::{invalid, Result};
use crate::vec3::Vec3;
use std::ops::Mul;

/// A quaternion `w + x i + y j + z k`; unit quaternions represent rotations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quat {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Default for Quat {
    fn default() -> Self {
        Quat::IDENTITY
    }
}

impl Quat {
    pub const IDENTITY: Quat = Quat::new(1.0, 0.0, 0.0, 0.0);

    pub const fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self { w, x, y, z }
    }

    /// Rotation by `angle` (radians, right-handed) about `axis` (any nonzero length).
    pub fn from_axis_angle(axis: Vec3, angle: f64) -> Quat {
        let n = axis.norm();
        if n == 0.0 {
            return Quat::IDENTITY;
        }
        let (s, c) = (0.5 * angle).sin_cos();
        let a = axis * (s / n);
        Quat::new(c, a.x, a.y, a.z)
    }

    pub fn vector(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    pub fn conjugate(self) -> Quat {
        Quat::new(self.w, -self.x, -self.y, -self.z)
    }

    pub fn norm(self) -> f64 {
        (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalized(self) -> Quat {
        let n = self.norm();
        Quat::new(self.w / n, self.x / n, self.y / n, self.z / n)
    }

    /// Rotates `v` (body to space for a body orientation).
    pub fn rotate(self, v: Vec3) -> Vec3 {
        let u = self.vector();
        let t = 2.0 * u.cross(v);
        v + self.w * t + u.cross(t)
    }

    /// Rotates `v` by the inverse rotation (space to body).
    pub fn rotate_inverse(self, v: Vec3) -> Vec3 {
        self.conjugate().rotate(v)
    }

    /// The rotation matrix, row major.
    pub fn to_matrix(self) -> [[f64; 3]; 3] {
        let ex = self.rotate(Vec3::new(1.0, 0.0, 0.0));
        let ey = self.rotate(Vec3::new(0.0, 1.0, 0.0));
        let ez = self.rotate(Vec3::new(0.0, 0.0, 1.0));
        [[ex.x, ey.x, ez.x], [ex.y, ey.y, ez.y], [ex.z, ey.z, ez.z]]
    }

    pub fn to_array(self) -> [f64; 4] {
        [self.w, self.x, self.y, self.z]
    }
}

impl Mul for Quat {
    type Output = Quat;
    fn mul(self, o: Quat) -> Quat {
        Quat::new(
            self.w * o.w - self.x * o.x - self.y * o.y - self.z * o.z,
            self.w * o.x + self.x * o.w + self.y * o.z - self.z * o.y,
            self.w * o.y - self.x * o.z + self.y * o.w + self.z * o.x,
            self.w * o.z + self.x * o.y - self.y * o.x + self.z * o.w,
        )
    }
}

/// Principal moments of inertia of common solids about their centre, along body x, y, z.
pub mod inertia {
    use crate::vec3::Vec3;

    /// Box with side lengths `a`, `b`, `c` along x, y, z.
    pub fn solid_box(mass: f64, a: f64, b: f64, c: f64) -> Vec3 {
        Vec3::new(b * b + c * c, a * a + c * c, a * a + b * b) * (mass / 12.0)
    }

    /// Ellipsoid with semi-axes `a`, `b`, `c` along x, y, z.
    pub fn solid_ellipsoid(mass: f64, a: f64, b: f64, c: f64) -> Vec3 {
        Vec3::new(b * b + c * c, a * a + c * c, a * a + b * b) * (mass / 5.0)
    }

    pub fn solid_sphere(mass: f64, radius: f64) -> Vec3 {
        let i = 0.4 * mass * radius * radius;
        Vec3::new(i, i, i)
    }

    /// Cylinder of `radius` and `length` with its axis along z.
    pub fn solid_cylinder(mass: f64, radius: f64, length: f64) -> Vec3 {
        let t = mass * (3.0 * radius * radius + length * length) / 12.0;
        Vec3::new(t, t, 0.5 * mass * radius * radius)
    }
}

/// One rigid body. Vectors are in the space frame unless they say otherwise.
#[derive(Debug, Clone, PartialEq)]
pub struct RigidBody {
    pub mass: f64,
    /// Principal moments of inertia about the reference point, along the body axes.
    pub inertia: Vec3,
    /// Reference point: the centre of mass, or the pivot of a pivoted body.
    pub pos: Vec3,
    /// Velocity of the reference point (zero for a pivoted body).
    pub vel: Vec3,
    /// Body-to-space rotation.
    pub orientation: Quat,
    /// Angular momentum about the reference point.
    pub ang_mom: Vec3,
    /// The reference point is a fixed pivot.
    pub pivot: bool,
    /// Centre of mass in body coordinates relative to the reference point (zero unless pivoted).
    pub com: Vec3,
}

impl RigidBody {
    /// A free body at rest at the origin with the given principal moments about its centre.
    pub fn new(mass: f64, inertia: Vec3) -> Self {
        Self {
            mass,
            inertia,
            pos: Vec3::ZERO,
            vel: Vec3::ZERO,
            orientation: Quat::IDENTITY,
            ang_mom: Vec3::ZERO,
            pivot: false,
            com: Vec3::ZERO,
        }
    }

    /// A body pivoted at `pivot`, its centre of mass at `com` in body coordinates relative
    /// to the pivot, with `inertia_about_com` its principal moments about the centre of mass.
    /// `com` must lie on a principal axis so the moments about the pivot stay principal
    /// (parallel-axis theorem).
    pub fn pivoted(mass: f64, inertia_about_com: Vec3, pivot: Vec3, com: Vec3) -> Result<Self> {
        let off_axis = [com.x, com.y, com.z].iter().filter(|c| **c != 0.0).count();
        if off_axis > 1 {
            return invalid(format!(
                "the centre of mass of a pivoted body must lie on a principal axis, got {com:?}"
            ));
        }
        let d2 = com.norm_squared();
        let shift = Vec3::new(d2 - com.x * com.x, d2 - com.y * com.y, d2 - com.z * com.z) * mass;
        Ok(Self {
            inertia: inertia_about_com + shift,
            pos: pivot,
            pivot: true,
            com,
            ..Self::new(mass, Vec3::ZERO)
        })
    }

    /// Converts a body-frame vector to the space frame.
    pub fn to_space(&self, v: Vec3) -> Vec3 {
        self.orientation.rotate(v)
    }

    /// Converts a space-frame vector to the body frame.
    pub fn to_body(&self, v: Vec3) -> Vec3 {
        self.orientation.rotate_inverse(v)
    }

    /// Angular velocity in the space frame, `R I⁻¹ Rᵀ L`.
    pub fn angular_velocity(&self) -> Vec3 {
        let l = self.to_body(self.ang_mom);
        self.to_space(Vec3::new(
            l.x / self.inertia.x,
            l.y / self.inertia.y,
            l.z / self.inertia.z,
        ))
    }

    /// Sets the angular momentum that gives the space-frame angular velocity `omega`.
    pub fn set_angular_velocity(&mut self, omega: Vec3) {
        let w = self.to_body(omega);
        self.ang_mom = self.to_space(Vec3::new(
            w.x * self.inertia.x,
            w.y * self.inertia.y,
            w.z * self.inertia.z,
        ));
    }

    /// Position of a body point given in body coordinates relative to the reference point.
    pub fn point(&self, body_point: Vec3) -> Vec3 {
        self.pos + self.to_space(body_point)
    }

    /// Velocity of a body point given in body coordinates.
    pub fn point_velocity(&self, body_point: Vec3) -> Vec3 {
        self.vel + self.angular_velocity().cross(self.to_space(body_point))
    }

    pub fn com_position(&self) -> Vec3 {
        self.point(self.com)
    }

    pub fn com_velocity(&self) -> Vec3 {
        self.point_velocity(self.com)
    }

    /// Rotational kinetic energy `½ L·ω` about the reference point.
    pub fn rotational_energy(&self) -> f64 {
        let l = self.to_body(self.ang_mom);
        0.5 * (l.x * l.x / self.inertia.x + l.y * l.y / self.inertia.y + l.z * l.z / self.inertia.z)
    }

    /// Total kinetic energy (translation plus rotation).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * self.vel.norm_squared() + self.rotational_energy()
    }

    fn validate(&self, index: usize) -> Result<()> {
        let i = self.inertia;
        let finite = [i.x, i.y, i.z, self.mass].iter().all(|v| v.is_finite());
        if !(finite && self.mass > 0.0 && i.x > 0.0 && i.y > 0.0 && i.z > 0.0) {
            return invalid(format!(
                "body {index}: mass and principal moments must be positive, got {} and {i:?}",
                self.mass
            ));
        }
        // Moments of a real mass distribution satisfy the triangle inequality.
        let tol = 1e-12 * (i.x + i.y + i.z);
        if i.x + i.y < i.z - tol || i.y + i.z < i.x - tol || i.x + i.z < i.y - tol {
            return invalid(format!(
                "body {index}: principal moments {i:?} violate the triangle inequality"
            ));
        }
        if !self.pivot && self.com != Vec3::ZERO {
            return invalid(format!(
                "body {index}: a free body's reference point is its centre of mass (com must be zero)"
            ));
        }
        let n = self.orientation.norm();
        if (n - 1.0).abs() >= 1e-6 || n.is_nan() {
            return invalid(format!(
                "body {index}: orientation must be a unit quaternion, |q| = {n}"
            ));
        }
        Ok(())
    }
}

/// Force and torque (about the reference point) on one body.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Wrench {
    pub force: Vec3,
    pub torque: Vec3,
}

impl Wrench {
    /// Adds a force `f` applied at `lever`: the space-frame offset of the point of application
    /// from the reference point.
    pub fn add_at(&mut self, f: Vec3, lever: Vec3) {
        self.force += f;
        self.torque += lever.cross(f);
    }
}

/// Something that pushes and twists rigid bodies.
pub trait BodyForce: Send + Sync {
    /// Adds the force and torque on each body to `out`.
    fn wrench(&self, t: f64, bodies: &[RigidBody], out: &mut [Wrench]) -> Result<()>;

    /// Potential energy, if the force is conservative.
    fn potential(&self, _t: f64, _bodies: &[RigidBody]) -> Result<Option<f64>> {
        Ok(None)
    }

    /// Checks body indices against a system of `n` bodies.
    fn validate(&self, _n: usize) -> Result<()> {
        Ok(())
    }

    fn name(&self) -> String;
}

/// Uniform gravity `g`, acting at each body's centre of mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyGravity {
    pub g: Vec3,
}

impl BodyForce for BodyGravity {
    fn wrench(&self, _t: f64, bodies: &[RigidBody], out: &mut [Wrench]) -> Result<()> {
        for (b, w) in bodies.iter().zip(out) {
            w.add_at(self.g * b.mass, b.to_space(b.com));
        }
        Ok(())
    }

    fn potential(&self, _t: f64, bodies: &[RigidBody]) -> Result<Option<f64>> {
        Ok(Some(
            bodies
                .iter()
                .map(|b| -b.mass * self.g.dot(b.com_position()))
                .sum(),
        ))
    }

    fn name(&self) -> String {
        "BodyGravity".into()
    }
}

/// End of a [`BodySpring`]: a point fixed to a body (body coordinates relative to its
/// reference point) or a fixed point in space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Attachment {
    Body { body: usize, point: Vec3 },
    Fixed(Vec3),
}

impl Attachment {
    fn position(&self, bodies: &[RigidBody]) -> Vec3 {
        match *self {
            Attachment::Body { body, point } => bodies[body].point(point),
            Attachment::Fixed(p) => p,
        }
    }

    fn apply(&self, bodies: &[RigidBody], out: &mut [Wrench], f: Vec3) {
        if let Attachment::Body { body, point } = *self {
            out[body].add_at(f, bodies[body].to_space(point));
        }
    }
}

/// Spring `k (d - rest)` between two attachment points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodySpring {
    pub a: Attachment,
    pub b: Attachment,
    pub k: f64,
    pub rest: f64,
}

impl BodyForce for BodySpring {
    fn wrench(&self, _t: f64, bodies: &[RigidBody], out: &mut [Wrench]) -> Result<()> {
        let d = self.b.position(bodies) - self.a.position(bodies);
        let len = d.norm();
        let f = if len > 0.0 {
            d * (self.k * (len - self.rest) / len)
        } else if self.rest == 0.0 {
            Vec3::ZERO
        } else {
            return invalid("BodySpring: the ends coincide, so the direction is undefined");
        };
        self.a.apply(bodies, out, f);
        self.b.apply(bodies, out, -f);
        Ok(())
    }

    fn potential(&self, _t: f64, bodies: &[RigidBody]) -> Result<Option<f64>> {
        let len = (self.b.position(bodies) - self.a.position(bodies)).norm();
        Ok(Some(0.5 * self.k * (len - self.rest).powi(2)))
    }

    fn validate(&self, n: usize) -> Result<()> {
        for end in [self.a, self.b] {
            if let Attachment::Body { body, .. } = end {
                if body >= n {
                    return invalid(format!("BodySpring: body {body} out of range ({n} bodies)"));
                }
            }
        }
        if !(self.k.is_finite() && self.rest >= 0.0) {
            return invalid("BodySpring: k must be finite and rest length non-negative");
        }
        Ok(())
    }

    fn name(&self) -> String {
        "BodySpring".into()
    }
}

/// A body force from a closure (no potential).
pub struct BodyClosure<F>
where
    F: Fn(f64, &[RigidBody], &mut [Wrench]) -> Result<()> + Send + Sync,
{
    pub f: F,
    pub name: String,
}

impl<F> BodyForce for BodyClosure<F>
where
    F: Fn(f64, &[RigidBody], &mut [Wrench]) -> Result<()> + Send + Sync,
{
    fn wrench(&self, t: f64, bodies: &[RigidBody], out: &mut [Wrench]) -> Result<()> {
        (self.f)(t, bodies, out)
    }

    fn name(&self) -> String {
        self.name.clone()
    }
}

/// Recorded frames of a [`RigidSystem`] run.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RigidTrajectory {
    pub n_bodies: usize,
    pub t: Vec<f64>,
    /// Frame-major: frame `f`, body `i` at `f * n_bodies + i`.
    pub pos: Vec<Vec3>,
    pub vel: Vec<Vec3>,
    pub orientation: Vec<Quat>,
    pub ang_mom: Vec<Vec3>,
    pub angular_velocity: Vec<Vec3>,
    pub kinetic: Vec<f64>,
    /// Potential energy of the forces that have one (others count as zero).
    pub potential: Vec<f64>,
}

impl RigidTrajectory {
    pub fn n_frames(&self) -> usize {
        self.t.len()
    }

    pub fn total_energy(&self) -> Vec<f64> {
        self.kinetic
            .iter()
            .zip(&self.potential)
            .map(|(k, u)| k + u)
            .collect()
    }
}

/// Rigid bodies and the forces between them.
pub struct RigidSystem {
    pub t: f64,
    bodies: Vec<RigidBody>,
    forces: Vec<Box<dyn BodyForce>>,
    order: u8,
}

impl Default for RigidSystem {
    fn default() -> Self {
        Self::new()
    }
}

impl RigidSystem {
    pub fn new() -> Self {
        Self {
            t: 0.0,
            bodies: Vec::new(),
            forces: Vec::new(),
            order: 2,
        }
    }

    /// Order of the time step: 2 (Strang splitting) or 4 (Yoshida composition).
    pub fn set_order(&mut self, order: u8) -> Result<()> {
        if order != 2 && order != 4 {
            return invalid(format!("order must be 2 or 4, got {order}"));
        }
        self.order = order;
        Ok(())
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn add_body(&mut self, mut body: RigidBody) -> Result<usize> {
        body.validate(self.bodies.len())?;
        body.orientation = body.orientation.normalized();
        if body.pivot {
            body.vel = Vec3::ZERO;
        }
        self.bodies.push(body);
        Ok(self.bodies.len() - 1)
    }

    pub fn bodies(&self) -> &[RigidBody] {
        &self.bodies
    }

    /// Replaces body `i`.
    pub fn set_body(&mut self, i: usize, body: RigidBody) -> Result<()> {
        if i >= self.bodies.len() {
            return invalid(format!("body {i} out of range"));
        }
        body.validate(i)?;
        self.bodies[i] = body;
        Ok(())
    }

    pub fn add_force(&mut self, force: Box<dyn BodyForce>) -> Result<()> {
        force.validate(self.bodies.len())?;
        self.forces.push(force);
        Ok(())
    }

    pub fn forces(&self) -> &[Box<dyn BodyForce>] {
        &self.forces
    }

    pub fn wrenches(&self) -> Result<Vec<Wrench>> {
        let mut out = vec![Wrench::default(); self.bodies.len()];
        for f in &self.forces {
            f.validate(self.bodies.len())?;
            f.wrench(self.t, &self.bodies, &mut out)?;
        }
        Ok(out)
    }

    pub fn kinetic_energy(&self) -> f64 {
        self.bodies.iter().map(RigidBody::kinetic_energy).sum()
    }

    /// Sum of the potentials of the forces that define one.
    pub fn potential_energy(&self) -> Result<f64> {
        let mut u = 0.0;
        for f in &self.forces {
            u += f.potential(self.t, &self.bodies)?.unwrap_or(0.0);
        }
        Ok(u)
    }

    pub fn total_energy(&self) -> Result<f64> {
        Ok(self.kinetic_energy() + self.potential_energy()?)
    }

    /// Total linear momentum (pivoted bodies count with their centre-of-mass motion).
    pub fn momentum(&self) -> Vec3 {
        self.bodies
            .iter()
            .fold(Vec3::ZERO, |p, b| p + b.com_velocity() * b.mass)
    }

    /// Total angular momentum about the origin: orbital plus spin.
    pub fn angular_momentum(&self) -> Vec3 {
        self.bodies.iter().fold(Vec3::ZERO, |l, b| {
            // About the reference point plus the reference point's own motion.
            l + b.ang_mom + b.pos.cross(b.com_velocity() * b.mass)
        })
    }

    fn kick(&mut self, h: f64) -> Result<()> {
        let w = self.wrenches()?;
        for (b, w) in self.bodies.iter_mut().zip(w) {
            if !b.pivot {
                b.vel += w.force * (h / b.mass);
            }
            b.ang_mom += w.torque * h;
        }
        Ok(())
    }

    fn drift(&mut self, h: f64) {
        for b in &mut self.bodies {
            if !b.pivot {
                b.pos += b.vel * h;
            }
            for (axis, frac) in [(0, 0.5), (1, 0.5), (2, 1.0), (1, 0.5), (0, 0.5)] {
                rotate_about_body_axis(b, axis, frac * h);
            }
            b.orientation = b.orientation.normalized();
        }
        self.t += h;
    }

    fn strang(&mut self, h: f64) -> Result<()> {
        self.kick(0.5 * h)?;
        self.drift(h);
        self.kick(0.5 * h)
    }

    /// Advances by one step `dt`.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        if !(dt.is_finite() && dt != 0.0) {
            return invalid(format!("dt must be finite and nonzero, got {dt}"));
        }
        if self.order == 4 {
            let w1 = 1.0 / (2.0 - 2f64.powf(1.0 / 3.0));
            let w0 = 1.0 - 2.0 * w1;
            for w in [w1, w0, w1] {
                self.strang(w * dt)?;
            }
            Ok(())
        } else {
            self.strang(dt)
        }
    }

    fn record(&self, tr: &mut RigidTrajectory) -> Result<()> {
        tr.t.push(self.t);
        for b in &self.bodies {
            tr.pos.push(b.pos);
            tr.vel.push(b.vel);
            tr.orientation.push(b.orientation);
            tr.ang_mom.push(b.ang_mom);
            tr.angular_velocity.push(b.angular_velocity());
        }
        tr.kinetic.push(self.kinetic_energy());
        tr.potential.push(self.potential_energy()?);
        Ok(())
    }

    /// Takes `steps` steps of `dt`, recording the start and every `record_every`-th step.
    pub fn run(&mut self, dt: f64, steps: usize, record_every: usize) -> Result<RigidTrajectory> {
        if record_every == 0 {
            return invalid("record_every must be at least 1");
        }
        let mut tr = RigidTrajectory {
            n_bodies: self.bodies.len(),
            ..Default::default()
        };
        self.record(&mut tr)?;
        for k in 1..=steps {
            self.step(dt)?;
            if k % record_every == 0 {
                self.record(&mut tr)?;
            }
        }
        Ok(tr)
    }
}

/// Exact flow of `H_k = L_k² / (2 I_k)`: the body turns about its own axis `k` at the
/// constant rate `L_k / I_k`, and the space-frame angular momentum does not change.
fn rotate_about_body_axis(b: &mut RigidBody, axis: usize, h: f64) {
    let l = b.to_body(b.ang_mom);
    let rate = l[axis] / b.inertia[axis];
    let mut e = [0.0; 3];
    e[axis] = 1.0;
    b.orientation = b.orientation * Quat::from_axis_angle(Vec3::from(e), rate * h);
}
