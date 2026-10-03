//! Python bindings for rigid bodies (`physim.RigidSystem`).

use numpy::{PyArray1, PyArray2, PyArray3, PyArray4, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use super::{extract_vec3, extract_vecs, flatten, vec3_to_array, vecs_to_array};
use crate::error::{Result as SimResult, SimError};
use crate::rigid::{
    Attachment, BodyForce, BodyGravity, BodySpring, Quat, RigidBody, RigidSystem, RigidTrajectory,
    Wrench,
};
use crate::vec3::Vec3;

/// Two `(N, 3)` arrays.
type Pair<'py> = (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>);

fn extract_quat(obj: &Bound<'_, PyAny>) -> PyResult<Quat> {
    let v: Vec<f64> = super::as_f64_array(obj)?
        .call_method0("ravel")?
        .call_method0("tolist")?
        .extract()?;
    match v[..] {
        [w, x, y, z] => {
            let q = Quat::new(w, x, y, z);
            let n = q.norm();
            if !(n > 0.0 && n.is_finite()) {
                return Err(PyValueError::new_err(
                    "orientation must be a nonzero quaternion",
                ));
            }
            Ok(q.normalized())
        }
        _ => Err(PyValueError::new_err(
            "orientation: expected a quaternion (w, x, y, z)",
        )),
    }
}

fn quats_flat(q: &[Quat]) -> Vec<f64> {
    q.iter().flat_map(|q| q.to_array()).collect()
}

fn matrices_flat(q: &[Quat]) -> Vec<f64> {
    q.iter()
        .flat_map(|q| q.to_matrix().into_iter().flatten())
        .collect()
}

/// Uniform gravity ``g`` acting at each body's centre of mass.
#[pyclass(frozen, name = "BodyGravity", module = "physim")]
struct PyBodyGravity(BodyGravity);

#[pymethods]
impl PyBodyGravity {
    #[new]
    fn new(g: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(BodyGravity {
            g: extract_vec3(g, "g")?,
        }))
    }
    /// The gravitational acceleration ``g``.
    #[getter]
    fn g<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.0.g)
    }
    fn __repr__(&self) -> String {
        let g = self.0.g;
        format!("BodyGravity(g=[{}, {}, {}])", g.x, g.y, g.z)
    }
}

/// ``(body, point)`` with ``point`` in body coordinates, or a fixed 3-vector in space.
fn extract_attachment(obj: &Bound<'_, PyAny>, what: &str) -> PyResult<Attachment> {
    if let Ok(t) = obj.cast::<PyTuple>() {
        if t.len() == 2 {
            if let Ok(body) = t.get_item(0)?.extract::<usize>() {
                return Ok(Attachment::Body {
                    body,
                    point: extract_vec3(&t.get_item(1)?, what)?,
                });
            }
        }
    }
    Ok(Attachment::Fixed(extract_vec3(obj, what).map_err(
        |_| {
            PyValueError::new_err(format!(
                "{what}: expected (body_index, body_point) or a fixed 3-vector"
            ))
        },
    )?))
}

/// Spring ``k (d - rest_length)`` between two attachment points, each either
/// ``(body_index, point_in_body_coordinates)`` or a fixed point in space.
#[pyclass(frozen, name = "BodySpring", module = "physim")]
struct PyBodySpring(BodySpring);

#[pymethods]
impl PyBodySpring {
    #[new]
    #[pyo3(signature = (a, b, k, rest_length = 0.0))]
    fn new(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, k: f64, rest_length: f64) -> PyResult<Self> {
        let s = BodySpring {
            a: extract_attachment(a, "a")?,
            b: extract_attachment(b, "b")?,
            k,
            rest: rest_length,
        };
        s.validate(usize::MAX)?;
        Ok(Self(s))
    }
    fn __repr__(&self) -> String {
        format!("BodySpring(k={}, rest_length={})", self.0.k, self.0.rest)
    }
}

/// A Python callable ``f(t, positions, orientations) -> (forces, torques)``.
struct PythonBodyForce {
    f: Py<PyAny>,
    name: String,
}

impl BodyForce for PythonBodyForce {
    fn wrench(&self, t: f64, bodies: &[RigidBody], out: &mut [Wrench]) -> SimResult<()> {
        Python::attach(|py| -> PyResult<()> {
            let n = bodies.len();
            let pos: Vec<Vec3> = bodies.iter().map(|b| b.pos).collect();
            let q: Vec<Quat> = bodies.iter().map(|b| b.orientation).collect();
            let res = self.f.bind(py).call1((
                t,
                vecs_to_array(py, &pos)?,
                PyArray1::from_vec(py, quats_flat(&q)).reshape([n, 4])?,
            ))?;
            let what = format!("{} return value", self.name);
            let pair = res
                .cast::<PyTuple>()
                .ok()
                .filter(|t| t.len() == 2)
                .ok_or_else(|| {
                    PyValueError::new_err(format!("{what}: expected a (forces, torques) tuple"))
                })?;
            let f = extract_vecs(&pair.get_item(0)?, n, &what)?;
            let tq = extract_vecs(&pair.get_item(1)?, n, &what)?;
            for ((w, f), tq) in out.iter_mut().zip(f).zip(tq) {
                w.force += f;
                w.torque += tq;
            }
            Ok(())
        })
        .map_err(SimError::Python)
    }

    fn name(&self) -> String {
        self.name.clone()
    }
}

/// Recorded frames of a :class:`RigidSystem` run. Arrays are ``(frames, bodies, ...)``.
#[pyclass(frozen, name = "RigidTrajectory", module = "physim")]
struct PyRigidTrajectory {
    /// Recorded times, shape ``(F,)``.
    #[pyo3(get)]
    t: Py<PyArray1<f64>>,
    /// Reference-point positions, shape ``(F, N, 3)``.
    #[pyo3(get)]
    pos: Py<PyArray3<f64>>,
    /// Reference-point velocities, shape ``(F, N, 3)``.
    #[pyo3(get)]
    vel: Py<PyArray3<f64>>,
    /// Unit quaternions ``(w, x, y, z)``, body to space.
    #[pyo3(get)]
    orientation: Py<PyArray3<f64>>,
    /// Rotation matrices, body to space: column ``k`` is body axis ``k``.
    #[pyo3(get)]
    rotation: Py<PyArray4<f64>>,
    /// Angular momentum about each body's reference point, space frame.
    #[pyo3(get)]
    angular_momentum: Py<PyArray3<f64>>,
    /// Angular velocity, space frame.
    #[pyo3(get)]
    angular_velocity: Py<PyArray3<f64>>,
    /// Kinetic energy per frame, shape ``(F,)``.
    #[pyo3(get)]
    kinetic: Py<PyArray1<f64>>,
    /// Potential energy per frame, shape ``(F,)``.
    #[pyo3(get)]
    potential: Py<PyArray1<f64>>,
    /// Total energy per frame, shape ``(F,)``.
    #[pyo3(get)]
    energy: Py<PyArray1<f64>>,
    /// Number of bodies.
    #[pyo3(get)]
    n_bodies: usize,
    /// Number of recorded frames.
    #[pyo3(get)]
    n_frames: usize,
}

#[pymethods]
impl PyRigidTrajectory {
    fn __len__(&self) -> usize {
        self.n_frames
    }
    fn __repr__(&self) -> String {
        format!(
            "RigidTrajectory(frames={}, bodies={})",
            self.n_frames, self.n_bodies
        )
    }
}

impl PyRigidTrajectory {
    fn from_rust(py: Python<'_>, tr: RigidTrajectory) -> PyResult<Self> {
        let (f, n) = (tr.n_frames(), tr.n_bodies);
        let v3 = |v: &[Vec3]| -> PyResult<Py<PyArray3<f64>>> {
            Ok(PyArray1::from_vec(py, flatten(v))
                .reshape([f, n, 3])?
                .unbind())
        };
        let energy = tr.total_energy();
        Ok(Self {
            pos: v3(&tr.pos)?,
            vel: v3(&tr.vel)?,
            angular_momentum: v3(&tr.ang_mom)?,
            angular_velocity: v3(&tr.angular_velocity)?,
            orientation: PyArray1::from_vec(py, quats_flat(&tr.orientation))
                .reshape([f, n, 4])?
                .unbind(),
            rotation: PyArray1::from_vec(py, matrices_flat(&tr.orientation))
                .reshape([f, n, 3, 3])?
                .unbind(),
            t: PyArray1::from_vec(py, tr.t).unbind(),
            kinetic: PyArray1::from_vec(py, tr.kinetic).unbind(),
            potential: PyArray1::from_vec(py, tr.potential).unbind(),
            energy: PyArray1::from_vec(py, energy).unbind(),
            n_bodies: n,
            n_frames: f,
        })
    }
}

/// Rigid bodies with orientation, angular momentum and torques, stepped by a symplectic
/// splitting (free rotation split into exact rotations about the body axes).
///
/// ``order`` is 2 (Strang splitting) or 4 (Yoshida composition).
#[pyclass(name = "RigidSystem", module = "physim", unsendable)]
struct PyRigidSystem {
    inner: RigidSystem,
}

fn opt_vec3(obj: Option<&Bound<'_, PyAny>>, what: &str) -> PyResult<Option<Vec3>> {
    obj.map(|o| extract_vec3(o, what)).transpose()
}

impl PyRigidSystem {
    fn apply_state(
        body: &mut RigidBody,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        orientation: Option<&Bound<'_, PyAny>>,
        angular_velocity: Option<&Bound<'_, PyAny>>,
        angular_momentum: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        if angular_velocity.is_some() && angular_momentum.is_some() {
            return Err(PyValueError::new_err(
                "give angular_velocity or angular_momentum, not both",
            ));
        }
        if let Some(p) = opt_vec3(position, "position")? {
            if body.pivot {
                return Err(PyValueError::new_err(
                    "a pivoted body's position is its pivot; it cannot move",
                ));
            }
            body.pos = p;
        }
        if let Some(v) = opt_vec3(velocity, "velocity")? {
            if body.pivot && v != Vec3::ZERO {
                return Err(PyValueError::new_err("a pivoted body cannot translate"));
            }
            body.vel = v;
        }
        // Keep the angular velocity unless a new one is given.
        let omega = body.angular_velocity();
        if let Some(q) = orientation {
            body.orientation = extract_quat(q)?;
            body.set_angular_velocity(omega);
        }
        if let Some(w) = opt_vec3(angular_velocity, "angular_velocity")? {
            body.set_angular_velocity(w);
        }
        if let Some(l) = opt_vec3(angular_momentum, "angular_momentum")? {
            body.ang_mom = l;
        }
        Ok(())
    }

    fn body(&self, i: usize) -> PyResult<&RigidBody> {
        self.inner
            .bodies()
            .get(i)
            .ok_or_else(|| PyValueError::new_err(format!("body {i} out of range")))
    }
}

#[pymethods]
impl PyRigidSystem {
    #[new]
    #[pyo3(signature = (order = 2))]
    fn new(order: u8) -> PyResult<Self> {
        let mut inner = RigidSystem::new();
        inner.set_order(order)?;
        Ok(Self { inner })
    }

    /// Adds a body and returns its index.
    ///
    /// ``inertia`` holds the principal moments about the centre of mass, along the body
    /// axes. A free body's ``position`` is its centre of mass. With ``pivot`` the body
    /// turns about that fixed point and ``center_of_mass`` (body coordinates, on a
    /// principal axis) places its centre of mass relative to the pivot. Spin is set by
    /// ``angular_velocity`` or ``angular_momentum`` (space frame).
    #[pyo3(signature = (mass, inertia, position = None, velocity = None, orientation = None,
        angular_velocity = None, angular_momentum = None, pivot = None, center_of_mass = None))]
    #[allow(clippy::too_many_arguments)]
    fn add_body(
        &mut self,
        mass: f64,
        inertia: &Bound<'_, PyAny>,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        orientation: Option<&Bound<'_, PyAny>>,
        angular_velocity: Option<&Bound<'_, PyAny>>,
        angular_momentum: Option<&Bound<'_, PyAny>>,
        pivot: Option<&Bound<'_, PyAny>>,
        center_of_mass: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<usize> {
        let inertia = extract_vec3(inertia, "inertia")?;
        let mut body = match opt_vec3(pivot, "pivot")? {
            Some(p) => {
                let com = opt_vec3(center_of_mass, "center_of_mass")?.unwrap_or(Vec3::ZERO);
                RigidBody::pivoted(mass, inertia, p, com)?
            }
            None => {
                if center_of_mass.is_some() {
                    return Err(PyValueError::new_err(
                        "center_of_mass needs a pivot; a free body's position is its centre of mass",
                    ));
                }
                RigidBody::new(mass, inertia)
            }
        };
        Self::apply_state(
            &mut body,
            position,
            velocity,
            orientation,
            angular_velocity,
            angular_momentum,
        )?;
        Ok(self.inner.add_body(body)?)
    }

    /// Changes the state of body ``i``; arguments left out keep their values (a new
    /// ``orientation`` keeps the angular velocity).
    #[pyo3(signature = (i, position = None, velocity = None, orientation = None,
        angular_velocity = None, angular_momentum = None))]
    fn set_body_state(
        &mut self,
        i: usize,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        orientation: Option<&Bound<'_, PyAny>>,
        angular_velocity: Option<&Bound<'_, PyAny>>,
        angular_momentum: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let mut body = self.body(i)?.clone();
        Self::apply_state(
            &mut body,
            position,
            velocity,
            orientation,
            angular_velocity,
            angular_momentum,
        )?;
        Ok(self.inner.set_body(i, body)?)
    }

    /// Adds :class:`BodyGravity`, :class:`BodySpring`, or a callable
    /// ``f(t, positions, orientations) -> (forces, torques)`` returning ``(N, 3)`` arrays
    /// (torques about each body's reference point, space frame).
    fn add_force(&mut self, force: &Bound<'_, PyAny>) -> PyResult<()> {
        let built: Box<dyn BodyForce> = if let Ok(g) = force.extract::<PyRef<'_, PyBodyGravity>>() {
            Box::new(g.0)
        } else if let Ok(s) = force.extract::<PyRef<'_, PyBodySpring>>() {
            Box::new(s.0)
        } else if force.is_callable() {
            let name = force
                .getattr("__name__")
                .and_then(|n| n.extract())
                .unwrap_or_else(|_| "custom body force".to_string());
            Box::new(PythonBodyForce {
                f: force.clone().unbind(),
                name,
            })
        } else {
            return Err(PyValueError::new_err(format!(
                "not a body force: {}",
                force.repr()?
            )));
        };
        Ok(self.inner.add_force(built)?)
    }

    /// Advances by one step ``dt``.
    fn step(&mut self, dt: f64) -> PyResult<()> {
        Ok(self.inner.step(dt)?)
    }

    /// Takes ``steps`` steps of ``dt``, recording the start and every ``record_every``-th step.
    #[pyo3(signature = (dt, steps, record_every = 1))]
    fn run(
        &mut self,
        py: Python<'_>,
        dt: f64,
        steps: usize,
        record_every: usize,
    ) -> PyResult<PyRigidTrajectory> {
        let tr = self.inner.run(dt, steps, record_every)?;
        PyRigidTrajectory::from_rust(py, tr)
    }

    /// Force and torque on each body, as two ``(N, 3)`` arrays.
    fn wrenches<'py>(&self, py: Python<'py>) -> PyResult<Pair<'py>> {
        let w = self.inner.wrenches()?;
        let f: Vec<Vec3> = w.iter().map(|w| w.force).collect();
        let t: Vec<Vec3> = w.iter().map(|w| w.torque).collect();
        Ok((vecs_to_array(py, &f)?, vecs_to_array(py, &t)?))
    }

    /// Translational plus rotational kinetic energy.
    fn kinetic_energy(&self) -> f64 {
        self.inner.kinetic_energy()
    }

    /// Potential energy of the forces that define one (callables count as zero).
    fn potential_energy(&self) -> PyResult<f64> {
        Ok(self.inner.potential_energy()?)
    }

    /// Kinetic plus potential energy.
    fn total_energy(&self) -> PyResult<f64> {
        Ok(self.inner.total_energy()?)
    }

    /// Total linear momentum.
    fn momentum<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.inner.momentum())
    }

    /// Total angular momentum about the origin (orbital plus spin).
    fn angular_momentum<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.inner.angular_momentum())
    }

    /// Current time; assignable.
    #[getter]
    fn t(&self) -> f64 {
        self.inner.t
    }
    #[setter]
    fn set_t(&mut self, t: f64) {
        self.inner.t = t;
    }
    /// Order of the time step (2 or 4).
    #[getter]
    fn order(&self) -> u8 {
        self.inner.order()
    }
    /// Number of bodies.
    #[getter]
    fn n_bodies(&self) -> usize {
        self.inner.bodies().len()
    }
    fn __len__(&self) -> usize {
        self.inner.bodies().len()
    }
    /// Masses, shape ``(N,)``.
    #[getter]
    fn masses<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.bodies().iter().map(|b| b.mass).collect())
    }
    /// Principal moments about each body's reference point (the pivot for pivoted bodies).
    #[getter]
    fn inertia<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let v: Vec<Vec3> = self.inner.bodies().iter().map(|b| b.inertia).collect();
        vecs_to_array(py, &v)
    }
    /// Reference points (centre of mass, or the pivot), shape ``(N, 3)``.
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let v: Vec<Vec3> = self.inner.bodies().iter().map(|b| b.pos).collect();
        vecs_to_array(py, &v)
    }
    /// Velocities of the reference points, shape ``(N, 3)``.
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let v: Vec<Vec3> = self.inner.bodies().iter().map(|b| b.vel).collect();
        vecs_to_array(py, &v)
    }
    /// Centre-of-mass positions (differ from ``positions`` for pivoted bodies).
    #[getter]
    fn centers_of_mass<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let v: Vec<Vec3> = self
            .inner
            .bodies()
            .iter()
            .map(|b| b.com_position())
            .collect();
        vecs_to_array(py, &v)
    }
    /// Body-to-space unit quaternions ``(w, x, y, z)``, shape ``(N, 4)``.
    #[getter]
    fn orientations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let q: Vec<Quat> = self.inner.bodies().iter().map(|b| b.orientation).collect();
        PyArray1::from_vec(py, quats_flat(&q)).reshape([q.len(), 4])
    }
    /// Body-to-space rotation matrices, shape ``(N, 3, 3)``; column ``k`` is body axis ``k``.
    #[getter]
    fn rotation_matrices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f64>>> {
        let q: Vec<Quat> = self.inner.bodies().iter().map(|b| b.orientation).collect();
        PyArray1::from_vec(py, matrices_flat(&q)).reshape([q.len(), 3, 3])
    }
    /// Angular momentum about each reference point (space frame), shape ``(N, 3)``.
    #[getter]
    fn angular_momenta<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let v: Vec<Vec3> = self.inner.bodies().iter().map(|b| b.ang_mom).collect();
        vecs_to_array(py, &v)
    }
    /// Angular velocities (space frame), shape ``(N, 3)``.
    #[getter]
    fn angular_velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let v: Vec<Vec3> = self
            .inner
            .bodies()
            .iter()
            .map(|b| b.angular_velocity())
            .collect();
        vecs_to_array(py, &v)
    }

    fn __repr__(&self) -> String {
        format!(
            "RigidSystem(bodies={}, forces={}, order={}, t={})",
            self.inner.bodies().len(),
            self.inner.forces().len(),
            self.inner.order(),
            self.inner.t
        )
    }
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBodyGravity>()?;
    m.add_class::<PyBodySpring>()?;
    m.add_class::<PyRigidSystem>()?;
    m.add_class::<PyRigidTrajectory>()?;
    Ok(())
}
