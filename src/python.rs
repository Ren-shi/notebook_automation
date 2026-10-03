//! Python bindings (`physim._core`). Built only with the `python` feature.

use numpy::ndarray::ArrayView2;
use numpy::{
    PyArray1, PyArray2, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::error::{Result as SimResult, SimError};
use crate::events::{self, Direction, Event, EventFunction};
use crate::forces::{self, Force, ForceId, Param};
use crate::integrators;
use crate::vec3::Vec3;
use crate::world::{RunFailure, Trajectory, World};

impl From<SimError> for PyErr {
    fn from(e: SimError) -> PyErr {
        match e {
            SimError::Invalid(msg) => PyValueError::new_err(msg),
            SimError::Python(err) => err,
        }
    }
}

// ---------------------------------------------------------------------------
// numpy <-> Vec3 conversion

fn flatten(v: &[Vec3]) -> Vec<f64> {
    v.iter().flat_map(|p| p.to_array()).collect()
}

fn vecs_to_array<'py>(py: Python<'py>, v: &[Vec3]) -> PyResult<Bound<'py, PyArray2<f64>>> {
    PyArray1::from_vec(py, flatten(v)).reshape([v.len(), 3])
}

fn view_to_vecs(a: ArrayView2<'_, f64>, n: usize, what: &str) -> PyResult<Vec<Vec3>> {
    if a.shape() != [n, 3] {
        return Err(PyValueError::new_err(format!(
            "{what}: expected shape ({n}, 3), got {:?}",
            a.shape()
        )));
    }
    Ok(a.rows()
        .into_iter()
        .map(|r| Vec3::new(r[0], r[1], r[2]))
        .collect())
}

/// Accepts anything numpy can turn into a float64 array (lists, tuples, int arrays...).
fn as_f64_array<'py>(obj: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
    obj.py()
        .import("numpy")?
        .call_method1("asarray", (obj, "float64"))
}

fn extract_vecs(obj: &Bound<'_, PyAny>, n: usize, what: &str) -> PyResult<Vec<Vec3>> {
    let arr = as_f64_array(obj)?;
    let arr: PyReadonlyArray2<'_, f64> = arr.extract().map_err(|_| {
        PyValueError::new_err(format!("{what}: expected a 2D array of shape ({n}, 3)"))
    })?;
    view_to_vecs(arr.as_array(), n, what)
}

fn extract_vec3(obj: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec3> {
    let arr = as_f64_array(obj)?;
    let arr: PyReadonlyArray1<'_, f64> = arr
        .extract()
        .map_err(|_| PyValueError::new_err(format!("{what}: expected a 3-vector")))?;
    match arr.as_slice()? {
        [x, y, z] => Ok(Vec3::new(*x, *y, *z)),
        s => Err(PyValueError::new_err(format!(
            "{what}: expected 3 components, got {}",
            s.len()
        ))),
    }
}

fn vec3_to_array<'py>(py: Python<'py>, v: Vec3) -> Bound<'py, PyArray1<f64>> {
    PyArray1::from_slice(py, &v.to_array())
}

// ---------------------------------------------------------------------------
// Forces

/// Builds a fresh Rust force from a Python force object.
trait PyForceSpec {
    fn build(&self, py: Python<'_>) -> Box<dyn Force>;
}

/// Constant acceleration field, e.g. ``UniformField([0, -9.81, 0])``.
#[pyclass(frozen, name = "UniformField", module = "physim")]
struct PyUniformField(forces::UniformField);

#[pymethods]
impl PyUniformField {
    #[new]
    fn new(g: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(forces::UniformField {
            g: extract_vec3(g, "g")?,
        }))
    }
    #[getter]
    fn g<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.0.g)
    }
    fn __repr__(&self) -> String {
        let g = self.0.g;
        format!("UniformField(g=[{}, {}, {}])", g.x, g.y, g.z)
    }
}

impl PyForceSpec for PyUniformField {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

/// Pairwise Newtonian gravity with Plummer softening.
#[pyclass(frozen, name = "NewtonianGravity", module = "physim")]
struct PyNewtonianGravity(forces::NewtonianGravity);

#[pymethods]
impl PyNewtonianGravity {
    #[new]
    #[pyo3(signature = (G = 1.0, softening = 0.0))]
    #[allow(non_snake_case)]
    fn new(G: f64, softening: f64) -> Self {
        Self(forces::NewtonianGravity { g: G, softening })
    }
    #[getter(G)]
    fn g(&self) -> f64 {
        self.0.g
    }
    #[getter]
    fn softening(&self) -> f64 {
        self.0.softening
    }
    fn __repr__(&self) -> String {
        format!(
            "NewtonianGravity(G={}, softening={})",
            self.0.g, self.0.softening
        )
    }
}

impl PyForceSpec for PyNewtonianGravity {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

/// Hookean spring between particles ``i`` and ``j``.
#[pyclass(frozen, name = "Spring", module = "physim")]
struct PySpring(forces::Spring);

#[pymethods]
impl PySpring {
    #[new]
    #[pyo3(signature = (i, j, k, rest_length))]
    fn new(i: usize, j: usize, k: f64, rest_length: f64) -> Self {
        Self(forces::Spring {
            i,
            j,
            k,
            rest_length,
        })
    }
    fn __repr__(&self) -> String {
        let s = &self.0;
        format!(
            "Spring(i={}, j={}, k={}, rest_length={})",
            s.i, s.j, s.k, s.rest_length
        )
    }
}

impl PyForceSpec for PySpring {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

/// Spring tying particle ``i`` to a fixed ``anchor`` point.
#[pyclass(frozen, name = "AnchorSpring", module = "physim")]
struct PyAnchorSpring(forces::AnchorSpring);

#[pymethods]
impl PyAnchorSpring {
    #[new]
    #[pyo3(signature = (i, anchor, k, rest_length = 0.0))]
    fn new(i: usize, anchor: &Bound<'_, PyAny>, k: f64, rest_length: f64) -> PyResult<Self> {
        Ok(Self(forces::AnchorSpring {
            i,
            anchor: extract_vec3(anchor, "anchor")?,
            k,
            rest_length,
        }))
    }
    fn __repr__(&self) -> String {
        let s = &self.0;
        let a = s.anchor;
        format!(
            "AnchorSpring(i={}, anchor=[{}, {}, {}], k={}, rest_length={})",
            s.i, a.x, a.y, a.z, s.k, s.rest_length
        )
    }
}

impl PyForceSpec for PyAnchorSpring {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

/// Linear drag ``a = -gamma v``.
#[pyclass(frozen, name = "LinearDrag", module = "physim")]
struct PyLinearDrag(forces::LinearDrag);

#[pymethods]
impl PyLinearDrag {
    #[new]
    fn new(gamma: f64) -> Self {
        Self(forces::LinearDrag { gamma })
    }
    fn __repr__(&self) -> String {
        format!("LinearDrag(gamma={})", self.0.gamma)
    }
}

impl PyForceSpec for PyLinearDrag {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

/// Quadratic drag force ``F = -c |v| v``.
#[pyclass(frozen, name = "QuadraticDrag", module = "physim")]
struct PyQuadraticDrag(forces::QuadraticDrag);

#[pymethods]
impl PyQuadraticDrag {
    #[new]
    fn new(c: f64) -> Self {
        Self(forces::QuadraticDrag { c })
    }
    fn __repr__(&self) -> String {
        format!("QuadraticDrag(c={})", self.0.c)
    }
}

impl PyForceSpec for PyQuadraticDrag {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

/// A force defined in Python, for prototyping new physics without recompiling.
///
/// ``acceleration(t, pos, vel, mass)`` receives ``pos``/``vel`` as ``(N, 3)`` arrays and
/// ``mass`` as ``(N,)``, and must return the ``(N, 3)`` accelerations. The optional
/// ``potential(t, pos, mass)`` returns the potential energy so it is included in energy
/// diagnostics. Set ``velocity_dependent=False`` when ``acceleration`` ignores ``vel``.
#[pyclass(frozen, name = "CustomForce", module = "physim")]
struct PyCustomForce {
    acceleration: Py<PyAny>,
    potential: Option<Py<PyAny>>,
    velocity_dependent: bool,
    name: String,
}

#[pymethods]
impl PyCustomForce {
    #[new]
    #[pyo3(signature = (acceleration, potential = None, velocity_dependent = true, name = None))]
    fn new(
        acceleration: &Bound<'_, PyAny>,
        potential: Option<&Bound<'_, PyAny>>,
        velocity_dependent: bool,
        name: Option<String>,
    ) -> PyResult<Self> {
        if !acceleration.is_callable() || potential.is_some_and(|p| !p.is_callable()) {
            return Err(PyValueError::new_err(
                "acceleration and potential must be callable",
            ));
        }
        let name = match name {
            Some(n) => n,
            None => acceleration
                .getattr("__name__")
                .and_then(|n| n.extract())
                .unwrap_or_else(|_| "CustomForce".to_string()),
        };
        Ok(Self {
            acceleration: acceleration.clone().unbind(),
            potential: potential.map(|p| p.clone().unbind()),
            velocity_dependent,
            name,
        })
    }
    fn __repr__(&self) -> String {
        format!("CustomForce({})", self.name)
    }
}

impl PyForceSpec for PyCustomForce {
    fn build(&self, py: Python<'_>) -> Box<dyn Force> {
        Box::new(PythonForce {
            acceleration: self.acceleration.clone_ref(py),
            potential: self.potential.as_ref().map(|p| p.clone_ref(py)),
            velocity_dependent: self.velocity_dependent,
            name: self.name.clone(),
        })
    }
}

struct PythonForce {
    acceleration: Py<PyAny>,
    potential: Option<Py<PyAny>>,
    velocity_dependent: bool,
    name: String,
}

impl Force for PythonForce {
    fn accumulate(
        &self,
        t: f64,
        pos: &[Vec3],
        vel: &[Vec3],
        mass: &[f64],
        acc: &mut [Vec3],
    ) -> SimResult<()> {
        Python::attach(|py| -> PyResult<()> {
            let out = self.acceleration.bind(py).call1((
                t,
                vecs_to_array(py, pos)?,
                vecs_to_array(py, vel)?,
                PyArray1::from_slice(py, mass),
            ))?;
            let what = format!("{} return value", self.name);
            for (a, extra) in acc.iter_mut().zip(extract_vecs(&out, pos.len(), &what)?) {
                *a += extra;
            }
            Ok(())
        })
        .map_err(SimError::Python)
    }

    fn potential(&self, t: f64, pos: &[Vec3], mass: &[f64]) -> SimResult<Option<f64>> {
        let Some(potential) = &self.potential else {
            return Ok(None);
        };
        Python::attach(|py| -> PyResult<Option<f64>> {
            let u = potential.bind(py).call1((
                t,
                vecs_to_array(py, pos)?,
                PyArray1::from_slice(py, mass),
            ))?;
            Ok(Some(u.extract()?))
        })
        .map_err(SimError::Python)
    }

    fn velocity_dependent(&self) -> bool {
        self.velocity_dependent
    }

    fn name(&self) -> String {
        self.name.clone()
    }
}

fn build_force(force: &Bound<'_, PyAny>) -> PyResult<Box<dyn Force>> {
    let py = force.py();
    macro_rules! try_build {
        ($($ty:ty),*) => {
            $(if let Ok(f) = force.extract::<PyRef<'_, $ty>>() {
                return Ok(f.build(py));
            })*
        };
    }
    try_build!(
        PyUniformField,
        PyNewtonianGravity,
        PySpring,
        PyAnchorSpring,
        PyLinearDrag,
        PyQuadraticDrag,
        PyCustomForce
    );
    Err(PyValueError::new_err(format!(
        "not a physim force: {}; wrap Python functions in physim.CustomForce",
        force.repr()?
    )))
}

// ---------------------------------------------------------------------------
// Trajectory

/// Recorded frames of a run: ``t`` (F,), ``pos``/``vel`` (F, N, 3),
/// ``kinetic``/``potential``/``energy`` (F,).
#[pyclass(frozen, name = "Trajectory", module = "physim")]
struct PyTrajectory {
    #[pyo3(get)]
    t: Py<PyArray1<f64>>,
    #[pyo3(get)]
    pos: Py<PyArray3<f64>>,
    #[pyo3(get)]
    vel: Py<PyArray3<f64>>,
    #[pyo3(get)]
    kinetic: Py<PyArray1<f64>>,
    #[pyo3(get)]
    potential: Py<PyArray1<f64>>,
    #[pyo3(get)]
    energy: Py<PyArray1<f64>>,
    #[pyo3(get)]
    n_particles: usize,
    /// Time of each detected event, (K,).
    #[pyo3(get)]
    event_t: Py<PyArray1<f64>>,
    /// Which event (index into the ``events`` list) fired, (K,).
    #[pyo3(get)]
    event_index: Py<PyArray1<i64>>,
    /// Positions and velocities at each event, (K, N, 3).
    #[pyo3(get)]
    event_pos: Py<PyArray3<f64>>,
    #[pyo3(get)]
    event_vel: Py<PyArray3<f64>>,
    /// Index of the terminal event that stopped the run, or None.
    #[pyo3(get)]
    terminated_by: Option<usize>,
}

impl PyTrajectory {
    fn from_rust(py: Python<'_>, tr: Trajectory) -> PyResult<Self> {
        let shape = [tr.n_frames(), tr.n_particles, 3];
        let energy = tr.total_energy();
        let event_shape = [tr.events.len(), tr.n_particles, 3];
        let event_pos: Vec<Vec3> = tr
            .events
            .iter()
            .flat_map(|h| h.pos.iter().copied())
            .collect();
        let event_vel: Vec<Vec3> = tr
            .events
            .iter()
            .flat_map(|h| h.vel.iter().copied())
            .collect();
        Ok(Self {
            event_t: PyArray1::from_vec(py, tr.events.iter().map(|h| h.t).collect()).unbind(),
            event_index: PyArray1::from_vec(py, tr.events.iter().map(|h| h.event as i64).collect())
                .unbind(),
            event_pos: PyArray1::from_vec(py, flatten(&event_pos))
                .reshape(event_shape)?
                .unbind(),
            event_vel: PyArray1::from_vec(py, flatten(&event_vel))
                .reshape(event_shape)?
                .unbind(),
            terminated_by: tr.terminated_by,
            pos: PyArray1::from_vec(py, flatten(&tr.pos))
                .reshape(shape)?
                .unbind(),
            vel: PyArray1::from_vec(py, flatten(&tr.vel))
                .reshape(shape)?
                .unbind(),
            t: PyArray1::from_vec(py, tr.t).unbind(),
            kinetic: PyArray1::from_vec(py, tr.kinetic).unbind(),
            potential: PyArray1::from_vec(py, tr.potential).unbind(),
            energy: PyArray1::from_vec(py, energy).unbind(),
            n_particles: tr.n_particles,
        })
    }
}

#[pymethods]
impl PyTrajectory {
    #[getter]
    fn n_frames(&self, py: Python<'_>) -> usize {
        self.t.bind(py).len()
    }
    fn __len__(&self, py: Python<'_>) -> usize {
        self.n_frames(py)
    }
    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "Trajectory(n_frames={}, n_particles={})",
            self.n_frames(py),
            self.n_particles
        )
    }
}

// ---------------------------------------------------------------------------
// World

/// A simulation: particles, forces and a time integrator.
#[pyclass(name = "World", module = "physim")]
struct PyWorld {
    inner: World,
}

#[pymethods]
impl PyWorld {
    #[new]
    #[pyo3(signature = (integrator = "verlet"))]
    fn new(integrator: &str) -> PyResult<Self> {
        Ok(Self {
            inner: World::with_integrator(integrator)?,
        })
    }

    /// Adds a particle and returns its index.
    #[pyo3(signature = (pos, vel = None, mass = 1.0))]
    fn add_particle(
        &mut self,
        pos: &Bound<'_, PyAny>,
        vel: Option<&Bound<'_, PyAny>>,
        mass: f64,
    ) -> PyResult<usize> {
        let pos = extract_vec3(pos, "pos")?;
        let vel = vel
            .map(|v| extract_vec3(v, "vel"))
            .transpose()?
            .unwrap_or(Vec3::ZERO);
        Ok(self.inner.add_particle(pos, vel, mass)?)
    }

    /// Removes particle ``i``; higher indices shift down by one and index-based forces
    /// (springs) are renumbered. Fails while a force still refers to particle ``i``.
    fn remove_particle(&mut self, i: usize) -> PyResult<()> {
        Ok(self.inner.remove_particle(i)?)
    }

    /// Pins particle ``i`` in place (it still exerts forces), or releases it with ``pinned=False``.
    #[pyo3(signature = (i, pinned = true))]
    fn pin(&mut self, i: usize, pinned: bool) -> PyResult<()> {
        Ok(self.inner.pin(i, pinned)?)
    }

    /// Boolean mask of pinned particles, shape (N,).
    #[getter]
    fn pinned<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<bool>> {
        PyArray1::from_slice(py, &self.inner.state.pinned)
    }

    /// Adds a force and returns its id, for use with ``remove_force``, ``replace_force``,
    /// ``force_params`` and ``set_force_params``.
    fn add_force(&mut self, force: &Bound<'_, PyAny>) -> PyResult<ForceId> {
        Ok(self.inner.forces.add(build_force(force)?))
    }

    fn remove_force(&mut self, id: ForceId) -> PyResult<()> {
        self.inner.remove_force(id)?;
        Ok(())
    }

    /// Replaces force ``id`` with ``force``, keeping its id.
    fn replace_force(&mut self, id: ForceId, force: &Bound<'_, PyAny>) -> PyResult<()> {
        self.inner.forces.replace(id, build_force(force)?)?;
        Ok(())
    }

    /// Current tunable parameters of force ``id`` as a dict.
    fn force_params<'py>(&self, py: Python<'py>, id: ForceId) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (name, value) in self.inner.forces.get(id)?.params() {
            match value {
                Param::Scalar(x) => dict.set_item(name, x)?,
                Param::Vector(v) => dict.set_item(name, vec3_to_array(py, v))?,
            }
        }
        Ok(dict)
    }

    /// Changes parameters of force ``id``, e.g. ``set_force_params(g_id, G=2.0)``.
    /// Either all changes apply or none do.
    #[pyo3(signature = (id, **params))]
    fn set_force_params(
        &mut self,
        id: ForceId,
        params: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<()> {
        let mut values = Vec::new();
        for (name, value) in params.into_iter().flat_map(|d| d.iter()) {
            let name: String = name.extract()?;
            let value = match value.extract::<f64>() {
                Ok(x) => Param::Scalar(x),
                Err(_) => Param::Vector(extract_vec3(&value, &name)?),
            };
            values.push((name, value));
        }
        Ok(self.inner.set_force_params(id, &values)?)
    }

    fn clear_forces(&mut self) {
        self.inner.forces.clear();
    }

    /// Forces acting on the system as ``{id: name}``.
    #[getter]
    fn forces(&self) -> std::collections::BTreeMap<ForceId, String> {
        self.inner
            .forces
            .iter()
            .map(|(id, f)| (id, f.name()))
            .collect()
    }

    #[getter]
    fn integrator(&self) -> &'static str {
        self.inner.integrator().name()
    }

    #[setter]
    fn set_integrator(&mut self, name: &str) -> PyResult<()> {
        self.inner.set_integrator(integrators::by_name(name)?);
        Ok(())
    }

    /// ``{"name": ..., "order": ..., "symplectic": ...}`` for the current integrator.
    #[getter]
    fn integrator_info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let i = self.inner.integrator();
        let dict = PyDict::new(py);
        dict.set_item("name", i.name())?;
        dict.set_item("order", i.order())?;
        dict.set_item("symplectic", i.symplectic())?;
        Ok(dict)
    }

    /// Takes ``n`` steps of size ``dt``. If a step fails, the world is left at the last
    /// completed step.
    #[pyo3(signature = (dt, n = 1))]
    fn step(&mut self, py: Python<'_>, dt: f64, n: usize) -> PyResult<()> {
        let world = &mut self.inner;
        py.detach(|| (0..n).try_for_each(|_| world.step(dt)))?;
        Ok(())
    }

    /// Takes ``steps`` steps of size ``dt`` and returns a Trajectory containing the
    /// initial state, every ``record_every``-th step, and the final state.
    ///
    /// If a step fails, the world is left at the last completed step and the exception
    /// carries the frames recorded so far as its ``trajectory`` attribute.
    ///
    /// ``events`` is a list of :class:`Event`; detected events are returned in the
    /// trajectory's ``event_*`` arrays, and a terminal event stops the run at that moment.
    #[pyo3(signature = (dt, steps, record_every = 1, events = None))]
    fn run(
        &mut self,
        py: Python<'_>,
        dt: f64,
        steps: usize,
        record_every: usize,
        events: Option<Vec<PyRef<'_, PyEvent>>>,
    ) -> PyResult<PyTrajectory> {
        let events: Vec<Event> = events
            .unwrap_or_default()
            .iter()
            .map(|e| e.build(py))
            .collect();
        let world = &mut self.inner;
        match py.detach(|| world.run_with_events(dt, steps, record_every, &events)) {
            Ok(traj) => PyTrajectory::from_rust(py, traj),
            Err(RunFailure { error, trajectory }) => {
                let err = PyErr::from(error);
                if let Ok(partial) =
                    PyTrajectory::from_rust(py, *trajectory).and_then(|t| Py::new(py, t))
                {
                    // Some exception types reject new attributes; the error still propagates.
                    let _ = err.value(py).setattr("trajectory", partial);
                }
                Err(err)
            }
        }
    }

    #[getter]
    fn t(&self) -> f64 {
        self.inner.state.t
    }
    #[setter]
    fn set_t(&mut self, t: f64) {
        self.inner.state.t = t;
    }

    #[getter]
    fn n_particles(&self) -> usize {
        self.inner.state.len()
    }

    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        vecs_to_array(py, &self.inner.state.pos)
    }
    #[setter]
    fn set_positions(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        self.inner.state.pos = extract_vecs(value, self.inner.state.len(), "positions")?;
        Ok(())
    }

    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        vecs_to_array(py, &self.inner.state.vel)
    }
    #[setter]
    fn set_velocities(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let vel = extract_vecs(value, self.inner.state.len(), "velocities")?;
        Ok(self.inner.set_velocities(vel)?)
    }

    #[getter]
    fn masses<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.inner.state.mass)
    }
    #[setter]
    fn set_masses(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let arr = as_f64_array(value)?;
        let arr: PyReadonlyArray1<'_, f64> = arr
            .extract()
            .map_err(|_| PyValueError::new_err("masses: expected a 1D array"))?;
        Ok(self.inner.set_masses(arr.as_array().to_vec())?)
    }

    /// Current total acceleration of every particle, shape (N, 3).
    fn accelerations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        vecs_to_array(py, &self.inner.accelerations()?)
    }

    fn kinetic_energy(&self) -> f64 {
        self.inner.kinetic_energy()
    }
    fn potential_energy(&self) -> PyResult<f64> {
        Ok(self.inner.potential_energy()?)
    }
    fn total_energy(&self) -> PyResult<f64> {
        Ok(self.inner.total_energy()?)
    }
    fn momentum<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.inner.state.momentum())
    }
    /// Total angular momentum about the origin.
    fn angular_momentum<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.inner.state.angular_momentum())
    }
    fn center_of_mass<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        vec3_to_array(py, self.inner.state.center_of_mass())
    }

    fn __repr__(&self) -> String {
        format!(
            "World(n_particles={}, forces={:?}, integrator={:?}, t={})",
            self.inner.state.len(),
            self.forces().into_values().collect::<Vec<_>>(),
            self.integrator(),
            self.inner.state.t
        )
    }
}

// ---------------------------------------------------------------------------
// Events

enum EventSpec {
    Custom { function: Py<PyAny>, name: String },
    RadialVelocity(events::RadialVelocity),
    Coordinate(events::CoordinateCrossing),
    Separation(events::Separation),
}

/// An event for :meth:`World.run`: a moment when a scalar function ``g`` of the state
/// crosses zero. ``direction`` is +1 (``g`` rising through zero), -1 (falling) or 0
/// (either); ``terminal=True`` stops the run there.
///
/// ``Event(g)`` takes a Python function ``g(t, pos, vel, mass) -> float``. The static
/// constructors build Rust-side events that cost nothing per step from Python.
#[pyclass(frozen, name = "Event", module = "physim")]
struct PyEvent {
    spec: EventSpec,
    direction: Direction,
    terminal: bool,
}

fn direction(sign: i32) -> PyResult<Direction> {
    Ok(Direction::from_sign(sign)?)
}

#[pymethods]
impl PyEvent {
    #[new]
    #[pyo3(signature = (function, direction = 0, terminal = false, name = None))]
    fn new(
        function: &Bound<'_, PyAny>,
        direction: i32,
        terminal: bool,
        name: Option<String>,
    ) -> PyResult<Self> {
        if !function.is_callable() {
            return Err(PyValueError::new_err("Event function must be callable"));
        }
        let name = match name {
            Some(n) => n,
            None => function
                .getattr("__name__")
                .and_then(|n| n.extract())
                .unwrap_or_else(|_| "Event".to_string()),
        };
        Ok(Self {
            spec: EventSpec::Custom {
                function: function.clone().unbind(),
                name,
            },
            direction: self::direction(direction)?,
            terminal,
        })
    }

    /// ``(r_i - r_j) · (v_i - v_j)`` (``j=None``: relative to the origin). Rising
    /// crossings (``direction=+1``) are periapses, falling ones apoapses.
    #[staticmethod]
    #[pyo3(signature = (i, j = None, direction = 0, terminal = false))]
    fn radial_velocity(
        i: usize,
        j: Option<usize>,
        direction: i32,
        terminal: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            spec: EventSpec::RadialVelocity(events::RadialVelocity { i, j }),
            direction: self::direction(direction)?,
            terminal,
        })
    }

    /// ``pos[i, axis] - value``; ``axis`` is 0/1/2 or "x"/"y"/"z".
    #[staticmethod]
    #[pyo3(signature = (i, axis, value = 0.0, direction = 0, terminal = false))]
    fn coordinate(
        i: usize,
        axis: &Bound<'_, PyAny>,
        value: f64,
        direction: i32,
        terminal: bool,
    ) -> PyResult<Self> {
        let axis = match axis.extract::<usize>() {
            Ok(a) if a < 3 => a,
            _ => match axis.extract::<String>().as_deref() {
                Ok("x") => 0,
                Ok("y") => 1,
                Ok("z") => 2,
                _ => {
                    return Err(PyValueError::new_err(
                        "axis must be 0, 1, 2, 'x', 'y' or 'z'",
                    ))
                }
            },
        };
        Ok(Self {
            spec: EventSpec::Coordinate(events::CoordinateCrossing { i, axis, value }),
            direction: self::direction(direction)?,
            terminal,
        })
    }

    /// ``|r_i - r_j| - distance``: falling for approach to contact, rising for escape.
    #[staticmethod]
    #[pyo3(signature = (i, j, distance, direction = 0, terminal = false))]
    fn separation(
        i: usize,
        j: usize,
        distance: f64,
        direction: i32,
        terminal: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            spec: EventSpec::Separation(events::Separation { i, j, distance }),
            direction: self::direction(direction)?,
            terminal,
        })
    }

    #[getter]
    fn name(&self) -> String {
        match &self.spec {
            EventSpec::Custom { name, .. } => name.clone(),
            EventSpec::RadialVelocity(e) => e.name(),
            EventSpec::Coordinate(e) => e.name(),
            EventSpec::Separation(e) => e.name(),
        }
    }

    #[getter(direction)]
    fn direction_sign(&self) -> i32 {
        match self.direction {
            Direction::Rising => 1,
            Direction::Falling => -1,
            Direction::Either => 0,
        }
    }

    #[getter]
    fn terminal(&self) -> bool {
        self.terminal
    }

    fn __repr__(&self) -> String {
        format!(
            "Event({}, direction={}, terminal={})",
            self.name(),
            self.direction_sign(),
            if self.terminal { "True" } else { "False" }
        )
    }
}

impl PyEvent {
    fn build(&self, py: Python<'_>) -> Event {
        let function: Box<dyn EventFunction> = match &self.spec {
            EventSpec::Custom { function, name } => Box::new(PythonEvent {
                function: function.clone_ref(py),
                name: name.clone(),
            }),
            EventSpec::RadialVelocity(e) => Box::new(e.clone()),
            EventSpec::Coordinate(e) => Box::new(e.clone()),
            EventSpec::Separation(e) => Box::new(e.clone()),
        };
        Event {
            function,
            direction: self.direction,
            terminal: self.terminal,
        }
    }
}

struct PythonEvent {
    function: Py<PyAny>,
    name: String,
}

impl EventFunction for PythonEvent {
    fn value(&self, t: f64, pos: &[Vec3], vel: &[Vec3], mass: &[f64]) -> SimResult<f64> {
        Python::attach(|py| -> PyResult<f64> {
            self.function
                .bind(py)
                .call1((
                    t,
                    vecs_to_array(py, pos)?,
                    vecs_to_array(py, vel)?,
                    PyArray1::from_slice(py, mass),
                ))?
                .extract()
        })
        .map_err(SimError::Python)
    }

    fn name(&self) -> String {
        self.name.clone()
    }
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyEvent>()?;
    m.add_class::<PyWorld>()?;
    m.add_class::<PyTrajectory>()?;
    m.add_class::<PyUniformField>()?;
    m.add_class::<PyNewtonianGravity>()?;
    m.add_class::<PySpring>()?;
    m.add_class::<PyAnchorSpring>()?;
    m.add_class::<PyLinearDrag>()?;
    m.add_class::<PyQuadraticDrag>()?;
    m.add_class::<PyCustomForce>()?;
    m.add("INTEGRATORS", integrators::NAMES.to_vec())?;
    Ok(())
}
