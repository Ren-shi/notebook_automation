//! Python bindings (`physim._core`). Built only with the `python` feature.

use numpy::ndarray::ArrayView2;
use numpy::{
    PyArray1, PyArray2, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyReadonlyArray3, PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::adaptive::{AdaptiveOptions, Output};
use crate::chaos::LyapunovOptions;
use crate::checkpoint::{Checkpoint, SavedForce};
use crate::collisions::{Collisions, Wall};
use crate::constraints::{Anchor, ConstraintId, Constraints, Rod};
use crate::error::{Result as SimResult, SimError};
use crate::events::{self, Direction, Event, EventFunction};
use crate::forces::{self, BuiltinForce, Force, ForceId, Param};
use crate::integrators::{self, Op, Scheme};
use crate::state::State;
use crate::vec3::Vec3;
use crate::world::{Frame, Recorder, RunOptions, Trajectory, World};

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

/// Runs each parameter of a newly built force through its own `set_param` checks.
fn validated<F: Force>(mut force: F) -> PyResult<F> {
    for (name, value) in force.params() {
        force.set_param(name, value)?;
    }
    Ok(force)
}

/// ``Type(key=value, ...)`` from the force's checkpoint description.
fn builtin_repr(py: Python<'_>, force: &dyn Force) -> PyResult<String> {
    let d = describe_force(py, 0, force)?;
    let kind: String = d
        .get_item("type")?
        .map_or(Ok(String::new()), |t| t.extract())?;
    let mut args = Vec::new();
    for (k, v) in d.iter() {
        let k: String = k.extract()?;
        if k != "id" && k != "type" {
            args.push(format!("{k}={}", v.repr()?));
        }
    }
    Ok(format!("{kind}({})", args.join(", ")))
}

/// Particle index or fixed point, as accepted by ``World.add_rod``.
fn extract_anchor(to: &Bound<'_, PyAny>) -> PyResult<Anchor> {
    Ok(match to.extract::<usize>() {
        Ok(j) => Anchor::Particle(j),
        Err(_) => Anchor::Point(extract_vec3(to, "to")?),
    })
}

/// A scalar broadcast to ``n`` values, or an array of ``n`` values.
fn scalar_or_array(obj: &Bound<'_, PyAny>, n: usize, what: &str) -> PyResult<Vec<f64>> {
    if let Ok(x) = obj.extract::<f64>() {
        return Ok(vec![x; n]);
    }
    let arr = as_f64_array(obj)?;
    let arr: PyReadonlyArray1<'_, f64> = arr
        .extract()
        .map_err(|_| PyValueError::new_err(format!("{what}: expected a scalar or a 1D array")))?;
    Ok(arr.as_array().to_vec())
}

macro_rules! simple_spec {
    ($($ty:ty),*) => {
        $(impl PyForceSpec for $ty {
            fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
                Box::new(self.0.clone())
            }
        })*
    };
}

/// Spring with a dashpot between particles ``i`` and ``j``: force on ``i`` along the unit
/// bond vector ``n`` is ``[k (r - rest_length) + c (v_j - v_i)·n] n``.
#[pyclass(frozen, name = "DampedSpring", module = "physim")]
struct PyDampedSpring(forces::DampedSpring);

#[pymethods]
impl PyDampedSpring {
    #[new]
    #[pyo3(signature = (i, j, k, rest_length, c))]
    fn new(i: usize, j: usize, k: f64, rest_length: f64, c: f64) -> PyResult<Self> {
        Ok(Self(validated(forces::DampedSpring {
            i,
            j,
            k,
            rest_length,
            c,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Spring from particle ``i`` to particle or point ``to`` whose stiffness is modulated:
/// ``k(t) = k (1 + depth cos(omega t + phase))`` (parametric driving, Mathieu equation).
#[pyclass(frozen, name = "ModulatedSpring", module = "physim")]
struct PyModulatedSpring(forces::ModulatedSpring);

#[pymethods]
impl PyModulatedSpring {
    #[new]
    #[pyo3(signature = (i, to, k, depth, omega, phase = 0.0, rest_length = 0.0))]
    fn new(
        i: usize,
        to: &Bound<'_, PyAny>,
        k: f64,
        depth: f64,
        omega: f64,
        phase: f64,
        rest_length: f64,
    ) -> PyResult<Self> {
        let to = extract_anchor(to)?;
        if to == Anchor::Particle(i) {
            return Err(PyValueError::new_err(
                "ModulatedSpring: i and to must differ",
            ));
        }
        Ok(Self(validated(forces::ModulatedSpring {
            i,
            to,
            k,
            depth,
            omega,
            phase,
            rest_length,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Many Hookean springs in one force: bond ``b`` joins ``i[b]`` and ``j[b]`` with stiffness
/// ``k[b]`` and rest length ``rest_length[b]`` (``k`` and ``rest_length`` may be scalars).
/// Equivalent to one ``Spring`` per bond, but much faster for lattices and polymers.
#[pyclass(frozen, name = "SpringNetwork", module = "physim")]
struct PySpringNetwork(forces::SpringNetwork);

#[pymethods]
impl PySpringNetwork {
    #[new]
    fn new(
        i: Vec<usize>,
        j: Vec<usize>,
        k: &Bound<'_, PyAny>,
        rest_length: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let n = i.len();
        let k = scalar_or_array(k, n, "k")?;
        let rest_length = scalar_or_array(rest_length, n, "rest_length")?;
        Ok(Self(forces::SpringNetwork::new(i, j, k, rest_length)?))
    }
    fn __len__(&self) -> usize {
        self.0.len()
    }
    fn __repr__(&self) -> String {
        format!("SpringNetwork({} bonds)", self.0.len())
    }
}

/// Power-law potential ``Φ = k r^n`` per unit mass about ``center`` (``k ln r`` for
/// ``n = 0``), acting on every particle.
#[pyclass(frozen, name = "PowerLaw", module = "physim")]
struct PyPowerLaw(forces::PowerLaw);

#[pymethods]
impl PyPowerLaw {
    #[new]
    #[pyo3(signature = (k, n, center = None))]
    fn new(k: f64, n: f64, center: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let center = center.map_or(Ok(Vec3::ZERO), |c| extract_vec3(c, "center"))?;
        Ok(Self(validated(forces::PowerLaw { center, k, n })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Yukawa (screened) potential ``Φ = -k exp(-r/length) / r`` per unit mass about ``center``.
#[pyclass(frozen, name = "Yukawa", module = "physim")]
struct PyYukawa(forces::Yukawa);

#[pymethods]
impl PyYukawa {
    #[new]
    #[pyo3(signature = (k, length, center = None))]
    fn new(k: f64, length: f64, center: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let center = center.map_or(Ok(Vec3::ZERO), |c| extract_vec3(c, "center"))?;
        Ok(Self(validated(forces::Yukawa { center, k, length })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Plummer sphere potential ``Φ = -GM / sqrt(r² + a²)`` per unit mass about ``center``.
#[pyclass(frozen, name = "PlummerPotential", module = "physim")]
struct PyPlummerPotential(forces::PlummerPotential);

#[pymethods]
impl PyPlummerPotential {
    #[new]
    #[pyo3(signature = (GM, a, center = None))]
    #[allow(non_snake_case)]
    fn new(GM: f64, a: f64, center: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let center = center.map_or(Ok(Vec3::ZERO), |c| extract_vec3(c, "center"))?;
        Ok(Self(validated(forces::PlummerPotential {
            center,
            gm: GM,
            a,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Hernquist potential ``Φ = -GM / (r + a)`` per unit mass about ``center``.
#[pyclass(frozen, name = "HernquistPotential", module = "physim")]
struct PyHernquistPotential(forces::HernquistPotential);

#[pymethods]
impl PyHernquistPotential {
    #[new]
    #[pyo3(signature = (GM, a, center = None))]
    #[allow(non_snake_case)]
    fn new(GM: f64, a: f64, center: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let center = center.map_or(Ok(Vec3::ZERO), |c| extract_vec3(c, "center"))?;
        Ok(Self(validated(forces::HernquistPotential {
            center,
            gm: GM,
            a,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Harmonic trap ``Φ = ½ Σ ω_k² (x_k - c_k)²`` per unit mass, acting on every particle.
/// ``omega`` is a scalar (isotropic) or one angular frequency per axis.
#[pyclass(frozen, name = "HarmonicTrap", module = "physim")]
struct PyHarmonicTrap(forces::HarmonicTrap);

#[pymethods]
impl PyHarmonicTrap {
    #[new]
    #[pyo3(signature = (omega, center = None))]
    fn new(omega: &Bound<'_, PyAny>, center: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let omega = match omega.extract::<f64>() {
            Ok(w) => Vec3::new(w, w, w),
            Err(_) => extract_vec3(omega, "omega")?,
        };
        let center = center.map_or(Ok(Vec3::ZERO), |c| extract_vec3(c, "center"))?;
        Ok(Self(validated(forces::HarmonicTrap { center, omega })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Sinusoidal force on particle ``i``: ``F(t) = amplitude cos(omega t + phase)``.
#[pyclass(frozen, name = "PeriodicForce", module = "physim")]
struct PyPeriodicForce(forces::PeriodicForce);

#[pymethods]
impl PyPeriodicForce {
    #[new]
    #[pyo3(signature = (i, amplitude, omega, phase = 0.0))]
    fn new(i: usize, amplitude: &Bound<'_, PyAny>, omega: f64, phase: f64) -> PyResult<Self> {
        let amplitude = extract_vec3(amplitude, "amplitude")?;
        Ok(Self(validated(forces::PeriodicForce {
            i,
            amplitude,
            omega,
            phase,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// First post-Newtonian (general-relativistic) correction from particle ``central``, in the
/// test-particle limit. Add it alongside ``NewtonianGravity``; ``c`` is the speed of light
/// in simulation units.
#[pyclass(frozen, name = "PostNewtonian", module = "physim")]
struct PyPostNewtonian(forces::PostNewtonian);

#[pymethods]
impl PyPostNewtonian {
    #[new]
    #[pyo3(signature = (central, c, G = 1.0))]
    #[allow(non_snake_case)]
    fn new(central: usize, c: f64, G: f64) -> PyResult<Self> {
        Ok(Self(validated(forces::PostNewtonian { central, g: G, c })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// J2 oblateness of particle ``central`` (equatorial radius ``radius``, spin ``axis``).
/// Add it alongside ``NewtonianGravity``.
#[pyclass(frozen, name = "J2Oblateness", module = "physim")]
struct PyJ2Oblateness(forces::J2Oblateness);

#[pymethods]
impl PyJ2Oblateness {
    #[new]
    #[pyo3(signature = (central, J2, radius, axis = None, G = 1.0))]
    #[allow(non_snake_case)]
    fn new(
        central: usize,
        J2: f64,
        radius: f64,
        axis: Option<&Bound<'_, PyAny>>,
        G: f64,
    ) -> PyResult<Self> {
        let axis = axis.map_or(Ok(Vec3::new(0.0, 0.0, 1.0)), |a| extract_vec3(a, "axis"))?;
        Ok(Self(validated(forces::J2Oblateness {
            central,
            g: G,
            j2: J2,
            radius,
            axis,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Hénon-Heiles potential ``Φ = ½(x² + y²) + lam (x² y - y³/3)`` per unit mass in the
/// xy-plane about ``center``, acting on every particle (``lam = 1`` is the classic system).
#[pyclass(frozen, name = "HenonHeiles", module = "physim")]
struct PyHenonHeiles(forces::HenonHeiles);

#[pymethods]
impl PyHenonHeiles {
    #[new]
    #[pyo3(signature = (lam = 1.0, center = None))]
    fn new(lam: f64, center: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        let center = center.map_or(Ok(Vec3::ZERO), |c| extract_vec3(c, "center"))?;
        Ok(Self(validated(forces::HenonHeiles {
            center,
            lambda: lam,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Uniform electric field: ``a = (q/m) E`` on every charged particle.
#[pyclass(frozen, name = "ElectricField", module = "physim")]
struct PyElectricField(forces::ElectricField);

#[pymethods]
impl PyElectricField {
    #[new]
    #[pyo3(signature = (E))]
    #[allow(non_snake_case)]
    fn new(E: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(forces::ElectricField {
            e: extract_vec3(E, "E")?,
        }))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Uniform magnetic field: ``a = (q/m) v × B`` on every charged particle. Use the ``boris``
/// integrator to conserve kinetic energy exactly.
#[pyclass(frozen, name = "MagneticField", module = "physim")]
struct PyMagneticField(forces::MagneticField);

#[pymethods]
impl PyMagneticField {
    #[new]
    #[pyo3(signature = (B))]
    #[allow(non_snake_case)]
    fn new(B: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(forces::MagneticField {
            b: extract_vec3(B, "B")?,
        }))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Pairwise Coulomb interaction ``U = k Σ q_i q_j / sqrt(r² + softening²)`` (like charges
/// repel); same O(N²) parallel pair loop as ``NewtonianGravity``.
#[pyclass(frozen, name = "Coulomb", module = "physim")]
struct PyCoulomb(forces::Coulomb);

#[pymethods]
impl PyCoulomb {
    #[new]
    #[pyo3(signature = (k = 1.0, softening = 0.0))]
    fn new(k: f64, softening: f64) -> PyResult<Self> {
        Ok(Self(validated(forces::Coulomb { k, softening })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Electric and/or magnetic fields given by Python functions ``E(t, pos)`` and
/// ``B(t, pos)``, each returning an ``(N, 3)`` array of the field at every particle position.
/// The Lorentz force ``a = (q/m)(E + v × B)`` follows, and ``boris`` rotates about ``B``
/// exactly. Like ``CustomForce`` it cannot be saved in checkpoints.
#[pyclass(frozen, name = "FieldForce", module = "physim")]
struct PyFieldForce {
    e: Option<Py<PyAny>>,
    b: Option<Py<PyAny>>,
    name: String,
}

#[pymethods]
impl PyFieldForce {
    #[new]
    #[pyo3(signature = (E = None, B = None, name = "FieldForce"))]
    #[allow(non_snake_case)]
    fn new(
        E: Option<&Bound<'_, PyAny>>,
        B: Option<&Bound<'_, PyAny>>,
        name: &str,
    ) -> PyResult<Self> {
        if E.is_none() && B.is_none() {
            return Err(PyValueError::new_err("FieldForce needs E, B or both"));
        }
        if E.is_some_and(|f| !f.is_callable()) || B.is_some_and(|f| !f.is_callable()) {
            return Err(PyValueError::new_err(
                "E and B must be callables E(t, pos) -> (N, 3)",
            ));
        }
        Ok(Self {
            e: E.map(|f| f.clone().unbind()),
            b: B.map(|f| f.clone().unbind()),
            name: name.to_string(),
        })
    }
    fn __repr__(&self) -> String {
        format!("FieldForce({})", self.name)
    }
}

/// Wraps a Python field function `f(t, pos) -> (N, 3)` as a Rust field closure.
fn python_field(
    f: Py<PyAny>,
    what: String,
) -> impl Fn(f64, &[Vec3], &mut [Vec3]) -> SimResult<()> + Send + Sync + 'static {
    move |t, pos, out| {
        Python::attach(|py| -> PyResult<()> {
            let value = f.bind(py).call1((t, vecs_to_array(py, pos)?))?;
            for (o, x) in out.iter_mut().zip(extract_vecs(&value, pos.len(), &what)?) {
                *o += x;
            }
            Ok(())
        })
        .map_err(SimError::Python)
    }
}

impl PyForceSpec for PyFieldForce {
    fn build(&self, py: Python<'_>) -> Box<dyn Force> {
        let mut f = forces::FieldFunctions::new(self.name.clone());
        if let Some(e) = &self.e {
            f = f.electric(python_field(e.clone_ref(py), format!("{} E", self.name)));
        }
        if let Some(b) = &self.b {
            f = f.magnetic(python_field(b.clone_ref(py), format!("{} B", self.name)));
        }
        Box::new(f)
    }
}

/// Barnes-Hut tree gravity: the same physics as ``NewtonianGravity`` at O(N log N) cost.
/// ``theta`` in (0, 1] sets the accuracy (rms force error ~1e-3 at 0.5, ~1e-2 at 1.0);
/// ``quadrupole=True`` adds each node's quadrupole (about 3x more accurate at the same
/// ``theta``, 2.5x slower). Momentum is conserved only to the force accuracy.
#[pyclass(frozen, name = "TreeGravity", module = "physim")]
struct PyTreeGravity(forces::TreeGravity);

#[pymethods]
impl PyTreeGravity {
    #[new]
    #[pyo3(signature = (G = 1.0, softening = 0.0, theta = 0.5, quadrupole = false))]
    #[allow(non_snake_case)]
    fn new(G: f64, softening: f64, theta: f64, quadrupole: bool) -> PyResult<Self> {
        Ok(Self(validated(forces::TreeGravity {
            g: G,
            softening,
            theta,
            quadrupole,
        })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Soft contact between particles with a radius: overlapping spheres (overlap ``δ``) repel with
/// ``k δ`` (``law="linear"``) or ``k δ^1.5`` (``law="hertz"``) plus a dashpot ``-damping * u_n``
/// on the normal approach speed. Set radii with ``add_particle(..., radius=r)`` or
/// ``w.radii``. Use steps well below the contact duration (``π sqrt(μ/k)`` for the linear law).
#[pyclass(frozen, name = "SoftContact", module = "physim")]
struct PySoftContact(forces::SoftContact);

#[pymethods]
impl PySoftContact {
    #[new]
    #[pyo3(signature = (k, damping = 0.0, law = "linear"))]
    fn new(k: f64, damping: f64, law: &str) -> PyResult<Self> {
        let law = match law {
            "linear" => forces::ContactLaw::Linear,
            "hertz" => forces::ContactLaw::Hertz,
            other => {
                return Err(PyValueError::new_err(format!(
                    "law must be \"linear\" or \"hertz\", got {other:?}"
                )))
            }
        };
        Ok(Self(validated(forces::SoftContact { k, damping, law })?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

fn period_from(box_: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Vec3>> {
    box_.map(|b| match b.extract::<f64>() {
        Ok(l) => Ok(Vec3::new(l, l, l)),
        Err(_) => extract_vec3(b, "box"),
    })
    .transpose()
}

fn checked_pair(p: forces::PairPotential) -> PyResult<forces::PairPotential> {
    // Validate by evaluating on an empty system.
    p.potential(0.0, &[], &[])?;
    Ok(p)
}

/// Lennard-Jones ``4 epsilon [(sigma/r)^12 - (sigma/r)^6]`` between every pair closer than
/// ``cutoff`` (shifted to zero there if ``shift``), optionally in the periodic box ``box``
/// (side lengths: a number or a 3-vector; minimum image, coordinates need not be wrapped).
/// O(N) via cell lists and Verlet neighbour lists.
#[pyclass(frozen, name = "LennardJones", module = "physim")]
struct PyLennardJones(forces::PairPotential);

#[pymethods]
impl PyLennardJones {
    #[new]
    #[pyo3(signature = (epsilon = 1.0, sigma = 1.0, cutoff = 2.5, shift = true, r#box = None))]
    fn new(
        epsilon: f64,
        sigma: f64,
        cutoff: f64,
        shift: bool,
        r#box: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self(checked_pair(forces::PairPotential::new(
            forces::PairKind::LennardJones { epsilon, sigma },
            cutoff,
            shift,
            period_from(r#box)?,
        ))?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// Morse ``depth [(1 - exp(-a (r - r0)))^2 - 1]`` between every pair closer than ``cutoff``;
/// ``shift`` and ``box`` as for :class:`LennardJones`.
#[pyclass(frozen, name = "Morse", module = "physim")]
struct PyMorse(forces::PairPotential);

#[pymethods]
impl PyMorse {
    #[new]
    #[pyo3(signature = (depth, a, r0, cutoff, shift = true, r#box = None))]
    fn new(
        depth: f64,
        a: f64,
        r0: f64,
        cutoff: f64,
        shift: bool,
        r#box: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self(checked_pair(forces::PairPotential::new(
            forces::PairKind::Morse { depth, a, r0 },
            cutoff,
            shift,
            period_from(r#box)?,
        ))?))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

/// A pair potential given as a table: ``V`` at equally spaced ``r`` (``r_min``, step ``dr``),
/// with ``dV`` its derivative at the same points (see :func:`physim.tabulate_pair` to build
/// one from any vectorised Python function). Cubic Hermite interpolation, so forces are
/// continuous; evaluated in Rust with no Python calls. ``cutoff``, ``shift``, ``box`` as for
/// :class:`LennardJones`.
#[pyclass(frozen, name = "TabulatedPair", module = "physim")]
struct PyTabulatedPair(forces::PairPotential);

#[pymethods]
impl PyTabulatedPair {
    #[new]
    #[pyo3(signature = (r_min, dr, V, dV, cutoff, shift = true, r#box = None))]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn new(
        r_min: f64,
        dr: f64,
        V: Vec<f64>,
        dV: Vec<f64>,
        cutoff: f64,
        shift: bool,
        r#box: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self(checked_pair(forces::PairPotential::new(
            forces::PairKind::Table {
                r_min,
                dr,
                v: V,
                dv: dV,
            },
            cutoff,
            shift,
            period_from(r#box)?,
        ))?))
    }
    fn __repr__(&self) -> String {
        format!("TabulatedPair(cutoff={})", self.0.cutoff)
    }
}

simple_spec!(
    PyDampedSpring,
    PyModulatedSpring,
    PySpringNetwork,
    PyPowerLaw,
    PyYukawa,
    PyPlummerPotential,
    PyHernquistPotential,
    PyHarmonicTrap,
    PyPeriodicForce,
    PyPostNewtonian,
    PyJ2Oblateness,
    PyHenonHeiles,
    PyElectricField,
    PyMagneticField,
    PyCoulomb,
    PyTreeGravity,
    PySoftContact,
    PyLennardJones,
    PyMorse,
    PyTabulatedPair
);

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
        PyDampedSpring,
        PyModulatedSpring,
        PySpringNetwork,
        PyPowerLaw,
        PyYukawa,
        PyPlummerPotential,
        PyHernquistPotential,
        PyHarmonicTrap,
        PyPeriodicForce,
        PyPostNewtonian,
        PyJ2Oblateness,
        PyHenonHeiles,
        PyElectricField,
        PyMagneticField,
        PyCoulomb,
        PyTreeGravity,
        PySoftContact,
        PyLennardJones,
        PyMorse,
        PyTabulatedPair,
        PyFieldForce,
        PyCustomForce
    );
    Err(PyValueError::new_err(format!(
        "not a physim force: {}; wrap Python functions in physim.CustomForce",
        force.repr()?
    )))
}

/// ``{"id": id, "type": "Spring", "i": 0, ...}``: the constructor name and keyword
/// arguments of a built-in force, or ``{"id": id, "type": "external", "name": ...}``.
fn describe_force<'py>(
    py: Python<'py>,
    id: ForceId,
    force: &dyn Force,
) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("id", id)?;
    let v = |v: Vec3| v.to_array().to_vec();
    match force.builtin() {
        None => {
            d.set_item("type", "external")?;
            d.set_item("name", force.name())?;
        }
        Some(BuiltinForce::UniformField(f)) => {
            d.set_item("type", "UniformField")?;
            d.set_item("g", v(f.g))?;
        }
        Some(BuiltinForce::NewtonianGravity(f)) => {
            d.set_item("type", "NewtonianGravity")?;
            d.set_item("G", f.g)?;
            d.set_item("softening", f.softening)?;
        }
        Some(BuiltinForce::Spring(f)) => {
            d.set_item("type", "Spring")?;
            d.set_item("i", f.i)?;
            d.set_item("j", f.j)?;
            d.set_item("k", f.k)?;
            d.set_item("rest_length", f.rest_length)?;
        }
        Some(BuiltinForce::AnchorSpring(f)) => {
            d.set_item("type", "AnchorSpring")?;
            d.set_item("i", f.i)?;
            d.set_item("anchor", v(f.anchor))?;
            d.set_item("k", f.k)?;
            d.set_item("rest_length", f.rest_length)?;
        }
        Some(BuiltinForce::LinearDrag(f)) => {
            d.set_item("type", "LinearDrag")?;
            d.set_item("gamma", f.gamma)?;
        }
        Some(BuiltinForce::QuadraticDrag(f)) => {
            d.set_item("type", "QuadraticDrag")?;
            d.set_item("c", f.c)?;
        }
        Some(BuiltinForce::DampedSpring(f)) => {
            d.set_item("type", "DampedSpring")?;
            d.set_item("i", f.i)?;
            d.set_item("j", f.j)?;
            d.set_item("k", f.k)?;
            d.set_item("rest_length", f.rest_length)?;
            d.set_item("c", f.c)?;
        }
        Some(BuiltinForce::ModulatedSpring(f)) => {
            d.set_item("type", "ModulatedSpring")?;
            d.set_item("i", f.i)?;
            match f.to {
                Anchor::Particle(j) => d.set_item("to", j)?,
                Anchor::Point(p) => d.set_item("to", v(p))?,
            }
            d.set_item("k", f.k)?;
            d.set_item("depth", f.depth)?;
            d.set_item("omega", f.omega)?;
            d.set_item("phase", f.phase)?;
            d.set_item("rest_length", f.rest_length)?;
        }
        Some(BuiltinForce::SpringNetwork(f)) => {
            d.set_item("type", "SpringNetwork")?;
            d.set_item("i", f.i().to_vec())?;
            d.set_item("j", f.j().to_vec())?;
            d.set_item("k", f.k().to_vec())?;
            d.set_item("rest_length", f.rest_length().to_vec())?;
        }
        Some(BuiltinForce::PowerLaw(f)) => {
            d.set_item("type", "PowerLaw")?;
            d.set_item("k", f.k)?;
            d.set_item("n", f.n)?;
            d.set_item("center", v(f.center))?;
        }
        Some(BuiltinForce::Yukawa(f)) => {
            d.set_item("type", "Yukawa")?;
            d.set_item("k", f.k)?;
            d.set_item("length", f.length)?;
            d.set_item("center", v(f.center))?;
        }
        Some(BuiltinForce::PlummerPotential(f)) => {
            d.set_item("type", "PlummerPotential")?;
            d.set_item("GM", f.gm)?;
            d.set_item("a", f.a)?;
            d.set_item("center", v(f.center))?;
        }
        Some(BuiltinForce::HernquistPotential(f)) => {
            d.set_item("type", "HernquistPotential")?;
            d.set_item("GM", f.gm)?;
            d.set_item("a", f.a)?;
            d.set_item("center", v(f.center))?;
        }
        Some(BuiltinForce::HarmonicTrap(f)) => {
            d.set_item("type", "HarmonicTrap")?;
            d.set_item("omega", v(f.omega))?;
            d.set_item("center", v(f.center))?;
        }
        Some(BuiltinForce::PeriodicForce(f)) => {
            d.set_item("type", "PeriodicForce")?;
            d.set_item("i", f.i)?;
            d.set_item("amplitude", v(f.amplitude))?;
            d.set_item("omega", f.omega)?;
            d.set_item("phase", f.phase)?;
        }
        Some(BuiltinForce::PostNewtonian(f)) => {
            d.set_item("type", "PostNewtonian")?;
            d.set_item("central", f.central)?;
            d.set_item("c", f.c)?;
            d.set_item("G", f.g)?;
        }
        Some(BuiltinForce::ElectricField(f)) => {
            d.set_item("type", "ElectricField")?;
            d.set_item("E", v(f.e))?;
        }
        Some(BuiltinForce::MagneticField(f)) => {
            d.set_item("type", "MagneticField")?;
            d.set_item("B", v(f.b))?;
        }
        Some(BuiltinForce::Coulomb(f)) => {
            d.set_item("type", "Coulomb")?;
            d.set_item("k", f.k)?;
            d.set_item("softening", f.softening)?;
        }
        Some(BuiltinForce::TreeGravity(f)) => {
            d.set_item("type", "TreeGravity")?;
            d.set_item("G", f.g)?;
            d.set_item("softening", f.softening)?;
            d.set_item("theta", f.theta)?;
            d.set_item("quadrupole", f.quadrupole)?;
        }
        Some(BuiltinForce::SoftContact(f)) => {
            d.set_item("type", "SoftContact")?;
            d.set_item("k", f.k)?;
            d.set_item("damping", f.damping)?;
            d.set_item(
                "law",
                match f.law {
                    forces::ContactLaw::Linear => "linear",
                    forces::ContactLaw::Hertz => "hertz",
                },
            )?;
        }
        Some(BuiltinForce::PairPotential(f)) => {
            match &f.kind {
                forces::PairKind::LennardJones { epsilon, sigma } => {
                    d.set_item("type", "LennardJones")?;
                    d.set_item("epsilon", epsilon)?;
                    d.set_item("sigma", sigma)?;
                }
                forces::PairKind::Morse { depth, a, r0 } => {
                    d.set_item("type", "Morse")?;
                    d.set_item("depth", depth)?;
                    d.set_item("a", a)?;
                    d.set_item("r0", r0)?;
                }
                forces::PairKind::Table {
                    r_min,
                    dr,
                    v: vt,
                    dv,
                } => {
                    d.set_item("type", "TabulatedPair")?;
                    d.set_item("r_min", r_min)?;
                    d.set_item("dr", dr)?;
                    d.set_item("V", vt.clone())?;
                    d.set_item("dV", dv.clone())?;
                }
            }
            d.set_item("cutoff", f.cutoff)?;
            d.set_item("shift", f.shift)?;
            match f.period {
                Some(l) => d.set_item("box", v(l))?,
                None => d.set_item("box", py.None())?,
            }
        }
        Some(BuiltinForce::HenonHeiles(f)) => {
            d.set_item("type", "HenonHeiles")?;
            d.set_item("lam", f.lambda)?;
            d.set_item("center", v(f.center))?;
        }
        Some(BuiltinForce::J2Oblateness(f)) => {
            d.set_item("type", "J2Oblateness")?;
            d.set_item("central", f.central)?;
            d.set_item("J2", f.j2)?;
            d.set_item("radius", f.radius)?;
            d.set_item("axis", v(f.axis))?;
            d.set_item("G", f.g)?;
        }
    }
    Ok(d)
}

fn describe_forces<'py>(py: Python<'py>, world: &World) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for (id, f) in world.forces.iter() {
        list.append(describe_force(py, id, f)?)?;
    }
    Ok(list)
}

/// Inverse of [`describe_force`] for built-in forces: calls the constructor named by "type".
fn builtin_from_description(desc: &Bound<'_, PyDict>) -> PyResult<BuiltinForce> {
    let py = desc.py();
    let kind: String = desc
        .get_item("type")?
        .ok_or_else(|| PyValueError::new_err("force description has no \"type\""))?
        .extract()?;
    const BUILTIN: [&str; 26] = [
        "LennardJones",
        "Morse",
        "TabulatedPair",
        "SoftContact",
        "TreeGravity",
        "HenonHeiles",
        "ElectricField",
        "MagneticField",
        "Coulomb",
        "UniformField",
        "NewtonianGravity",
        "Spring",
        "AnchorSpring",
        "LinearDrag",
        "QuadraticDrag",
        "DampedSpring",
        "ModulatedSpring",
        "SpringNetwork",
        "PowerLaw",
        "Yukawa",
        "PlummerPotential",
        "HernquistPotential",
        "HarmonicTrap",
        "PeriodicForce",
        "PostNewtonian",
        "J2Oblateness",
    ];
    if !BUILTIN.contains(&kind.as_str()) {
        return Err(PyValueError::new_err(format!(
            "unknown force type {kind:?}"
        )));
    }
    let kwargs = desc.copy()?;
    kwargs.del_item("type")?;
    if kwargs.contains("id")? {
        kwargs.del_item("id")?;
    }
    let obj = py
        .import("physim._core")?
        .getattr(kind.as_str())?
        .call((), Some(&kwargs))?;
    build_force(&obj)?
        .builtin()
        .ok_or_else(|| PyValueError::new_err(format!("{kind} is not a built-in force")))
}

/// ``[{"id": id, "i": i, "j": j, "length": L}, ...]``; a rod to a fixed point has
/// ``"anchor": [x, y, z]`` in place of ``"j"``.
fn describe_constraints<'py>(py: Python<'py>, c: &Constraints) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for (id, r) in c.iter() {
        let d = PyDict::new(py);
        d.set_item("id", id)?;
        d.set_item("i", r.i)?;
        match r.anchor {
            Anchor::Particle(j) => d.set_item("j", j)?,
            Anchor::Point(p) => d.set_item("anchor", p.to_array().to_vec())?,
        }
        d.set_item("length", r.length)?;
        list.append(d)?;
    }
    Ok(list)
}

fn rod_from_description(desc: &Bound<'_, PyDict>) -> PyResult<(ConstraintId, Rod)> {
    let field = |k: &str| -> PyResult<Bound<'_, PyAny>> {
        desc.get_item(k)?
            .ok_or_else(|| PyValueError::new_err(format!("constraint description has no {k:?}")))
    };
    let anchor = match desc.get_item("j")? {
        Some(j) => Anchor::Particle(j.extract()?),
        None => Anchor::Point(extract_vec3(&field("anchor")?, "anchor")?),
    };
    let rod = Rod {
        i: field("i")?.extract()?,
        anchor,
        length: field("length")?.extract()?,
    };
    Ok((field("id")?.extract()?, rod))
}

/// `d[key]` if present and not None.
fn optional<'py, T: for<'a> FromPyObject<'a, 'py>>(
    d: &Bound<'py, PyDict>,
    key: &str,
) -> PyResult<Option<T>> {
    match d.get_item(key)? {
        Some(v) if !v.is_none() => Ok(Some(v.extract().map_err(Into::into)?)),
        _ => Ok(None),
    }
}

// ---------------------------------------------------------------------------
// Trajectory

/// Recorded frames of a run: ``t`` (F,), ``pos``/``vel`` (F, N, 3),
/// ``kinetic``/``potential``/``energy`` (F,) or None if the run skipped energies.
///
/// ``metadata`` describes the run (engine version, integrator, dt, forces...). It is
/// saved with the trajectory by :func:`physim.save_trajectory`.
#[pyclass(frozen, name = "Trajectory", module = "physim")]
struct PyTrajectory {
    #[pyo3(get)]
    t: Py<PyArray1<f64>>,
    #[pyo3(get)]
    pos: Py<PyArray3<f64>>,
    #[pyo3(get)]
    vel: Py<PyArray3<f64>>,
    #[pyo3(get)]
    kinetic: Option<Py<PyArray1<f64>>>,
    #[pyo3(get)]
    potential: Option<Py<PyArray1<f64>>>,
    #[pyo3(get)]
    energy: Option<Py<PyArray1<f64>>>,
    /// Tension in each constraint, (F, C): the force pulling each rod's ends together.
    #[pyo3(get)]
    tension: Py<PyArray2<f64>>,
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
    #[pyo3(get)]
    metadata: Py<PyDict>,
}

impl PyTrajectory {
    fn from_rust(py: Python<'_>, tr: Trajectory, metadata: Py<PyDict>) -> PyResult<Self> {
        let shape = [tr.n_frames(), tr.n_particles, 3];
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
        let (kinetic, potential, energy) = if tr.energies {
            let energy = tr.total_energy();
            (
                Some(PyArray1::from_vec(py, tr.kinetic).unbind()),
                Some(PyArray1::from_vec(py, tr.potential).unbind()),
                Some(PyArray1::from_vec(py, energy).unbind()),
            )
        } else {
            (None, None, None)
        };
        let tension = PyArray1::from_vec(py, tr.tension)
            .reshape([tr.t.len(), tr.n_constraints])?
            .unbind();
        Ok(Self {
            tension,
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
            kinetic,
            potential,
            energy,
            n_particles: tr.n_particles,
            metadata,
        })
    }
}

fn extract_f64s(obj: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<f64>> {
    let arr = as_f64_array(obj)?;
    let arr: PyReadonlyArray1<'_, f64> = arr
        .extract()
        .map_err(|_| PyValueError::new_err(format!("{what}: expected a 1D array")))?;
    Ok(arr.as_array().to_vec())
}

/// A ``(frames, n_particles, 3)`` array as frame-major vectors.
fn extract_frames(
    obj: &Bound<'_, PyAny>,
    frames: Option<usize>,
    n_particles: Option<usize>,
    what: &str,
) -> PyResult<(usize, usize, Vec<Vec3>)> {
    let arr = as_f64_array(obj)?;
    let arr: PyReadonlyArray3<'_, f64> = arr.extract().map_err(|_| {
        PyValueError::new_err(format!(
            "{what}: expected a 3D array (frames, particles, 3)"
        ))
    })?;
    let shape = arr.shape().to_vec();
    if shape[2] != 3
        || frames.is_some_and(|f| f != shape[0])
        || n_particles.is_some_and(|n| n != shape[1])
    {
        return Err(PyValueError::new_err(format!(
            "{what}: expected shape ({}, {}, 3), got {shape:?}",
            frames.map_or("F".into(), |f| f.to_string()),
            n_particles.map_or("N".into(), |n| n.to_string()),
        )));
    }
    let vecs = arr
        .as_array()
        .rows()
        .into_iter()
        .map(|r| Vec3::new(r[0], r[1], r[2]))
        .collect();
    Ok((shape[0], shape[1], vecs))
}

#[pymethods]
impl PyTrajectory {
    /// Builds a trajectory from arrays, e.g. ones loaded from a file.
    #[new]
    #[pyo3(signature = (
        t, pos, vel, kinetic = None, potential = None, *, tension = None, event_t = None,
        event_index = None, event_pos = None, event_vel = None, terminated_by = None, metadata = None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn py_new(
        py: Python<'_>,
        t: &Bound<'_, PyAny>,
        pos: &Bound<'_, PyAny>,
        vel: &Bound<'_, PyAny>,
        kinetic: Option<&Bound<'_, PyAny>>,
        potential: Option<&Bound<'_, PyAny>>,
        tension: Option<&Bound<'_, PyAny>>,
        event_t: Option<&Bound<'_, PyAny>>,
        event_index: Option<Vec<usize>>,
        event_pos: Option<&Bound<'_, PyAny>>,
        event_vel: Option<&Bound<'_, PyAny>>,
        terminated_by: Option<usize>,
        metadata: Option<Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let t = extract_f64s(t, "t")?;
        let (_, n, pos) = extract_frames(pos, Some(t.len()), None, "pos")?;
        let (_, _, vel) = extract_frames(vel, Some(t.len()), Some(n), "vel")?;
        let mut tr = Trajectory::new(n, kinetic.is_some() || potential.is_some());
        tr.t = t;
        tr.pos = pos;
        tr.vel = vel;
        if tr.energies {
            let (Some(k), Some(u)) = (kinetic, potential) else {
                return Err(PyValueError::new_err(
                    "give both kinetic and potential, or neither",
                ));
            };
            tr.kinetic = extract_f64s(k, "kinetic")?;
            tr.potential = extract_f64s(u, "potential")?;
            if tr.kinetic.len() != tr.n_frames() || tr.potential.len() != tr.n_frames() {
                return Err(PyValueError::new_err(
                    "kinetic and potential need one value per frame",
                ));
            }
        }
        if let Some(tension) = tension {
            let arr = as_f64_array(tension)?;
            let arr: PyReadonlyArray2<'_, f64> = arr.extract().map_err(|_| {
                PyValueError::new_err("tension: expected a 2D array (frames, constraints)")
            })?;
            if arr.shape()[0] != tr.n_frames() {
                return Err(PyValueError::new_err(format!(
                    "tension: expected {} rows, got {}",
                    tr.n_frames(),
                    arr.shape()[0]
                )));
            }
            tr.n_constraints = arr.shape()[1];
            tr.tension = arr.as_array().iter().copied().collect();
        }
        let event_t = event_t
            .map(|e| extract_f64s(e, "event_t"))
            .transpose()?
            .unwrap_or_default();
        let k = event_t.len();
        let event_index = event_index.unwrap_or_default();
        let empty = || vec![Vec3::ZERO; 0];
        let event_pos = match event_pos {
            Some(e) => extract_frames(e, Some(k), Some(n), "event_pos")?.2,
            None => empty(),
        };
        let event_vel = match event_vel {
            Some(e) => extract_frames(e, Some(k), Some(n), "event_vel")?.2,
            None => empty(),
        };
        if event_index.len() != k || event_pos.len() != k * n || event_vel.len() != k * n {
            return Err(PyValueError::new_err(
                "event_t, event_index, event_pos and event_vel must all be given, one entry per event",
            ));
        }
        tr.events = (0..k)
            .map(|e| crate::events::EventHit {
                event: event_index[e],
                t: event_t[e],
                pos: event_pos[e * n..(e + 1) * n].to_vec(),
                vel: event_vel[e * n..(e + 1) * n].to_vec(),
            })
            .collect();
        tr.terminated_by = terminated_by;
        let metadata = metadata.unwrap_or_else(|| PyDict::new(py)).unbind();
        Self::from_rust(py, tr, metadata)
    }

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

/// Collects frames into chunks of `chunk_size` and passes each to a Python callable.
struct ChunkSink {
    sink: Py<PyAny>,
    buffer: Trajectory,
    chunk_size: usize,
    metadata: Py<PyDict>,
}

impl ChunkSink {
    fn flush(&mut self) -> SimResult<()> {
        let mut empty = Trajectory::new(self.buffer.n_particles, self.buffer.energies);
        empty.reserve(self.chunk_size);
        let chunk = std::mem::replace(&mut self.buffer, empty);
        Python::attach(|py| -> PyResult<()> {
            let chunk = Py::new(
                py,
                PyTrajectory::from_rust(py, chunk, self.metadata.clone_ref(py))?,
            )?;
            self.sink.bind(py).call1((chunk,))?;
            Ok(())
        })
        .map_err(SimError::Python)
    }
}

impl Recorder for ChunkSink {
    fn frame(&mut self, frame: &Frame<'_>) -> SimResult<()> {
        self.buffer.frame(frame)?;
        if self.buffer.n_frames() >= self.chunk_size {
            self.flush()?;
        }
        Ok(())
    }

    fn event(&mut self, hit: crate::events::EventHit) -> SimResult<()> {
        self.buffer.event(hit)
    }
}

/// Runs `drive` with a recorder that keeps frames in memory, or with one that streams
/// them to `sink` in chunks, and returns the Trajectory. On error the exception carries the
/// frames recorded so far as its ``trajectory`` attribute.
#[allow(clippy::too_many_arguments)]
fn record_run(
    py: Python<'_>,
    n: usize,
    energies: bool,
    reserve: usize,
    sink: Option<Py<PyAny>>,
    chunk_size: usize,
    metadata: Py<PyDict>,
    drive: impl FnOnce(&mut dyn Recorder) -> SimResult<Option<usize>> + Send,
) -> PyResult<PyTrajectory> {
    let (result, trajectory) = match sink {
        None => {
            let mut traj = Trajectory::new(n, energies);
            traj.reserve(reserve);
            let result = py.detach(|| drive(&mut traj));
            (result.map(|term| traj.terminated_by = term), traj)
        }
        Some(sink) => {
            if !sink.bind(py).is_callable() || chunk_size == 0 {
                return Err(PyValueError::new_err(
                    "sink must be callable and chunk_size at least 1",
                ));
            }
            let mut buffer = Trajectory::new(n, energies);
            buffer.reserve(chunk_size);
            let mut recorder = ChunkSink {
                sink,
                buffer,
                chunk_size,
                metadata: metadata.clone_ref(py),
            };
            let result = py.detach(|| {
                let outcome = drive(&mut recorder);
                // Hand over what is buffered even if the run failed.
                if let Ok(term) = outcome {
                    recorder.buffer.terminated_by = term;
                }
                let b = &recorder.buffer;
                let flushed =
                    if b.n_frames() > 0 || !b.events.is_empty() || b.terminated_by.is_some() {
                        recorder.flush()
                    } else {
                        Ok(())
                    };
                outcome.and_then(|term| flushed.map(|()| term))
            });
            let mut summary = Trajectory::new(n, energies);
            (result.map(|term| summary.terminated_by = term), summary)
        }
    };
    let trajectory = PyTrajectory::from_rust(py, trajectory, metadata)?;
    match result {
        Ok(()) => Ok(trajectory),
        Err(error) => {
            let err = PyErr::from(error);
            if let Ok(partial) = Py::new(py, trajectory) {
                // Some exception types reject new attributes; the error still propagates.
                let _ = err.value(py).setattr("trajectory", partial);
            }
            Err(err)
        }
    }
}

// ---------------------------------------------------------------------------
// User-defined integrator schemes

fn parse_ops(ops: &[(String, f64)]) -> PyResult<Vec<Op>> {
    ops.iter()
        .map(|(kind, c)| match kind.as_str() {
            "kick" => Ok(Op::Kick(*c)),
            "drift" => Ok(Op::Drift(*c)),
            other => Err(PyValueError::new_err(format!(
                "splitting steps are (\"kick\" or \"drift\", coefficient), got {other:?}"
            ))),
        })
        .collect()
}

/// ``{"kind": "composition", "name", "order", "weights"}`` or
/// ``{"kind": "splitting", "name", "order", "ops": [["kick", c], ["drift", c], ...]}``.
fn describe_scheme<'py>(py: Python<'py>, scheme: &Scheme) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    match scheme {
        Scheme::Composition {
            name,
            order,
            weights,
        } => {
            d.set_item("kind", "composition")?;
            d.set_item("name", name)?;
            d.set_item("order", order)?;
            d.set_item("weights", weights.clone())?;
        }
        Scheme::Splitting { name, order, ops } => {
            d.set_item("kind", "splitting")?;
            d.set_item("name", name)?;
            d.set_item("order", order)?;
            let ops: Vec<(&str, f64)> = ops
                .iter()
                .map(|op| match *op {
                    Op::Kick(c) => ("kick", c),
                    Op::Drift(c) => ("drift", c),
                })
                .collect();
            d.set_item("ops", ops)?;
        }
        Scheme::Langevin {
            temperature,
            friction,
            seed,
            counter,
        } => {
            d.set_item("kind", "langevin")?;
            d.set_item("temperature", temperature)?;
            d.set_item("friction", friction)?;
            d.set_item("seed", seed)?;
            d.set_item("counter", counter)?;
        }
        Scheme::NoseHoover {
            temperature,
            tau,
            xi,
            eta,
        } => {
            d.set_item("kind", "nose_hoover")?;
            d.set_item("temperature", temperature)?;
            d.set_item("tau", tau)?;
            d.set_item("xi", xi)?;
            d.set_item("eta", eta)?;
        }
    }
    Ok(d)
}

fn scheme_from_description(d: &Bound<'_, PyDict>) -> PyResult<Scheme> {
    let field = |k: &str| -> PyResult<Bound<'_, PyAny>> {
        d.get_item(k)?
            .ok_or_else(|| PyValueError::new_err(format!("integrator scheme has no {k:?}")))
    };
    let kind: String = field("kind")?.extract()?;
    match kind.as_str() {
        "langevin" => {
            return Ok(Scheme::Langevin {
                temperature: field("temperature")?.extract()?,
                friction: field("friction")?.extract()?,
                seed: field("seed")?.extract()?,
                counter: field("counter")?.extract()?,
            })
        }
        "nose_hoover" => {
            return Ok(Scheme::NoseHoover {
                temperature: field("temperature")?.extract()?,
                tau: field("tau")?.extract()?,
                xi: field("xi")?.extract()?,
                eta: field("eta")?.extract()?,
            })
        }
        _ => {}
    }
    let name: String = field("name")?.extract()?;
    let order: u32 = field("order")?.extract()?;
    match kind.as_str() {
        "composition" => Ok(Scheme::Composition {
            name,
            order,
            weights: field("weights")?.extract()?,
        }),
        "splitting" => {
            let ops: Vec<(String, f64)> = field("ops")?.extract()?;
            Ok(Scheme::Splitting {
                name,
                order,
                ops: parse_ops(&ops)?,
            })
        }
        _ => Err(PyValueError::new_err(format!(
            "unknown integrator scheme kind {kind:?}"
        ))),
    }
}

/// ``{"restitution", "walls": [[nx, ny, nz, offset], ...], "between_particles", "max_per_step"}``.
fn describe_collisions<'py>(py: Python<'py>, c: &Collisions) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("restitution", c.restitution)?;
    let walls: Vec<[f64; 4]> = c
        .walls
        .iter()
        .map(|w| [w.normal.x, w.normal.y, w.normal.z, w.offset])
        .collect();
    d.set_item("walls", walls)?;
    d.set_item("between_particles", c.between_particles)?;
    d.set_item("max_per_step", c.max_per_step)?;
    Ok(d)
}

fn collisions_from_description(d: &Bound<'_, PyDict>) -> PyResult<Collisions> {
    let field = |k: &str| -> PyResult<Bound<'_, PyAny>> {
        d.get_item(k)?
            .ok_or_else(|| PyValueError::new_err(format!("collisions description has no {k:?}")))
    };
    let walls: Vec<[f64; 4]> = field("walls")?.extract()?;
    Ok(Collisions {
        restitution: field("restitution")?.extract()?,
        walls: walls
            .iter()
            .map(|w| Wall {
                normal: Vec3::new(w[0], w[1], w[2]),
                offset: w[3],
            })
            .collect(),
        between_particles: field("between_particles")?.extract()?,
        max_per_step: field("max_per_step")?.extract()?,
    })
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

    /// Adds a particle (optionally charged) and returns its index.
    #[pyo3(signature = (pos, vel = None, mass = 1.0, charge = 0.0, radius = 0.0))]
    fn add_particle(
        &mut self,
        pos: &Bound<'_, PyAny>,
        vel: Option<&Bound<'_, PyAny>>,
        mass: f64,
        charge: f64,
        radius: f64,
    ) -> PyResult<usize> {
        let pos = extract_vec3(pos, "pos")?;
        let vel = vel
            .map(|v| extract_vec3(v, "vel"))
            .transpose()?
            .unwrap_or(Vec3::ZERO);
        if !charge.is_finite() {
            return Err(PyValueError::new_err(format!(
                "charge must be finite, got {charge}"
            )));
        }
        if !(radius.is_finite() && radius >= 0.0) {
            return Err(PyValueError::new_err(format!(
                "radius must be finite and non-negative, got {radius}"
            )));
        }
        let i = self.inner.add_particle(pos, vel, mass)?;
        self.inner.set_charge(i, charge)?;
        self.inner.set_radius(i, radius)?;
        Ok(i)
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

    /// Adds a rigid rod (RATTLE constraint) from particle ``i`` to ``to``: another particle's
    /// index, or a fixed point ``[x, y, z]``. ``length`` defaults to the current distance.
    /// Returns the constraint's id. Needs the ``verlet`` or ``yoshida4`` integrator.
    #[pyo3(signature = (i, to, length = None))]
    fn add_rod(
        &mut self,
        i: usize,
        to: &Bound<'_, PyAny>,
        length: Option<f64>,
    ) -> PyResult<ConstraintId> {
        let anchor = match to.extract::<usize>() {
            Ok(j) => Anchor::Particle(j),
            Err(_) => Anchor::Point(extract_vec3(to, "to")?),
        };
        Ok(self.inner.add_rod(i, anchor, length)?)
    }

    fn remove_constraint(&mut self, id: ConstraintId) -> PyResult<()> {
        self.inner.remove_constraint(id)?;
        Ok(())
    }

    fn clear_constraints(&mut self) {
        self.inner.clear_constraints();
    }

    /// Constraints as ``{id: name}``.
    #[getter]
    fn constraints(&self) -> std::collections::BTreeMap<ConstraintId, String> {
        self.inner
            .constraints
            .iter()
            .map(|(id, r)| (id, r.name()))
            .collect()
    }

    /// Tension in each constraint (in ``constraints`` order) at the end of the last step:
    /// the force pulling the rod's ends together, negative when it pushes them apart.
    fn constraint_tensions<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, self.inner.constraint_tensions())
    }

    /// Relative tolerance of the constraint solver (default 1e-10): rod lengths are held to
    /// ``length * (1 ± tolerance)``.
    #[getter]
    fn constraint_tolerance(&self) -> f64 {
        self.inner.constraints.tolerance
    }
    #[setter]
    fn set_constraint_tolerance(&mut self, tol: f64) -> PyResult<()> {
        if !(tol > 0.0 && tol.is_finite()) {
            return Err(PyValueError::new_err("tolerance must be positive"));
        }
        self.inner.constraints.tolerance = tol;
        Ok(())
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
    fn integrator(&self) -> String {
        self.inner.integrator().name().to_string()
    }

    /// Switches to Langevin dynamics at ``temperature`` with friction rate ``friction``
    /// (BAOAB splitting; k_B = 1). Random numbers are reproducible from ``seed`` and restart
    /// exactly from checkpoints.
    #[pyo3(signature = (temperature, friction = 1.0, seed = 1))]
    fn use_langevin(&mut self, temperature: f64, friction: f64, seed: u64) -> PyResult<()> {
        self.inner
            .set_integrator(Box::new(integrators::Langevin::new(
                temperature,
                friction,
                seed,
            )?));
        Ok(())
    }

    /// Switches to a Nosé-Hoover thermostat at ``temperature`` with relaxation time ``tau``.
    /// ``total_energy() + thermostat_energy()`` is conserved.
    #[pyo3(signature = (temperature, tau = 1.0))]
    fn use_nose_hoover(&mut self, temperature: f64, tau: f64) -> PyResult<()> {
        self.inner
            .set_integrator(Box::new(integrators::NoseHoover::new(temperature, tau)?));
        Ok(())
    }

    /// Instantaneous temperature ``Σ m v² / (3 N)`` (k_B = 1) over unpinned, massive particles.
    fn temperature(&self) -> f64 {
        self.inner.temperature()
    }

    /// Pressure ``(2 K + W) / (3 volume)``, with ``W`` the virial of the pair potentials.
    fn pressure(&self, volume: f64) -> PyResult<f64> {
        Ok(self.inner.pressure(volume)?)
    }

    /// Energy stored in a thermostat's variables (Nosé-Hoover); zero otherwise.
    fn thermostat_energy(&self) -> f64 {
        self.inner.thermostat_energy()
    }

    /// Switches to a composition of velocity Verlet substeps of ``weights[k] * dt`` (the
    /// weights must sum to 1; symmetric weights give a time-reversible symplectic scheme).
    /// ``order`` is what you claim for it (reported in ``integrator_info``). Supports rods.
    /// Saved in checkpoints.
    #[pyo3(signature = (weights, order, name = "composition"))]
    fn use_composition(&mut self, weights: Vec<f64>, order: u32, name: &str) -> PyResult<()> {
        let scheme = Scheme::Composition {
            name: name.to_string(),
            order,
            weights,
        };
        self.inner.set_integrator(scheme.build()?);
        Ok(())
    }

    /// Switches to a general splitting: ``ops`` is a sequence of ``("kick", c)`` (``v += c dt
    /// a(x)``) and ``("drift", c)`` (``x += c dt v``) applied in order; kick and drift
    /// coefficients must each sum to 1. Saved in checkpoints.
    #[pyo3(signature = (ops, order, name = "splitting"))]
    fn use_splitting(&mut self, ops: Vec<(String, f64)>, order: u32, name: &str) -> PyResult<()> {
        let scheme = Scheme::Splitting {
            name: name.to_string(),
            order,
            ops: parse_ops(&ops)?,
        };
        self.inner.set_integrator(scheme.build()?);
        Ok(())
    }

    #[setter]
    fn set_integrator(&mut self, name: &str) -> PyResult<()> {
        self.inner.set_integrator(integrators::by_name(name)?);
        Ok(())
    }

    /// Takes ``steps`` steps of size ``dt`` while integrating the variational equations, and
    /// returns Lyapunov exponents and MEGNO as a dict:
    ///
    /// - ``exponents`` (n,): final estimates, largest first, per unit time;
    /// - ``t`` (F,) and ``running`` (F, n): estimates every ``record_every`` steps;
    /// - ``megno`` and ``mean_megno`` (F,): MEGNO ``Y(t)`` and its average ``<Y>(t)``, which
    ///   tends to 2 for quasi-periodic motion and grows like ``λ t / 2`` for chaos.
    ///
    /// The tangent vectors are integrated by the world's integrator together with the state
    /// (the world advances as in :meth:`run`) and re-orthonormalised every
    /// ``renormalize_every`` steps (Benettin's method). Not available with rods or
    /// ``wisdom_holman``.
    #[pyo3(signature = (dt, steps, n = 1, *, renormalize_every = 1, record_every = 100, seed = 1))]
    #[allow(clippy::too_many_arguments)]
    fn lyapunov<'py>(
        &mut self,
        py: Python<'py>,
        dt: f64,
        steps: usize,
        n: usize,
        renormalize_every: usize,
        record_every: usize,
        seed: u64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let options = LyapunovOptions {
            n_exponents: n,
            renormalize_every,
            record_every,
            seed,
        };
        let world = &mut self.inner;
        let run = py.detach(|| world.lyapunov(dt, steps, &options))?;
        let d = PyDict::new(py);
        let k = run.exponents.len().max(n);
        let flat: Vec<f64> = run.running.iter().flatten().copied().collect();
        d.set_item("exponents", PyArray1::from_vec(py, run.exponents))?;
        d.set_item("t", PyArray1::from_vec(py, run.t))?;
        d.set_item(
            "running",
            PyArray1::from_vec(py, flat).reshape([run.running.len(), k])?,
        )?;
        d.set_item("megno", PyArray1::from_vec(py, run.megno))?;
        d.set_item("mean_megno", PyArray1::from_vec(py, run.mean_megno))?;
        Ok(d)
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
    ///
    /// ``energies=False`` skips the per-frame kinetic and potential energy (the potential is
    /// a full force pass, O(N²) for gravity), which makes recording every step cheap.
    ///
    /// With ``sink``, recorded frames and events are not kept: every ``chunk_size`` frames
    /// they are passed to ``sink(chunk)`` as a Trajectory (e.g. a
    /// :class:`physim.TrajectoryWriter`), so memory stays flat however long the run. The
    /// returned Trajectory then holds no frames, only ``terminated_by`` and ``metadata``.
    #[pyo3(signature = (dt, steps, record_every = 1, events = None, *, energies = true, sink = None, chunk_size = 1024))]
    #[allow(clippy::too_many_arguments)]
    fn run(
        &mut self,
        py: Python<'_>,
        dt: f64,
        steps: usize,
        record_every: usize,
        events: Option<Vec<PyRef<'_, PyEvent>>>,
        energies: bool,
        sink: Option<Py<PyAny>>,
        chunk_size: usize,
    ) -> PyResult<PyTrajectory> {
        let events = events.unwrap_or_default();
        let metadata = self.run_metadata(py, energies, &events)?;
        {
            let d = metadata.bind(py);
            d.set_item("dt", dt)?;
            d.set_item("steps", steps)?;
            d.set_item("record_every", record_every)?;
        }
        let events: Vec<Event> = events.iter().map(|e| e.build(py)).collect();
        let options = RunOptions {
            record_every,
            events: &events,
            energies,
        };
        let world = &mut self.inner;
        let reserve = steps.checked_div(record_every).map_or(0, |f| f + 2);
        let n = world.state.len();
        record_run(
            py,
            n,
            energies,
            reserve,
            sink,
            chunk_size,
            metadata,
            |recorder| world.run_into(dt, steps, &options, recorder),
        )
    }

    /// Integrates to ``t_end`` with adaptive Dormand-Prince 5(4) steps, keeping the
    /// estimated local error of every position and velocity component below
    /// ``atol + rtol * |value|``. The world's own integrator is not used, and constraints
    /// are not supported. ``t_end`` may lie before the current time (integrating backwards).
    ///
    /// Frames are recorded at the initial state and every accepted step, or, if ``times``
    /// is given, exactly at those times (from the method's continuous extension, at no
    /// extra cost; they must be ordered and lie within the run). ``events``, ``energies``,
    /// ``sink`` and ``chunk_size`` work as in :meth:`run`; events are located on the
    /// continuous extension.
    ///
    /// The returned trajectory's ``metadata["adaptive"]`` holds the settings and the work
    /// done: accepted and rejected steps and force evaluations.
    ///
    /// Adaptive Runge-Kutta is not symplectic: on conservative problems the energy error
    /// grows slowly with time instead of staying bounded. It pays off when the needed step
    /// varies a lot (eccentric orbits, close encounters, transients).
    #[pyo3(signature = (t_end, *, rtol = 1e-9, atol = 1e-12, times = None, events = None, energies = true, first_step = None, max_step = None, max_steps = 10_000_000, sink = None, chunk_size = 1024))]
    #[allow(clippy::too_many_arguments)]
    fn run_adaptive(
        &mut self,
        py: Python<'_>,
        t_end: f64,
        rtol: f64,
        atol: f64,
        times: Option<Vec<f64>>,
        events: Option<Vec<PyRef<'_, PyEvent>>>,
        energies: bool,
        first_step: Option<f64>,
        max_step: Option<f64>,
        max_steps: usize,
        sink: Option<Py<PyAny>>,
        chunk_size: usize,
    ) -> PyResult<PyTrajectory> {
        let events = events.unwrap_or_default();
        let metadata = self.run_metadata(py, energies, &events)?;
        let settings = PyDict::new(py);
        settings.set_item("t_end", t_end)?;
        settings.set_item("rtol", rtol)?;
        settings.set_item("atol", atol)?;
        settings.set_item("first_step", first_step)?;
        settings.set_item("max_step", max_step)?;
        settings.set_item("max_steps", max_steps)?;
        settings.set_item("times", times.is_some())?;
        metadata.bind(py).set_item("integrator", "dopri5")?;
        metadata.bind(py).set_item("adaptive", &settings)?;
        let events: Vec<Event> = events.iter().map(|e| e.build(py)).collect();
        let has_times = times.is_some();
        let times = times.unwrap_or_default();
        let options = AdaptiveOptions {
            rtol,
            atol,
            first_step,
            max_step: max_step.unwrap_or(f64::INFINITY),
            max_steps,
            output: if has_times {
                Output::Times(&times)
            } else {
                Output::Steps
            },
            events: &events,
            energies,
        };
        let reserve = match options.output {
            Output::Times(t) => t.len(),
            Output::Steps => 0,
        };
        let world = &mut self.inner;
        let n = world.state.len();
        let mut stats = None;
        let result = record_run(
            py,
            n,
            energies,
            reserve,
            sink,
            chunk_size,
            metadata.clone_ref(py),
            |recorder| {
                let outcome = world.run_adaptive_into(t_end, &options, recorder)?;
                stats = Some(outcome.stats);
                Ok(outcome.terminated_by)
            },
        );
        if let Some(s) = stats {
            settings.set_item("accepted", s.accepted)?;
            settings.set_item("rejected", s.rejected)?;
            settings.set_item("evaluations", s.evaluations)?;
        }
        result
    }

    /// Everything needed to continue this simulation later, bit for bit, as a dict of
    /// arrays and plain values. Save it with :func:`physim.save_checkpoint`.
    /// ``CustomForce`` functions cannot be stored: they appear as
    /// ``{"type": "external", "name": ...}`` and must be passed again on restore.
    fn checkpoint<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let s = &self.inner.state;
        let d = PyDict::new(py);
        d.set_item("format", "physim-checkpoint")?;
        d.set_item("format_version", 1)?;
        d.set_item("engine_version", env!("CARGO_PKG_VERSION"))?;
        d.set_item("t", s.t)?;
        d.set_item("positions", vecs_to_array(py, &s.pos)?)?;
        d.set_item("velocities", vecs_to_array(py, &s.vel)?)?;
        d.set_item("masses", PyArray1::from_slice(py, &s.mass))?;
        d.set_item("charges", PyArray1::from_slice(py, &s.charge))?;
        d.set_item("radii", PyArray1::from_slice(py, &s.radius))?;
        d.set_item("pinned", PyArray1::from_slice(py, &s.pinned))?;
        d.set_item("integrator", self.inner.integrator().name())?;
        let scheme = self.inner.integrator().scheme();
        match &scheme {
            Some(s) => d.set_item("integrator_scheme", describe_scheme(py, s)?)?,
            None => d.set_item("integrator_scheme", py.None())?,
        }
        d.set_item("forces", describe_forces(py, &self.inner)?)?;
        d.set_item("next_force_id", self.inner.forces.next_id())?;
        let c = &self.inner.constraints;
        d.set_item("constraints", describe_constraints(py, c)?)?;
        d.set_item("next_constraint_id", c.next_id())?;
        d.set_item("constraint_tolerance", c.tolerance)?;
        d.set_item("constraint_max_iterations", c.max_iterations)?;
        match self.inner.collisions() {
            Some(c) => d.set_item("collisions", describe_collisions(py, c)?)?,
            None => d.set_item("collisions", py.None())?,
        }
        Ok(d)
    }

    /// Rebuilds a World from :meth:`checkpoint` output. ``custom_forces`` maps the name or
    /// id of each force stored as ``"external"`` to the force to use, e.g.
    /// ``{"precession": CustomForce(...)}``.
    #[staticmethod]
    #[pyo3(signature = (checkpoint, custom_forces = None))]
    fn from_checkpoint(
        py: Python<'_>,
        checkpoint: &Bound<'_, PyDict>,
        custom_forces: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Self> {
        let get = |key: &str| -> PyResult<Bound<'_, PyAny>> {
            checkpoint
                .get_item(key)?
                .ok_or_else(|| PyValueError::new_err(format!("checkpoint has no {key:?}")))
        };
        if let Some(f) = checkpoint.get_item("format")? {
            if f.extract::<String>()? != "physim-checkpoint" {
                return Err(PyValueError::new_err("not a physim checkpoint"));
            }
        }
        let positions = get("positions")?;
        let n = as_f64_array(&positions)?.len()?;
        let pinned: PyReadonlyArray1<'_, bool> = py
            .import("numpy")?
            .call_method1("asarray", (get("pinned")?, "bool"))?
            .extract()?;
        let state = State {
            t: get("t")?.extract()?,
            pos: extract_vecs(&positions, n, "positions")?,
            vel: extract_vecs(&get("velocities")?, n, "velocities")?,
            mass: extract_f64s(&get("masses")?, "masses")?,
            // Checkpoints from before charges existed have none: all neutral.
            charge: match checkpoint.get_item("charges")? {
                Some(c) => extract_f64s(&c, "charges")?,
                None => vec![0.0; n],
            },
            radius: match checkpoint.get_item("radii")? {
                Some(r) => extract_f64s(&r, "radii")?,
                None => vec![0.0; n],
            },
            pinned: pinned.as_array().to_vec(),
        };
        let mut forces = Vec::new();
        let mut external = Vec::new();
        for desc in get("forces")?.try_iter()? {
            let desc = desc?;
            let desc = desc
                .cast::<PyDict>()
                .map_err(|_| PyValueError::new_err("each force in a checkpoint must be a dict"))?;
            let id: ForceId = desc
                .get_item("id")?
                .ok_or_else(|| PyValueError::new_err("force description has no \"id\""))?
                .extract()?;
            let kind: String = desc
                .get_item("type")?
                .map(|k| k.extract())
                .transpose()?
                .unwrap_or_default();
            let saved = if kind == "external" {
                let name: String = desc
                    .get_item("name")?
                    .map(|k| k.extract())
                    .transpose()?
                    .unwrap_or_default();
                // Look the force up now, while we hold the GIL.
                let supplied = match custom_forces {
                    None => None,
                    Some(c) => match c.get_item(id)? {
                        Some(f) => Some(f),
                        None => c.get_item(&name)?,
                    },
                };
                let Some(supplied) = supplied else {
                    return Err(PyValueError::new_err(format!(
                        "the checkpoint uses force {name:?} (id {id}), which cannot be saved; \
                         pass it as custom_forces={{{name:?}: ...}}"
                    )));
                };
                external.push((id, build_force(&supplied)?));
                SavedForce::External { name }
            } else {
                SavedForce::Builtin(builtin_from_description(desc)?)
            };
            forces.push((id, saved));
        }
        let mut rods = Vec::new();
        if let Some(list) = checkpoint.get_item("constraints")? {
            for desc in list.try_iter()? {
                let desc = desc?;
                let desc = desc.cast::<PyDict>().map_err(|_| {
                    PyValueError::new_err("each constraint in a checkpoint must be a dict")
                })?;
                rods.push(rod_from_description(desc)?);
            }
        }
        let default_next = rods.iter().map(|(id, _)| id + 1).max().unwrap_or(0);
        let mut constraints = Constraints::from_parts(
            rods,
            optional(checkpoint, "next_constraint_id")?.unwrap_or(default_next),
        )?;
        if let Some(tol) = optional(checkpoint, "constraint_tolerance")? {
            constraints.tolerance = tol;
        }
        if let Some(n) = optional(checkpoint, "constraint_max_iterations")? {
            constraints.max_iterations = n;
        }
        let integrator_scheme = match checkpoint.get_item("integrator_scheme")? {
            Some(d) if !d.is_none() => Some(scheme_from_description(d.cast::<PyDict>()?)?),
            _ => None,
        };
        let collisions = match checkpoint.get_item("collisions")? {
            Some(c) if !c.is_none() => Some(collisions_from_description(c.cast::<PyDict>()?)?),
            _ => None,
        };
        let checkpoint = Checkpoint {
            state,
            integrator: get("integrator")?.extract()?,
            integrator_scheme,
            collisions,
            forces,
            next_force_id: get("next_force_id")?.extract()?,
            constraints,
        };
        let inner = World::from_checkpoint(checkpoint, |id, _| {
            let k = external
                .iter()
                .position(|(i, _)| *i == id)
                .expect("supplied above");
            Ok(external.swap_remove(k).1)
        })?;
        Ok(Self { inner })
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

    /// Electric charge of every particle, shape (N,).
    #[getter]
    fn charges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.inner.state.charge)
    }
    #[setter]
    fn set_charges(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let arr = as_f64_array(value)?;
        let arr: PyReadonlyArray1<'_, f64> = arr
            .extract()
            .map_err(|_| PyValueError::new_err("charges: expected a 1D array"))?;
        Ok(self.inner.set_charges(arr.as_array().to_vec())?)
    }

    /// Enables hard collisions: particles with a radius bounce off each other (if
    /// ``between_particles``) and off ``walls`` with coefficient of restitution
    /// ``restitution``. Each wall is ``(normal, offset)``: the plane ``normal·x = offset``,
    /// with particles kept on the side the unit ``normal`` points to (see
    /// :func:`physim.box_walls`). Contacts are found and resolved event by event inside each
    /// step, exactly so when no forces act between collisions. Not available with
    /// :meth:`run_adaptive` or :meth:`lyapunov`.
    #[pyo3(signature = (restitution = 1.0, walls = None, *, between_particles = true, max_per_step = 100_000))]
    fn set_collisions(
        &mut self,
        restitution: f64,
        walls: Option<Vec<(Bound<'_, PyAny>, f64)>>,
        between_particles: bool,
        max_per_step: usize,
    ) -> PyResult<()> {
        let walls = walls
            .unwrap_or_default()
            .iter()
            .map(|(n, offset)| {
                Ok(Wall {
                    normal: extract_vec3(n, "wall normal")?,
                    offset: *offset,
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
        self.inner.set_collisions(Some(Collisions {
            restitution,
            walls,
            between_particles,
            max_per_step,
        }))?;
        Ok(())
    }

    /// Disables hard collisions.
    fn clear_collisions(&mut self) -> PyResult<()> {
        self.inner.set_collisions(None)?;
        Ok(())
    }

    /// The hard-collision settings as a dict, or None.
    #[getter]
    fn collisions<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        self.inner
            .collisions()
            .map(|c| describe_collisions(py, c))
            .transpose()
    }

    /// Number of hard collisions resolved so far.
    #[getter]
    fn collision_count(&self) -> u64 {
        self.inner.collision_count()
    }

    /// Radius of every particle (for collisions and contact forces), shape (N,).
    #[getter]
    fn radii<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &self.inner.state.radius)
    }
    #[setter]
    fn set_radii(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let arr = as_f64_array(value)?;
        let arr: PyReadonlyArray1<'_, f64> = arr
            .extract()
            .map_err(|_| PyValueError::new_err("radii: expected a 1D array"))?;
        Ok(self.inner.set_radii(arr.as_array().to_vec())?)
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

impl PyWorld {
    /// Describes a run for ``Trajectory.metadata``.
    fn run_metadata(
        &self,
        py: Python<'_>,
        energies: bool,
        events: &[PyRef<'_, PyEvent>],
    ) -> PyResult<Py<PyDict>> {
        let w = &self.inner;
        let d = PyDict::new(py);
        d.set_item("engine_version", env!("CARGO_PKG_VERSION"))?;
        d.set_item("integrator", w.integrator().name())?;
        d.set_item("energies", energies)?;
        d.set_item("t0", w.state.t)?;
        d.set_item("n_particles", w.state.len())?;
        d.set_item("masses", w.state.mass.clone())?;
        d.set_item("charges", w.state.charge.clone())?;
        d.set_item("radii", w.state.radius.clone())?;
        d.set_item("pinned", w.state.pinned.clone())?;
        d.set_item("forces", describe_forces(py, w)?)?;
        d.set_item("constraints", describe_constraints(py, &w.constraints)?)?;
        let names: Vec<String> = events.iter().map(|e| e.name()).collect();
        d.set_item("events", names)?;
        Ok(d.unbind())
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
    m.add_class::<PyDampedSpring>()?;
    m.add_class::<PyModulatedSpring>()?;
    m.add_class::<PySpringNetwork>()?;
    m.add_class::<PyPowerLaw>()?;
    m.add_class::<PyYukawa>()?;
    m.add_class::<PyPlummerPotential>()?;
    m.add_class::<PyHernquistPotential>()?;
    m.add_class::<PyHarmonicTrap>()?;
    m.add_class::<PyPeriodicForce>()?;
    m.add_class::<PyPostNewtonian>()?;
    m.add_class::<PyJ2Oblateness>()?;
    m.add_class::<PyHenonHeiles>()?;
    m.add_class::<PyElectricField>()?;
    m.add_class::<PyMagneticField>()?;
    m.add_class::<PyCoulomb>()?;
    m.add_class::<PyTreeGravity>()?;
    m.add_class::<PySoftContact>()?;
    m.add_class::<PyLennardJones>()?;
    m.add_class::<PyMorse>()?;
    m.add_class::<PyTabulatedPair>()?;
    m.add_class::<PyFieldForce>()?;
    m.add_class::<PyCustomForce>()?;
    m.add("INTEGRATORS", integrators::NAMES.to_vec())?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
