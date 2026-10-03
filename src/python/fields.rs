//! Python bindings for grid fields: the wave, heat and Poisson solvers and particle-mesh
//! gravity.

use numpy::{PyArray1, PyArrayDyn, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::{as_f64_array, builtin_repr, extract_vec3, PyForceSpec};
use crate::fields::pm::ParticleMesh;
use crate::fields::solvers::{poisson, Heat, HeatMethod, Wave};
use crate::fields::{Boundary, Grid};
use crate::forces::Force;
use crate::vec3::Vec3;

/// A float64 array's flat values and shape.
fn flat(obj: &Bound<'_, PyAny>, what: &str) -> PyResult<(Vec<f64>, Vec<usize>)> {
    let arr = as_f64_array(obj)?;
    let shape: Vec<usize> = arr.getattr("shape")?.extract()?;
    if shape.is_empty() || shape.len() > 3 {
        return Err(PyValueError::new_err(format!(
            "{what}: expected a 1, 2 or 3 dimensional array, got shape {shape:?}"
        )));
    }
    let values: Vec<f64> = arr
        .call_method0("ravel")?
        .call_method0("tolist")?
        .extract()?;
    Ok((values, shape))
}

fn to_array<'py>(
    py: Python<'py>,
    values: &[f64],
    grid: &Grid,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    PyArray1::from_slice(py, values).reshape(grid.shape())
}

fn grid_from(shape: Vec<usize>, spacing: f64, boundary: &str) -> PyResult<Grid> {
    Ok(Grid::new(&shape, spacing, Boundary::parse(boundary)?)?)
}

/// 1D coordinate arrays along each axis (cell centres for Neumann grids).
fn axes<'py>(py: Python<'py>, grid: &Grid) -> Vec<Bound<'py, PyArray1<f64>>> {
    let shift = if grid.boundary == Boundary::Neumann {
        0.5
    } else {
        0.0
    };
    grid.shape()
        .iter()
        .map(|&n| PyArray1::from_vec(py, (0..n).map(|i| (i as f64 + shift) * grid.h).collect()))
        .collect()
}

fn set_field(grid: &Grid, obj: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<f64>> {
    let (values, shape) = flat(obj, what)?;
    if shape != grid.shape() {
        return Err(PyValueError::new_err(format!(
            "{what}: expected shape {:?}, got {shape:?}",
            grid.shape()
        )));
    }
    Ok(values)
}

/// The wave equation ``u_tt = c² ∇²u`` on a grid of ``shape`` (1 to 3 sizes) with spacing
/// ``spacing``, stepped by leapfrog (second order in space and time, stable for
/// ``c dt / h ≤ 1/√d``).
///
/// ``boundary`` is ``"dirichlet"`` (the edge points keep their values: a clamped string or
/// membrane), ``"neumann"`` (free edges; points are cell centres) or ``"periodic"``.
/// Set ``u`` and ``v`` (displacement and velocity arrays of ``shape``), then ``step``.
#[pyclass(name = "WaveEquation", module = "physim")]
struct PyWave(Wave);

#[pymethods]
impl PyWave {
    #[new]
    #[pyo3(signature = (shape, spacing, c = 1.0, boundary = "dirichlet"))]
    fn new(shape: Vec<usize>, spacing: f64, c: f64, boundary: &str) -> PyResult<Self> {
        Ok(Self(Wave::new(grid_from(shape, spacing, boundary)?, c)?))
    }
    /// Displacement field; assignable.
    #[getter]
    fn u<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        to_array(py, &self.0.u, &self.0.grid)
    }
    #[setter]
    fn set_u(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let u = set_field(&self.0.grid, value, "u")?;
        let v = self.0.v.clone();
        Ok(self.0.set_state(u, v)?)
    }
    /// Velocity field ``u_t``; assignable.
    #[getter]
    fn v<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        to_array(py, &self.0.v, &self.0.grid)
    }
    #[setter]
    fn set_v(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let v = set_field(&self.0.grid, value, "v")?;
        let u = self.0.u.clone();
        Ok(self.0.set_state(u, v)?)
    }
    /// Current time; assignable.
    #[getter]
    fn t(&self) -> f64 {
        self.0.t
    }
    #[setter]
    fn set_t(&mut self, t: f64) {
        self.0.t = t;
    }
    /// Coordinates along each axis, a list of 1D arrays.
    #[getter]
    fn axes<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyArray1<f64>>> {
        axes(py, &self.0.grid)
    }
    /// Largest stable time step ``h / (c √d)``.
    #[getter]
    fn max_stable_dt(&self) -> f64 {
        self.0.max_stable_dt()
    }
    /// Takes ``n`` steps of ``dt``.
    #[pyo3(signature = (dt, n = 1))]
    fn step(&mut self, dt: f64, n: usize) -> PyResult<()> {
        Ok(self.0.step(dt, n)?)
    }
    /// Discrete energy ``½ ∫ (u_t² + c² |∇u|²)``, conserved to O(dt²) without drift.
    fn energy(&self) -> f64 {
        self.0.energy()
    }
    fn __repr__(&self) -> String {
        format!(
            "WaveEquation(shape={:?}, spacing={}, c={}, boundary={:?}, t={})",
            self.0.grid.shape(),
            self.0.grid.h,
            self.0.c,
            self.0.grid.boundary.name(),
            self.0.t
        )
    }
}

/// The heat (diffusion) equation ``u_t = D ∇²u`` on a grid of ``shape`` with spacing
/// ``spacing``. ``method="crank_nicolson"`` (default) is second order in time and stable for
/// any step (conjugate-gradient solve each step); ``"explicit"`` is first order and needs
/// ``dt ≤ h² / (2 d D)``. Boundaries as for :class:`WaveEquation`; with ``"neumann"`` or
/// ``"periodic"`` the total ``Σ u h^d`` is conserved.
#[pyclass(name = "HeatEquation", module = "physim")]
struct PyHeat(Heat);

#[pymethods]
impl PyHeat {
    #[new]
    #[pyo3(signature = (shape, spacing, diffusivity = 1.0, boundary = "neumann", method = "crank_nicolson"))]
    fn new(
        shape: Vec<usize>,
        spacing: f64,
        diffusivity: f64,
        boundary: &str,
        method: &str,
    ) -> PyResult<Self> {
        let method = match method {
            "crank_nicolson" => HeatMethod::CrankNicolson,
            "explicit" => HeatMethod::Explicit,
            m => {
                return Err(PyValueError::new_err(format!(
                    "method must be 'crank_nicolson' or 'explicit', got {m:?}"
                )))
            }
        };
        Ok(Self(Heat::new(
            grid_from(shape, spacing, boundary)?,
            diffusivity,
            method,
        )?))
    }
    /// The field; assignable.
    #[getter]
    fn u<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        to_array(py, &self.0.u, &self.0.grid)
    }
    #[setter]
    fn set_u(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let u = set_field(&self.0.grid, value, "u")?;
        Ok(self.0.set_u(u)?)
    }
    /// Current time; assignable.
    #[getter]
    fn t(&self) -> f64 {
        self.0.t
    }
    #[setter]
    fn set_t(&mut self, t: f64) {
        self.0.t = t;
    }
    /// Coordinates along each axis, a list of 1D arrays.
    #[getter]
    fn axes<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyArray1<f64>>> {
        axes(py, &self.0.grid)
    }
    /// Largest stable explicit step ``h² / (2 d D)``.
    #[getter]
    fn max_explicit_dt(&self) -> f64 {
        self.0.max_explicit_dt()
    }
    /// Conjugate-gradient iterations used by the last Crank-Nicolson step.
    #[getter]
    fn last_iterations(&self) -> usize {
        self.0.last_iterations
    }
    /// Takes ``n`` steps of ``dt``.
    #[pyo3(signature = (dt, n = 1))]
    fn step(&mut self, dt: f64, n: usize) -> PyResult<()> {
        Ok(self.0.step(dt, n)?)
    }
    /// ``Σ u h^d``, conserved with Neumann and periodic boundaries.
    fn total(&self) -> f64 {
        self.0.total()
    }
    fn __repr__(&self) -> String {
        format!(
            "HeatEquation(shape={:?}, spacing={}, diffusivity={}, boundary={:?}, t={})",
            self.0.grid.shape(),
            self.0.grid.h,
            self.0.diffusivity,
            self.0.grid.boundary.name(),
            self.0.t
        )
    }
}

/// Solves ``∇²φ = f`` for an array ``f`` (1 to 3 dimensions) with grid spacing ``spacing``.
///
/// - ``"periodic"`` (every size a power of two): FFT, exact for the discrete Laplacian; the
///   mean of ``f`` is removed and ``φ`` has zero mean.
/// - ``"dirichlet"``: conjugate gradients; ``phi`` (optional, same shape) supplies the
///   boundary values on the edge points (zero otherwise).
/// - ``"neumann"`` (cell-centred, zero normal derivative): conjugate gradients; the mean of
///   ``f`` is removed and ``φ`` has zero mean.
#[pyfunction]
#[pyo3(signature = (f, spacing, boundary = "periodic", phi = None))]
fn solve_poisson<'py>(
    py: Python<'py>,
    f: &Bound<'py, PyAny>,
    spacing: f64,
    boundary: &str,
    phi: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
    let (values, shape) = flat(f, "f")?;
    let grid = grid_from(shape, spacing, boundary)?;
    let mut out = match phi {
        Some(p) => set_field(&grid, p, "phi")?,
        None => vec![0.0; grid.len()],
    };
    poisson(&grid, &values, &mut out)?;
    to_array(py, &out, &grid)
}

/// Particle-mesh gravity: the masses are spread on a ``cells³`` grid (cloud-in-cell), the field
/// is found with FFTs, and accelerations are interpolated back. Cost O(N + M log M) for M grid
/// points, independent of how the particles are arranged.
///
/// - ``periodic=False`` (default): isolated (open) boundaries via a zero-padded grid; matches
///   direct summation with Plummer softening ``softening`` (default: one cell) to ~1% beyond
///   about 8 cells. Every particle must stay inside the cube of side ``box_size`` around
///   ``center``.
/// - ``periodic=True``: the box repeats in every direction (positions are wrapped);
///   ``∇²φ = 4πG(ρ - ρ̄)``. No potential energy is reported.
///
/// Momentum is conserved to round-off. ``cells`` is a power of two (memory and time grow
/// as ``cells³``; 32-64 is typical).
#[pyclass(frozen, name = "ParticleMesh", module = "physim")]
pub(super) struct PyParticleMesh(pub(super) ParticleMesh);

#[pymethods]
impl PyParticleMesh {
    #[new]
    #[pyo3(signature = (box_size, cells = 64, G = 1.0, center = None, periodic = false, softening = None))]
    #[allow(non_snake_case)]
    fn new(
        box_size: f64,
        cells: usize,
        G: f64,
        center: Option<&Bound<'_, PyAny>>,
        periodic: bool,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        let center = match center {
            Some(c) => extract_vec3(c, "center")?,
            None => Vec3::ZERO,
        };
        let pm = ParticleMesh::new(G, cells, box_size, center, periodic, softening);
        pm.validate()?;
        Ok(Self(pm))
    }
    /// Grid spacing ``box_size / cells``.
    #[getter]
    fn spacing(&self) -> f64 {
        self.0.spacing()
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        builtin_repr(py, &self.0)
    }
}

impl PyForceSpec for PyParticleMesh {
    fn build(&self, _py: Python<'_>) -> Box<dyn Force> {
        Box::new(self.0.clone())
    }
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyWave>()?;
    m.add_class::<PyHeat>()?;
    m.add_class::<PyParticleMesh>()?;
    m.add_function(wrap_pyfunction!(solve_poisson, m)?)?;
    Ok(())
}
