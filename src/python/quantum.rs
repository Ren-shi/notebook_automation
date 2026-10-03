//! Python bindings for the Schrödinger equation solver.

use numpy::{
    Complex64, PyArray1, PyArrayDyn, PyArrayMethods, PyReadonlyArrayDyn, PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use super::as_f64_array;
use super::fields::{flat, grid_from, set_field, to_array};
use crate::error::SimError;
use crate::fields::quantum::{QuantumMethod, Schrodinger};
use crate::fields::{Boundary, Grid};

/// `(t, psi)` recorded by [`PySchrodinger::run`].
type Frames<'py> = (Bound<'py, PyArray1<f64>>, Bound<'py, PyArrayDyn<Complex64>>);
/// `(energies, states)` from [`PySchrodinger::eigenstates`].
type Eigen<'py> = (Bound<'py, PyArray1<f64>>, Bound<'py, PyArrayDyn<f64>>);

fn complex_array<'py>(
    py: Python<'py>,
    values: &[(f64, f64)],
    shape: &[usize],
) -> PyResult<Bound<'py, PyArrayDyn<Complex64>>> {
    let data: Vec<Complex64> = values
        .iter()
        .map(|&(re, im)| Complex64::new(re, im))
        .collect();
    PyArray1::from_vec(py, data).reshape(shape.to_vec())
}

/// Per-axis values from a scalar or a sequence with one entry per grid axis.
fn per_axis(obj: &Bound<'_, PyAny>, dims: usize, what: &str) -> PyResult<[f64; 3]> {
    let (values, _) = match flat(obj, what) {
        Ok(v) => v,
        Err(_) => (vec![obj.extract::<f64>()?], vec![1]),
    };
    let mut out = [0.0; 3];
    match values.len() {
        1 => out[..dims].iter_mut().for_each(|o| *o = values[0]),
        n if n == dims => out[..dims].copy_from_slice(&values),
        n => {
            return Err(PyValueError::new_err(format!(
                "{what}: expected a number or {dims} values (one per axis), got {n}"
            )))
        }
    }
    Ok(out)
}

/// Default origin: the grid is centred on zero.
fn centred_origin(grid: &Grid) -> [f64; 3] {
    let mut origin = [0.0; 3];
    for (a, o) in origin.iter_mut().enumerate().take(grid.dims()) {
        let n = grid.n[a] as f64;
        *o = match grid.boundary {
            Boundary::Dirichlet => -0.5 * (n - 1.0) * grid.h,
            // Periodic: [-L/2, L/2) with a point at zero; Neumann points are cell centres.
            _ => -0.5 * n * grid.h,
        };
    }
    origin
}

/// The time-dependent Schrödinger equation ``iħ ∂ψ/∂t = [-ħ²/(2m) ∇² + V(x, t) - iW(x)] ψ``
/// for one particle on a grid of ``shape`` (1 to 3 sizes) with spacing ``spacing``.
///
/// - ``method="split_step"`` (the default on periodic grids; sizes must be powers of two):
///   Strang splitting with the kinetic step exact in Fourier space, second order in time, or
///   fourth order with ``order=4``. Exactly unitary.
/// - ``method="crank_nicolson"`` (the default otherwise; any boundary): second order in space
///   and time, unitary; on ``"dirichlet"`` grids ψ vanishes on the edge points (hard walls).
///
/// The grid is centred on zero unless ``origin`` (the coordinates of the first point) is
/// given. ``potential`` is an array of ``shape`` or a function ``V(t, x[, y[, z]])`` of
/// coordinate arrays (as from ``numpy.meshgrid(*axes, indexing="ij")``) returning an array
/// broadcastable to ``shape``. ``absorber`` (``W ≥ 0``, see :meth:`absorbing_layer`) removes
/// outgoing waves. ψ is normalised so that ``Σ |ψ|² h^d = 1``.
///
/// Split-step is accurate for any step on smooth potentials; with sharp ones (steps,
/// barriers) keep ``dt ≲ 2 m h² / (π ħ)`` so the shortest waves the grid holds are resolved
/// in time, or grid-scale noise appears.
#[pyclass(name = "Schrodinger", module = "physim")]
struct PySchrodinger {
    inner: Schrodinger,
    /// The potential function, if any (kept to return it from the getter).
    potential_fn: Option<Py<PyAny>>,
}

impl PySchrodinger {
    fn shape(&self) -> Vec<usize> {
        self.inner.grid.shape()
    }

    fn axes_vec(&self) -> Vec<Vec<f64>> {
        let g = &self.inner.grid;
        let shift = if g.boundary == Boundary::Neumann {
            0.5
        } else {
            0.0
        };
        (0..g.dims())
            .map(|a| {
                (0..g.n[a])
                    .map(|i| self.inner.origin[a] + (i as f64 + shift) * g.h)
                    .collect()
            })
            .collect()
    }

    fn vector<'py>(&self, py: Python<'py>, v: [f64; 3]) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_slice(py, &v[..self.inner.grid.dims()])
    }

    fn set_potential_function(&mut self, py: Python<'_>, f: &Bound<'_, PyAny>) -> PyResult<()> {
        let np = py.import("numpy")?;
        let axes: Vec<Bound<'_, PyArray1<f64>>> = self
            .axes_vec()
            .into_iter()
            .map(|a| PyArray1::from_vec(py, a))
            .collect();
        let kwargs = pyo3::types::PyDict::new(py);
        kwargs.set_item("indexing", "ij")?;
        let mesh = np.call_method("meshgrid", PyTuple::new(py, axes)?, Some(&kwargs))?;
        let coords: Vec<Py<PyAny>> = mesh
            .try_iter()?
            .map(|c| c.map(|c| c.unbind()))
            .collect::<PyResult<_>>()?;
        let func = f.clone().unbind();
        let shape = self.shape();
        let callback = move |t: f64, v: &mut [f64]| {
            Python::attach(|py| -> PyResult<()> {
                let mut args: Vec<Bound<'_, PyAny>> = vec![t.into_pyobject(py)?.into_any()];
                args.extend(coords.iter().map(|c| c.bind(py).clone()));
                let out = func.bind(py).call1(PyTuple::new(py, args)?)?;
                let np = py.import("numpy")?;
                let arr = np.call_method1(
                    "broadcast_to",
                    (as_f64_array(&out)?, PyTuple::new(py, &shape)?),
                )?;
                let arr: PyReadonlyArrayDyn<'_, f64> = arr.extract()?;
                for (vi, x) in v.iter_mut().zip(arr.as_array().iter()) {
                    *vi = *x;
                }
                Ok(())
            })
            .map_err(SimError::Python)
        };
        self.inner.set_potential_fn(Box::new(callback))?;
        self.potential_fn = Some(f.clone().unbind());
        Ok(())
    }
}

#[pymethods]
impl PySchrodinger {
    #[new]
    #[pyo3(signature = (shape, spacing, hbar = 1.0, mass = 1.0, boundary = "periodic", method = None, order = 2, origin = None, potential = None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        shape: Vec<usize>,
        spacing: f64,
        hbar: f64,
        mass: f64,
        boundary: &str,
        method: Option<&str>,
        order: usize,
        origin: Option<&Bound<'_, PyAny>>,
        potential: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let grid = grid_from(shape, spacing, boundary)?;
        let method = match method {
            Some(m) => QuantumMethod::parse(m)?,
            None if grid.boundary == Boundary::Periodic
                && grid.n.iter().all(|&s| s.is_power_of_two()) =>
            {
                QuantumMethod::SplitStep
            }
            None => QuantumMethod::CrankNicolson,
        };
        let mut inner = Schrodinger::new(grid, hbar, mass, method)?;
        inner.set_order(order)?;
        inner.origin = match origin {
            Some(o) => per_axis(o, grid.dims(), "origin")?,
            None => centred_origin(&grid),
        };
        let mut s = Self {
            inner,
            potential_fn: None,
        };
        if let Some(p) = potential {
            s.set_potential(py, p)?;
        }
        Ok(s)
    }

    /// The wavefunction (complex array of ``shape``); assignable (values on Dirichlet edge
    /// points are set to zero). Not normalised automatically: see :meth:`normalize`.
    #[getter]
    fn psi<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<Complex64>>> {
        complex_array(py, &self.inner.psi, &self.shape())
    }
    #[setter]
    fn set_psi(&mut self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let arr = py
            .import("numpy")?
            .call_method1("asarray", (value, "complex128"))?;
        let arr: PyReadonlyArrayDyn<'_, Complex64> = arr.extract()?;
        if arr.shape() != self.shape().as_slice() {
            return Err(PyValueError::new_err(format!(
                "psi: expected shape {:?}, got {:?}",
                self.shape(),
                arr.shape()
            )));
        }
        let psi = arr.as_array().iter().map(|c| (c.re, c.im)).collect();
        Ok(self.inner.set_psi(psi)?)
    }
    /// Probability density ``|ψ|²``.
    #[getter]
    fn density<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        let d: Vec<f64> = self
            .inner
            .psi
            .iter()
            .map(|(re, im)| re * re + im * im)
            .collect();
        to_array(py, &d, &self.inner.grid)
    }
    /// The potential at the current time (an array of ``shape``). Assign an array, a
    /// function ``V(t, x[, y[, z]])``, or ``None`` for zero.
    #[getter]
    fn potential<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        to_array(py, self.inner.potential(), &self.inner.grid)
    }
    #[setter]
    fn set_potential(&mut self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        if value.is_none() {
            self.potential_fn = None;
            return Ok(self.inner.set_potential(vec![0.0; self.inner.grid.len()])?);
        }
        if value.is_callable() {
            return self.set_potential_function(py, value);
        }
        let v = set_field(&self.inner.grid, value, "potential")?;
        self.inner.set_potential(v)?;
        self.potential_fn = None;
        Ok(())
    }
    /// The potential function, or ``None`` for a static potential.
    #[getter]
    fn potential_function(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.potential_fn.as_ref().map(|f| f.clone_ref(py))
    }
    /// The absorbing potential ``W ≥ 0`` (zero: none); assignable.
    #[getter]
    fn absorber<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArrayDyn<f64>>> {
        to_array(py, self.inner.absorber(), &self.inner.grid)
    }
    #[setter]
    fn set_absorber(&mut self, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let w = set_field(&self.inner.grid, value, "absorber")?;
        Ok(self.inner.set_absorber(w)?)
    }
    /// Current time; assignable.
    #[getter]
    fn t(&self) -> f64 {
        self.inner.t
    }
    #[setter]
    fn set_t(&mut self, t: f64) -> PyResult<()> {
        Ok(self.inner.set_time(t)?)
    }
    /// Split-step composition order (2 or 4); assignable.
    #[getter]
    fn order(&self) -> usize {
        self.inner.order()
    }
    #[setter]
    fn set_order(&mut self, order: usize) -> PyResult<()> {
        Ok(self.inner.set_order(order)?)
    }
    #[getter]
    fn method(&self) -> &'static str {
        self.inner.method.name()
    }
    #[getter]
    fn boundary(&self) -> &'static str {
        self.inner.grid.boundary.name()
    }
    #[getter(shape)]
    fn grid_shape(&self) -> Vec<usize> {
        self.shape()
    }
    #[getter]
    fn spacing(&self) -> f64 {
        self.inner.grid.h
    }
    #[getter]
    fn hbar(&self) -> f64 {
        self.inner.hbar
    }
    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass
    }
    /// Volume element ``h^d``.
    #[getter]
    fn cell_volume(&self) -> f64 {
        self.inner.cell_volume()
    }
    /// Coordinates along each axis, a list of 1D arrays.
    #[getter]
    fn axes<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyArray1<f64>>> {
        self.axes_vec()
            .into_iter()
            .map(|a| PyArray1::from_vec(py, a))
            .collect()
    }
    /// Linear-solver iterations of the last Crank-Nicolson step (0 for the direct 1D solve).
    #[getter]
    fn last_iterations(&self) -> usize {
        self.inner.last_iterations
    }
    /// Sets ψ to a normalised Gaussian wave packet ``Π exp(-(x - c)²/4σ²) e^{ip·x/ħ}`` with
    /// centre ``center``, position spread ``width`` (σ) and mean momentum ``momentum``; each
    /// a number or one value per axis (defaults: 0, 1 and 0).
    #[pyo3(signature = (center = None, width = None, momentum = None))]
    fn set_gaussian(
        &mut self,
        center: Option<&Bound<'_, PyAny>>,
        width: Option<&Bound<'_, PyAny>>,
        momentum: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let d = self.inner.grid.dims();
        let get = |v: Option<&Bound<'_, PyAny>>, default: f64, what: &str| match v {
            Some(v) => per_axis(v, d, what),
            None => Ok([default; 3]),
        };
        Ok(self.inner.set_gaussian(
            get(center, 0.0, "center")?,
            get(width, 1.0, "width")?,
            get(momentum, 0.0, "momentum")?,
        )?)
    }
    /// Sets the absorber to ``strength ((width - d) / width)²`` within ``width`` of the
    /// grid's edges (``d``: distance to the nearest edge). A layer several wavelengths wide with
    /// ``strength`` comparable to the packet's kinetic energy absorbs with little reflection.
    fn absorbing_layer(&mut self, width: f64, strength: f64) -> PyResult<()> {
        Ok(self.inner.set_absorbing_layer(width, strength)?)
    }
    /// Takes ``n`` steps of ``dt``. If a step fails, ψ and ``t`` stay at the last completed step.
    #[pyo3(signature = (dt, n = 1))]
    fn step(&mut self, dt: f64, n: usize) -> PyResult<()> {
        Ok(self.inner.step(dt, n)?)
    }
    /// Takes ``steps`` steps of ``dt``, recording ψ every ``record_every`` steps (and at the
    /// start). Returns ``(t, psi)``: times ``(F,)`` and wavefunctions ``(F, *shape)``.
    #[pyo3(signature = (dt, steps, record_every = 1))]
    fn run<'py>(
        &mut self,
        py: Python<'py>,
        dt: f64,
        steps: usize,
        record_every: usize,
    ) -> PyResult<Frames<'py>> {
        if record_every == 0 {
            return Err(PyValueError::new_err("record_every must be at least 1"));
        }
        let mut times = vec![self.inner.t];
        let mut frames = self.inner.psi.clone();
        let mut done = 0;
        while done < steps {
            let k = record_every.min(steps - done);
            self.inner.step(dt, k)?;
            done += k;
            times.push(self.inner.t);
            frames.extend_from_slice(&self.inner.psi);
        }
        let mut shape = vec![times.len()];
        shape.extend(self.shape());
        Ok((
            PyArray1::from_vec(py, times),
            complex_array(py, &frames, &shape)?,
        ))
    }
    /// Scales ψ to unit norm.
    fn normalize(&mut self) -> PyResult<()> {
        Ok(self.inner.normalize()?)
    }
    /// ``Σ |ψ|² h^d``: 1 for a normalised state; decreases as an absorber removes probability.
    fn norm(&self) -> f64 {
        self.inner.norm()
    }
    /// ``⟨x⟩`` along each axis (normalised by the current norm).
    fn position<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self.vector(py, self.inner.position()?))
    }
    /// ``⟨x²⟩ - ⟨x⟩²`` along each axis.
    fn position_variance<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self.vector(py, self.inner.position_variance()?))
    }
    /// ``⟨p⟩`` along each axis (spectral for split-step, central differences otherwise).
    fn momentum<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self.vector(py, self.inner.momentum()?))
    }
    /// ``⟨T⟩`` with the method's kinetic operator.
    fn kinetic_energy(&self) -> PyResult<f64> {
        Ok(self.inner.kinetic_energy()?)
    }
    /// ``⟨V⟩`` at the current time.
    fn potential_energy(&self) -> PyResult<f64> {
        Ok(self.inner.potential_energy()?)
    }
    /// ``⟨H⟩ = ⟨T⟩ + ⟨V⟩`` (absorber excluded). Conserved for a static potential: exactly by
    /// Crank-Nicolson, to O(dt²) (or O(dt⁴)) without drift by split-step.
    fn energy(&self) -> PyResult<f64> {
        Ok(self.inner.energy()?)
    }
    /// The ``count`` lowest eigenstates of ``H`` (potential at the current time, absorber
    /// ignored), of the same discretisation the time stepping uses. Returns
    /// ``(energies, states)``: energies ascending ``(count,)`` and real, normalised states
    /// ``(count, *shape)``. Converged when every residual ``‖Hψ - Eψ‖`` is below
    /// ``tol (E - V_min)``.
    #[pyo3(signature = (count = 1, tol = 1e-9, max_iter = 1000))]
    fn eigenstates<'py>(
        &mut self,
        py: Python<'py>,
        count: usize,
        tol: f64,
        max_iter: usize,
    ) -> PyResult<Eigen<'py>> {
        let (energies, states) = self.inner.eigenstates(count, tol, max_iter)?;
        let mut shape = vec![count];
        shape.extend(self.shape());
        let flat: Vec<f64> = states.into_iter().flatten().collect();
        Ok((
            PyArray1::from_vec(py, energies),
            PyArray1::from_vec(py, flat).reshape(shape)?,
        ))
    }
    fn __repr__(&self) -> String {
        format!(
            "Schrodinger(shape={:?}, spacing={}, hbar={}, mass={}, boundary={:?}, method={:?}, t={})",
            self.shape(),
            self.inner.grid.h,
            self.inner.hbar,
            self.inner.mass,
            self.inner.grid.boundary.name(),
            self.inner.method.name(),
            self.inner.t
        )
    }
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySchrodinger>()?;
    Ok(())
}
