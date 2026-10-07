//! Python bindings for the nuclear planner's per-event engine (`physim.nuclear.events` builds the
//! configuration from a setup file; these functions are its private back end).

use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::nuclear::{Channel, Face, Generator, Layer, Shape, Table, TwoBody};
use crate::vec3::Vec3;

fn item<'py>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<Bound<'py, PyAny>> {
    d.get_item(key)?
        .ok_or_else(|| PyValueError::new_err(format!("event generator config: missing '{key}'")))
}

fn f64_of(d: &Bound<'_, PyDict>, key: &str) -> PyResult<f64> {
    item(d, key)?.extract()
}

fn usize_of(d: &Bound<'_, PyDict>, key: &str) -> PyResult<usize> {
    item(d, key)?.extract()
}

fn vec_of(d: &Bound<'_, PyDict>, key: &str) -> PyResult<Vec<f64>> {
    item(d, key)?.extract()
}

fn vec3_of(d: &Bound<'_, PyDict>, key: &str) -> PyResult<Vec3> {
    let v = vec_of(d, key)?;
    match v.as_slice() {
        [x, y, z] => Ok(Vec3::new(*x, *y, *z)),
        _ => Err(PyValueError::new_err(format!(
            "'{key}' must have 3 components"
        ))),
    }
}

fn dicts<'py>(d: &Bound<'py, PyDict>, key: &str) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let list = item(d, key)?.cast_into::<PyList>()?;
    list.iter().map(|x| Ok(x.cast_into::<PyDict>()?)).collect()
}

fn face(d: &Bound<'_, PyDict>) -> PyResult<Face> {
    let shape: String = item(d, "shape")?.extract()?;
    let shape = match shape.as_str() {
        "rectangle" => Shape::Rectangle {
            width: f64_of(d, "width")?,
            height: f64_of(d, "height")?,
            strips_x: item(d, "strips_x")?.extract()?,
            strips_y: item(d, "strips_y")?.extract()?,
        },
        "annular" | "circle" => Shape::Annular {
            inner: f64_of(d, "inner_radius")?,
            outer: f64_of(d, "outer_radius")?,
            rings: item(d, "rings")?.extract()?,
            sectors: item(d, "sectors")?.extract()?,
        },
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown detector shape '{other}'"
            )))
        }
    };
    Ok(Face {
        centre: vec3_of(d, "centre")?,
        n: vec3_of(d, "n")?,
        u: vec3_of(d, "u")?,
        v: vec3_of(d, "v")?,
        shape,
        material: usize_of(d, "material")?,
        dead_layer: f64_of(d, "dead_layer")?,
        thickness: f64_of(d, "thickness")?,
        fwhm: f64_of(d, "fwhm")?,
        threshold: f64_of(d, "threshold")?,
    })
}

fn table(d: &Bound<'_, PyDict>) -> PyResult<Table> {
    Table::new(
        &vec_of(d, "energy")?,
        &vec_of(d, "range")?,
        &vec_of(d, "stopping")?,
        &vec_of(d, "w")?,
    )
    .map_err(PyValueError::new_err)
}

fn generator(cfg: &Bound<'_, PyDict>) -> PyResult<Generator> {
    let layers = dicts(cfg, "layers")?
        .iter()
        .map(|l| {
            Ok(Layer {
                thickness: f64_of(l, "thickness")?,
                material: usize_of(l, "material")?,
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    let channels = dicts(cfg, "channels")?
        .iter()
        .map(|c| {
            Ok(Channel {
                layer: usize_of(c, "layer")?,
                target: usize_of(c, "target")?,
                atoms_per_cm2: f64_of(c, "atoms_per_cm2")?,
                k: f64_of(c, "k")?,
                u_min: f64_of(c, "u_min")?,
                u_max: f64_of(c, "u_max")?,
                probability: f64_of(c, "probability")?,
                excitation: match c.get_item("excitation")? {
                    Some(v) => v.extract()?,
                    None => 0.0,
                },
                excite_recoil: match c.get_item("excite_recoil")? {
                    Some(v) => v.extract()?,
                    None => true,
                },
                p_table: match c.get_item("p_table")? {
                    Some(v) => v.extract()?,
                    None => Vec::new(),
                },
            })
        })
        .collect::<PyResult<Vec<_>>>()?;
    Ok(Generator {
        masses: vec_of(cfg, "masses")?,
        beam: usize_of(cfg, "beam")?,
        beam_energy: f64_of(cfg, "beam_energy")?,
        energy_sigma: f64_of(cfg, "energy_sigma")?,
        spot_sigma: f64_of(cfg, "spot_sigma")?,
        particles_per_second: f64_of(cfg, "particles_per_second")?,
        tilt: f64_of(cfg, "tilt")?,
        layers,
        channels,
        materials: usize_of(cfg, "materials")?,
        tables: dicts(cfg, "tables")?
            .iter()
            .map(table)
            .collect::<PyResult<Vec<_>>>()?,
        faces: dicts(cfg, "faces")?
            .iter()
            .map(face)
            .collect::<PyResult<Vec<_>>>()?,
        max_path_factor: f64_of(cfg, "max_path_factor")?,
        skip_misses: match cfg.get_item("skip_misses")? {
            Some(v) => v.extract::<bool>()?,
            None => false,
        },
    })
}

/// Runs the event generator described by ``config`` (built by :mod:`physim.nuclear.events`):
/// events ``first .. first + n`` of a run of ``total`` events with ``seed``. Returns a dict of
/// arrays, one entry per particle that reached a detector face.
#[pyfunction]
#[pyo3(signature = (config, n, seed, total = None, first = 0))]
fn nuclear_events<'py>(
    py: Python<'py>,
    config: &Bound<'py, PyDict>,
    n: u64,
    seed: u64,
    total: Option<u64>,
    first: u64,
) -> PyResult<Bound<'py, PyDict>> {
    let gen = generator(config)?;
    let records = py
        .detach(|| gen.run(n, total.unwrap_or(n), seed, first))
        .map_err(PyValueError::new_err)?;
    let out = PyDict::new(py);
    macro_rules! column {
        ($name:literal, $t:ty, $f:expr) => {
            let v: Vec<$t> = records.iter().map($f).collect();
            out.set_item($name, PyArray1::from_vec(py, v))?;
        };
    }
    column!("event", u64, |r| r.event);
    column!("detector", u32, |r| r.detector);
    column!("segment_i", u32, |r| r.segment.0);
    column!("segment_j", u32, |r| r.segment.1);
    column!("recoil", bool, |r| r.recoil);
    column!("channel", u32, |r| r.channel);
    column!("depth", f64, |r| r.depth);
    column!("beam_energy", f64, |r| r.beam_energy);
    column!("energy", f64, |r| r.energy);
    column!("energy_face", f64, |r| r.energy_face);
    column!("deposited", f64, |r| r.deposited);
    column!("measured", f64, |r| r.measured);
    column!("counted", bool, |r| r.counted);
    column!("theta", f64, |r| r.theta.to_degrees());
    column!("phi", f64, |r| r.phi.to_degrees());
    column!("theta_cm", f64, |r| r.theta_cm.to_degrees());
    column!("weight", f64, |r| r.weight);
    Ok(out)
}

/// Two float64 arrays.
type ArrayPair<'py> = (Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>);

/// Lab angles (degrees) and kinetic energies (MeV) from the Rust two-body kinematics, for
/// checking it against :class:`physim.nuclear.kinematics.TwoBody`.
#[pyfunction]
#[pyo3(signature = (masses, beam_energy, theta_cm, recoil = false))]
fn nuclear_two_body<'py>(
    py: Python<'py>,
    masses: [f64; 4],
    beam_energy: f64,
    theta_cm: PyReadonlyArray1<'py, f64>,
    recoil: bool,
) -> PyResult<ArrayPair<'py>> {
    let [m1, m2, m3, m4] = masses;
    let tb = TwoBody::new(m1, m2, m3, m4, beam_energy)
        .ok_or_else(|| PyValueError::new_err("below threshold"))?;
    let (th, e): (Vec<f64>, Vec<f64>) = theta_cm
        .as_slice()?
        .iter()
        .map(|t| {
            let p = tb.at_cm(t.to_radians(), recoil);
            (p.theta.to_degrees(), p.energy)
        })
        .unzip();
    Ok((PyArray1::from_vec(py, th), PyArray1::from_vec(py, e)))
}

/// Mean energy after ``path`` (mg/cm²) and the straggling σ added on the way, MeV, from a
/// transport table (a dict with ``energy``, ``range``, ``stopping``, ``w``).
#[pyfunction]
fn nuclear_cross(table_dict: &Bound<'_, PyDict>, energy: f64, path: f64) -> PyResult<(f64, f64)> {
    Ok(table(table_dict)?.cross(energy, path))
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(nuclear_events, m)?)?;
    m.add_function(wrap_pyfunction!(nuclear_two_body, m)?)?;
    m.add_function(wrap_pyfunction!(nuclear_cross, m)?)?;
    Ok(())
}
