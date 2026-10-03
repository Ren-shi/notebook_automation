use std::fmt;

/// Errors raised while building or stepping a simulation.
#[derive(Debug)]
pub enum SimError {
    /// Bad input: unknown integrator, out-of-range particle index, wrong array shape, ...
    Invalid(String),
    /// An exception raised by a user-supplied Python callback, kept intact so the
    /// original exception type and traceback reach the Python caller.
    #[cfg(feature = "python")]
    Python(pyo3::PyErr),
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SimError::Invalid(msg) => write!(f, "{msg}"),
            #[cfg(feature = "python")]
            SimError::Python(err) => write!(f, "Python callback failed: {err}"),
        }
    }
}

impl std::error::Error for SimError {}

pub type Result<T> = std::result::Result<T, SimError>;

pub(crate) fn invalid<T>(msg: impl Into<String>) -> Result<T> {
    Err(SimError::Invalid(msg.into()))
}
