use pyo3::exceptions::PyEOFError;

use pyo3::exceptions::PyOSError;
use pyo3::exceptions::PyRuntimeError;

use pyo3::exceptions::PyValueError;
use pyo3::PyErr;

// ---------------------------------------------------------------------------

#[macro_export]
macro_rules! raise(
    ($py:expr, $error_type:ident ($msg:expr) from $inner:expr ) => ({
        let err = $error_type::new_err($msg).to_object($py);
        err.call_method1(
            $py,
            "__setattr__",
            ("__cause__".to_object($py), $inner.to_object($py)),
        )?;
        return Err(PyErr::from_value(err.as_ref($py)))
    })
);

/// A wrapper to convert all errors from the different libraries into a `PyErr`.
#[derive(Debug)]
pub enum Error {
    Papasmurf(papasmurf::Error),
    Io(std::io::Error, String),
    SerdeJson(serde_json::error::Error),
}

impl From<Error> for PyErr {
    fn from(error: Error) -> PyErr {
        match error {
            // PAPASMURF value errors
            Error::Papasmurf(papasmurf::Error::InvalidDna) => {
                PyValueError::new_err("invalid DNA symbols")
            }
            Error::Papasmurf(papasmurf::Error::InvalidDimensions) => {
                PyValueError::new_err("invalid dimensions")
            }
            // I/O errors
            Error::Io(io_error, path) => {
                if let Some(n) = io_error.raw_os_error() {
                    PyOSError::new_err((n, path))
                } else {
                    PyRuntimeError::new_err(io_error.to_string())
                }
            }
            // Serde JSON error
            Error::SerdeJson(err) => {
                if err.is_io() {
                    let io_error: std::io::Error = err.into();
                    if let Some(n) = io_error.raw_os_error() {
                        PyOSError::new_err((n, io_error.to_string()))
                    } else {
                        PyRuntimeError::new_err(io_error.to_string())
                    }
                } else if err.is_eof() {
                    PyEOFError::new_err(err.to_string())
                } else {
                    PyValueError::new_err(err.to_string())
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------

impl From<papasmurf::Error> for Error {
    fn from(error: papasmurf::Error) -> Self {
        Error::Papasmurf(error)
    }
}

impl From<serde_json::error::Error> for Error {
    fn from(error: serde_json::error::Error) -> Self {
        Error::SerdeJson(error)
    }
}
