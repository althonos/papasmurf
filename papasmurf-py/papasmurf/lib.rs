use pyo3::exceptions::PyBufferError;
use pyo3::exceptions::PyIndexError;
use pyo3::exceptions::PyTypeError;
use pyo3::exceptions::PyValueError;
use pyo3::ffi::Py_ssize_t;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::types::PyList;
use pyo3::types::PyString;
use pyo3::AsPyPointer;

#[pyfunction]
pub fn add(left: usize, right: usize) -> usize {
    left + right
}

/// PyO3 bindings to ``papasmurf``, a library for 16S multiple region analysis.
#[pymodule]
#[pyo3(name = "lib")]
pub fn init(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__package__", "papasmurf")?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS").replace(':', "\n"))?;

    // m.add_function(wrap_pyfunction!(create, m)?)?;
    m.add_function(wrap_pyfunction!(add, m)?)?;

    Ok(())
}

