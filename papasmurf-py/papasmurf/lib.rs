use std::io::ErrorKind;
use std::sync::Arc;

use pyo3::exceptions::PyBufferError;
use pyo3::exceptions::PyFileNotFoundError;
use pyo3::exceptions::PyIndexError;
use pyo3::exceptions::PyIsADirectoryError;
use pyo3::exceptions::PyOSError;
use pyo3::exceptions::PyRuntimeError;
use pyo3::exceptions::PyTypeError;
use pyo3::exceptions::PyValueError;
use pyo3::ffi::Py_ssize_t;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::types::PyList;
use pyo3::types::PyString;
use pyo3::AsPyPointer;

// --- Builder -----------------------------------------------------------------

/// A builder to generate a database of reference k-mers from 16S genes.
#[pyclass(module = "lightmotif.lib")]
#[derive(Debug)]
pub struct Builder {
    builder: papasmurf::Builder,
}

impl From<papasmurf::Builder> for Builder {
    fn from(builder: papasmurf::Builder) -> Self {
        Self { builder }
    }
}

impl From<Builder> for papasmurf::Builder {
    fn from(builder: Builder) -> Self {
        builder.builder
    }
}

#[pymethods]
impl Builder {
    /// Create a new database builder with the given parameters.
    #[new]
    pub fn __init__<'py>(primers: &'py PyAny) -> PyResult<PyClassInitializer<Self>> {
        unimplemented!()
    }

    /// Add a new sequence to the builder, extracting k-mer regions.
    pub fn add<'py>(&self, id: &'py PyString, sequence: &'py PyString) -> PyResult<()> {
        unimplemented!()
    }

    /// Build and index the database from the k-mers stored in the builder.
    pub fn to_database(&self) -> PyResult<Database> {
        unimplemented!()
    }
}

// --- Database ----------------------------------------------------------------

/// A database, storing k-mer regions extracted from reference organisms.
#[pyclass(module = "lightmotif.lib")]
#[derive(Debug)]
pub struct Database {
    db: Arc<papasmurf::Database>,
}

impl From<papasmurf::Database> for Database {
    fn from(db: papasmurf::Database) -> Self {
        Self { db: Arc::new(db) }
    }
}

#[pymethods]
impl Database {
    /// Load a database serialized at the given path.
    #[staticmethod]
    pub fn load<'py>(filename: &'py PyString) -> PyResult<Self> {
        let name = filename.to_str()?;
        let f = match std::fs::File::open(name) {
            Ok(file) => std::io::BufReader::new(file),
            Err(e) => {
                if let Some(n) = e.raw_os_error() {
                    return Err(PyOSError::new_err((n, name.to_string())));
                } else {
                    return Err(PyRuntimeError::new_err(e.to_string()));
                }
            }
        };
        match serde_json::from_reader::<_, papasmurf::Database>(f) {
            Ok(database) => Ok(Database::from(database)),
            Err(e) => return Err(PyValueError::new_err(e.to_string())),
        }
    }

    /// Store the database to the given path.
    pub fn dump<'py>(&self, filename: &'py PyString) -> PyResult<()> {
        unimplemented!()
    }
}

// --- Mapper ------------------------------------------------------------------

#[pyclass(module = "lightmotif.lib")]
#[derive(Debug)]
pub struct Mapper {
    mapper: papasmurf::Mapper<Arc<papasmurf::Database>>,
}

#[pymethods]
impl Mapper {
    /// Create a new mapper for the given database.
    #[new]
    pub fn __init__<'py>(database: &'py Database) -> PyResult<PyClassInitializer<Self>> {
        let db = database.db.clone();
        let mapper = papasmurf::Mapper::new(db);
        Ok(Self { mapper }.into())
    }

    /// Add a new read to the mapper.
    pub fn add<'py>(&self, forward: &'py PyString, backward: &'py PyString) -> PyResult<()> {
        unimplemented!()
    }

    /// Finish mapping and get the mapper results.
    pub fn finish(&mut self) {
        let db = self.mapper.as_ref().clone();
        let mapper = std::mem::replace(&mut self.mapper, papasmurf::Mapper::new(db));
        unreachable!()
    }
}

/// PyO3 bindings to ``papasmurf``, a library for 16S multiple region analysis.
#[pymodule]
#[pyo3(name = "lib")]
pub fn init(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__package__", "papasmurf")?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS").replace(':', "\n"))?;

    m.add_class::<Database>()?;
    m.add_class::<Builder>()?;

    Ok(())
}
