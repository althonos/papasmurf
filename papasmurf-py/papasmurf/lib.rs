use std::sync::Arc;

use pyo3::exceptions::PyRuntimeError;
use pyo3::exceptions::PyValueError;
use pyo3::exceptions::PyIndexError;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::types::PyString;

mod error;

use self::error::Error;

// --- Builder -----------------------------------------------------------------

/// A builder to generate a database of reference k-mers from 16S genes.
#[pyclass(module = "papasmurf.lib")]
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
        let mut p = Vec::new();
        for result in primers.iter()? {
            let item = result?;
            if item.len()? != 2 {
                return Err(PyValueError::new_err("expected pair of strings"));
            }
            let forward = item.get_item(0)?.downcast::<PyString>()?;
            let backward = item.get_item(1)?.downcast::<PyString>()?;
            let f = papasmurf::Primer::new(forward.to_str()?).map_err(Error::from)?;
            let b = papasmurf::Primer::new(backward.to_str()?).map_err(Error::from)?;
            p.push(papasmurf::Paired::new(f, b))
        }
        Ok(Self {
            builder: papasmurf::Builder::new(p),
        }
        .into())
    }

    /// Add a new sequence to the builder, extracting k-mer regions.
    pub fn add<'py>(&self, name: &'py PyString, sequence: &'py PyString) -> PyResult<()> {
        let name_ = name.to_str()?;
        let seq_ = sequence.to_str()?;
        match self.builder.add(name_, seq_) {
            Ok(_) => Ok(()),
            Err(e) => Err(Error::from(e).into()),
        }
    }

    /// Build and index the database from the k-mers stored in the builder.
    pub fn to_database(&self) -> PyResult<Database> {
        Ok(Database::from(self.builder.to_database()))
    }
}

// --- Database ----------------------------------------------------------------

/// A database, storing k-mer regions extracted from reference organisms.
#[pyclass(module = "papasmurf.lib")]
#[derive(Debug)]
pub struct Database {
    db: Arc<papasmurf::Database>,
    #[pyo3(get)]
    names: DatabaseNames,
}

impl From<papasmurf::Database> for Database {
    fn from(db: papasmurf::Database) -> Self {
        let db = Arc::new(db);
        let names = DatabaseNames::new(db.clone());
        Self { db, names }
    }
}

#[pymethods]
impl Database {
    /// Load a database serialized at the given path.
    #[staticmethod]
    #[pyo3(signature = (filename, format = "messagepack"))]
    pub fn load<'py>(filename: &'py PyString, format: &str) -> PyResult<Self> {
        let name = filename.to_str()?;
        let f = std::fs::File::open(name)
            .map(std::io::BufReader::new)
            .map_err(|e| Error::Io(e, name.to_string()))?;
        match format {
            "json" => match serde_json::from_reader::<_, papasmurf::Database>(f) {
                Ok(db) => Ok(Database::from(db)),
                Err(e) => Err(Error::from(e).into()),
            },
            "messagepack" => match rmp_serde::from_read::<_, papasmurf::Database>(f) {
                Ok(database) => Ok(Database::from(database)),
                Err(e) => Err(PyRuntimeError::new_err(e.to_string())),
            },
            _ => Err(PyValueError::new_err(format!(
                "invalid format: {:?}",
                format
            ))),
        }
    }

    /// Store the database to the given path.
    #[pyo3(signature = (filename, format = "messagepack"))]
    pub fn dump<'py>(&self, filename: &'py PyString, format: &str) -> PyResult<()> {
        let name = filename.to_str()?;
        let mut f = match std::fs::File::create(name) {
            Ok(file) => std::io::BufWriter::new(file),
            Err(e) => return Err(Error::Io(e, name.to_string()).into()),
        };
        match format {
            "json" => match serde_json::to_writer(f, self.db.as_ref()) {
                Ok(_) => Ok(()),
                Err(e) => Err(Error::from(e).into()),
            },
            "messagepack" => match rmp_serde::encode::write(&mut f, self.db.as_ref()) {
                Ok(_) => Ok(()),
                Err(e) => Err(PyRuntimeError::new_err(e.to_string())),
            },
            _ => Err(PyValueError::new_err(format!(
                "invalid format: {:?}",
                format
            ))),
        }
    }
}

#[pyclass(module = "papasmurf.lib")]
#[derive(Debug, Clone)]
pub struct DatabaseNames {
    db: Arc<papasmurf::Database>,
}

impl DatabaseNames {
    pub fn new(db: Arc<papasmurf::Database>) -> Self {
        Self { db }
    }
}

#[pymethods]
impl DatabaseNames {
    pub fn __len__(&self) -> usize {
        self.db.names().len()
    }

    pub fn __getitem__(&self, i: usize) -> PyResult<PyObject> {

        let names = self.db.names();
        let mut i_ = i as isize;

        if i_ < 0 {
            i_ += names.len() as isize;
        }
        if i_ < 0 || i_ >= names.len() as isize {
            return Err(PyIndexError::new_err("list index out of range"));
        }

        let name = &names[i_ as usize];
        Ok(Python::with_gil(|py| PyString::new(py, &*name).to_object(py)))
    }
}


// --- Mapper ------------------------------------------------------------------

#[pyclass(module = "papasmurf.lib")]
#[derive(Debug)]
pub struct Mapper {
    mapper: papasmurf::Mapper<Arc<papasmurf::Database>>,
}

#[pymethods]
impl Mapper {
    /// Create a new mapper for the given database.
    #[new]
    #[pyo3(signature = (database, *, primer_mismatches=2, kmer_mismatches=2, error_probability=0.005, partial_hits=false))]
    pub fn __init__<'py>(
        database: &'py Database,
        primer_mismatches: usize,
        kmer_mismatches: usize,
        error_probability: f32,
        partial_hits: bool,
    ) -> PyResult<PyClassInitializer<Self>> {
        let db = database.db.clone();
        let mapper = papasmurf::Mapper::new(db)
            .with_primer_mismatches(primer_mismatches)
            .with_kmer_mismatches(kmer_mismatches)
            .with_error_probability(error_probability)
            .with_partial_hits(partial_hits);
        Ok(Self { mapper }.into())
    }

    /// Add a new read to the mapper.
    pub fn add<'py>(&self, forward: &'py PyString, backward: &'py PyString) -> PyResult<bool> {
        let py = forward.py();
        let read = papasmurf::Paired::new(forward.to_str()?, backward.to_str()?);
        py.allow_threads(|| self.mapper.add(read))
            .map_err(|e| Error::from(e).into())
    }

    /// Finish mapping and get the mapper results.
    ///
    /// The mapper is reset and can be used to map a new sample after calling
    /// this method.
    pub fn finish(&mut self) -> PyResult<MapperResult> {
        let db = AsRef::<Arc<papasmurf::Database>>::as_ref(&self.mapper).clone();
        let mapper = std::mem::replace(&mut self.mapper, papasmurf::Mapper::new(db));
        let result = mapper.finish();
        Ok(MapperResult::from(result))
    }
}

// --- MapperResult ------------------------------------------------------------

#[pyclass(module = "papasmurf.lib")]
#[derive(Debug)]
pub struct MapperResult {
    result: papasmurf::MapperResult<Arc<papasmurf::Database>>,
    frequencies: Option<PyObject>,
    proportions: Option<PyObject>,
}

impl From<papasmurf::MapperResult<Arc<papasmurf::Database>>> for MapperResult {
    fn from(result: papasmurf::MapperResult<Arc<papasmurf::Database>>) -> Self {
        Self {
            result,
            frequencies: None,
            proportions: None,
        }
    }
}

#[pymethods]
impl MapperResult {
    #[pyo3(signature = (n = 1))]
    pub fn refine(&mut self, n: usize) -> PyResult<()> {
        self.frequencies = None;
        self.proportions = None;
        for i in 0..n {
            self.result.refine();
        }
        Ok(())
    }

    #[getter]
    pub fn names(&self) -> DatabaseNames {
        let db: &Arc<papasmurf::Database> = self.result.as_ref();
        DatabaseNames::new(db.clone())
    }

    #[getter]
    pub fn frequencies(&mut self) -> PyResult<PyObject> {
        if let Some(freq) = &self.frequencies {
            return Ok(freq.clone());
        }
        let a = Python::with_gil(|py| {
            let f = py.allow_threads(|| self.result.frequencies());
            let l = PyList::new(py, f);
            py.import(intern!(py, "array"))?
                .call_method1(intern!(py, "array"), (intern!(py, "f"), l))
                .map(|a| a.to_object(py))
        })?;
        self.frequencies = Some(a.clone());
        Ok(a)
    }

    #[getter]
    pub fn proportions(&mut self) -> PyResult<PyObject> {
        if let Some(prop) = &self.proportions {
            return Ok(prop.clone());
        }
        let a = Python::with_gil(|py| {
            let p = py.allow_threads(|| self.result.proportions());
            let l = PyList::new(py, p);
            py.import(intern!(py, "array"))?
                .call_method1(intern!(py, "array"), (intern!(py, "f"), l))
                .map(|a| a.to_object(py))
        })?;
        self.proportions = Some(a.clone());
        Ok(a)
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
    m.add_class::<Mapper>()?;
    m.add_class::<MapperResult>()?;

    Ok(())
}
