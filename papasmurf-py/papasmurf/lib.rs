extern crate papasmurf;
#[macro_use]
extern crate pyo3_built;
extern crate pyo3;

use std::io::Read;
use std::io::Write;
use std::sync::Arc;

use pyo3::exceptions::PyIndexError;
use pyo3::exceptions::PyRuntimeError;
use pyo3::exceptions::PyValueError;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::types::PyString;
use pyo3::types::PyTuple;

mod error;
mod pyfile;

use self::error::Error;

#[allow(dead_code)]
mod build {
    include!(concat!(env!("OUT_DIR"), "/built.rs"));
}

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
    ///
    /// Arguments:
    ///     primers (iterable of `str` pairs): The primers to use to
    ///         extract the k-mers from the reference sequences. *Both
    ///         primer sequences must be given in 5'-3' direction*.
    ///
    /// Raises:
    ///     `ValueError`: When given a primer that is not a valid
    ///         IUPAC DNA sequence.
    ///
    /// Example:
    ///     >>> builder = papasmurf.Builder([
    ///     ...     ("TGGCGGACGGGTGAGTAA", "CTGCTGCCTCCCGTAGGA"),
    ///     ...     ("TCCTACGGGAGGCAGCAG", "TATTACCGCGGCTGCTGG"),
    ///     ...     ("CAGCAGCCGCGGTAATAC", "CGCATTTCACCGCTACAC"),
    ///     ...     ("AGGATTAGATACCCTGGT", "GAATTAAACCACATGCTC"),
    ///     ...     ("GCACAAGCGGTGGAGCAT", "CGCTCGTTGCGGGACTTA"),
    ///     ...     ("AGGAAGGTGGGGATGACG", "CCCGGGAACGTATTCACC"),
    ///     ... ])
    ///
    #[new]
    pub fn __new__<'py>(primers: &Bound<'py, PyAny>) -> PyResult<PyClassInitializer<Self>> {
        let mut p = Vec::new();
        for result in primers.try_iter()? {
            let item = result?;
            if item.len()? != 2 {
                return Err(PyValueError::new_err("expected pair of strings"));
            }
            let forward = item.get_item(0)?;
            let backward = item.get_item(1)?;
            let f = papasmurf::Primer::new(forward.cast::<PyString>()?.to_str()?)
                .map_err(Error::from)?;
            let b = papasmurf::Primer::new(backward.cast::<PyString>()?.to_str()?)
                .map_err(Error::from)?;
            p.push(papasmurf::Paired::new(f, b))
        }
        Ok(Self {
            builder: papasmurf::Builder::new(p),
        }
        .into())
    }

    /// Add a new sequence to the builder, extracting k-mer regions.
    ///
    /// Arguments:
    ///     name (`str`): The name of the reference bacterium being added.
    ///     sequence (`str`): The 16S gene sequence of the reference
    ///         bacterium being added.
    ///
    /// Raises:
    ///     `ValueError`: When given a sequence that is not a valid IUPAC
    ///         DNA sequence.
    ///
    pub fn add<'py>(
        slf: PyRef<'py, Self>,
        name: &Bound<'py, PyString>,
        sequence: &Bound<'py, PyString>,
    ) -> PyResult<()> {
        let name_ = name.to_str()?;
        let seq_ = sequence.to_str()?;
        let py = slf.py();
        let builder = &slf.builder;
        match py.detach(|| builder.add(name_, seq_)) {
            Ok(_) => Ok(()),
            Err(e) => Err(Error::from(e).into()),
        }
    }

    /// Build and index the database from the k-mers stored in the builder.
    pub fn to_database<'py>(slf: PyRef<'py, Self>) -> PyResult<Database> {
        let py = slf.py();
        let builder = &slf.builder;
        Ok(Database::from(py.detach(|| builder.to_database())))
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
    /// Load a database serialized at the given file.
    #[staticmethod]
    #[pyo3(signature = (file, format = "messagepack"))]
    pub fn load<'py>(file: &Bound<'py, PyAny>, format: &str) -> PyResult<Self> {
        let f: Box<dyn Read> = if let Ok(name) = file.cast::<PyString>() {
            std::fs::File::open(name.to_str()?)
                .map(std::io::BufReader::new)
                .map_err(|e| Error::Io(e, name.to_string()))
                .map(Box::new)?
        } else {
            pyfile::PyFileRead::from_ref(file.clone())
                .map(std::io::BufReader::new)
                .map(Box::new)?
        };
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

    /// Store the database to the given file.
    #[pyo3(signature = (file, format = "messagepack"))]
    pub fn dump<'py>(
        slf: PyRef<'py, Self>,
        file: &Bound<'py, PyAny>,
        format: &str,
    ) -> PyResult<()> {
        let mut f: Box<dyn Write> = if let Ok(name) = file.cast::<PyString>() {
            std::fs::File::open(name.to_str()?)
                .map(std::io::BufWriter::new)
                .map_err(|e| Error::Io(e, name.to_string()))
                .map(Box::new)?
        } else {
            pyfile::PyFileWrite::from_ref(file.clone())
                .map(std::io::BufWriter::new)
                .map(Box::new)?
        };
        match format {
            "json" => match serde_json::to_writer(f, slf.db.as_ref()) {
                Ok(_) => Ok(()),
                Err(e) => Err(Error::from(e).into()),
            },
            "messagepack" => match rmp_serde::encode::write(&mut f, slf.db.as_ref()) {
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

/// An immutable view over the names of the reference bacteria in a database.
#[pyclass(module = "papasmurf.lib", from_py_object)]
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
    pub fn __len__<'py>(slf: PyRef<'py, Self>) -> usize {
        slf.db.names().len()
    }

    pub fn __getitem__<'py>(slf: PyRef<'py, Self>, i: usize) -> PyResult<Bound<'py, PyString>> {
        let py = slf.py();
        let names = slf.db.names();
        let mut i_ = i as isize;

        if i_ < 0 {
            i_ += names.len() as isize;
        }
        if i_ < 0 || i_ >= names.len() as isize {
            return Err(PyIndexError::new_err("list index out of range"));
        }

        let name = &names[i_ as usize];
        Ok(PyString::new(py, &*name))
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
    ///
    /// Arguments:
    ///     database (`~papasmurf.Database`): The database against which to
    ///         map the 16S sequencing reads.
    ///
    /// Keyword Arguments:
    ///     primer_mismatches (`int`): The maximum number of allowed
    ///         mismatches between the forward or backward primer and
    ///         the read sequence.
    ///     kmer_mismatches (`int`): The maximum number of allowed
    ///         mismatches between the forwards or backward database k-mer
    ///         and the read sequence.
    ///     error_probability (`float`): The *a priori* error probability
    ///         per nucleotide to use to compute the probability of origin
    ///         for each read.
    ///     partial_hits (`bool`): Whether or not to enable partial read
    ///         matching for reads shorter than the k-mers.
    ///
    #[new]
    #[pyo3(signature = (database, *, primer_mismatches=2, kmer_mismatches=2, error_probability=0.005, partial_hits=false))]
    pub fn __new__<'py>(
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

    /// Add a new paired read to the mapper.
    ///
    /// Arguments:
    ///     forward (`str`): The forward read to add to the mapper.
    ///     backward (`str`): The backward read to add to the mapper.
    ///
    /// Returns:
    ///     `bool`: Whether the read passed quality filtering and was
    ///     mapped to any database region.
    ///     
    pub fn add<'py>(
        slf: PyRef<'py, Self>,
        forward: &Bound<'py, PyString>,
        backward: &Bound<'py, PyString>,
    ) -> PyResult<bool> {
        let read = papasmurf::Paired::new(forward.to_str()?, backward.to_str()?);
        let py = slf.py();
        let mapper = &slf.mapper;
        py.detach(|| mapper.add(read))
            .map_err(|e| Error::from(e).into())
    }

    /// Finish mapping and get the mapper results.
    ///
    /// The mapper is reset and can be used to map a new sample after calling
    /// this method.
    ///
    /// Returns:
    ///     `~papasmurf.MapperResult`: The result of the read mapping, which
    ///     can be further refined to approximate the unknown read proportion
    ///     vector.
    ///
    pub fn finish<'py>(mut slf: PyRefMut<'py, Self>) -> PyResult<MapperResult> {
        let py = slf.py();
        let db = AsRef::<Arc<papasmurf::Database>>::as_ref(&slf.mapper).clone();
        let mapper = std::mem::replace(&mut slf.mapper, papasmurf::Mapper::new(db));
        let result = py.detach(|| mapper.finish());
        Ok(MapperResult::from(result))
    }
}

// --- MapperResult ------------------------------------------------------------

#[pyclass(module = "papasmurf.lib")]
#[derive(Debug)]
pub struct MapperResult {
    result: papasmurf::MapperResult<Arc<papasmurf::Database>>,
    frequencies: Option<Py<PyAny>>,
    proportions: Option<Py<PyAny>>,
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
    /// sequence of `str`: The name of each bacterium in the database.
    #[getter]
    pub fn names<'py>(slf: PyRef<'py, Self>) -> DatabaseNames {
        let db: &Arc<papasmurf::Database> = slf.result.as_ref();
        DatabaseNames::new(db.clone())
    }

    /// `array` of `float`: The bacterium frequency vector, :math:`X`.
    #[getter]
    pub fn frequencies<'py>(mut slf: PyRefMut<'py, Self>) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        if let Some(freq) = &slf.frequencies {
            return Ok(freq.clone_ref(py));
        }
        slf.frequencies = None;
        let result = &slf.result;
        let a = {
            let f = py.detach(|| result.frequencies());
            let l = PyList::new(py, f)?;
            py.import(intern!(py, "array"))?
                .call_method1(intern!(py, "array"), (intern!(py, "f"), l))
                .map(|a| a.into_pyobject(py))??
                .unbind()
        };
        slf.frequencies = Some(a.clone_ref(py));
        Ok(a)
    }

    /// `array` of `float`: The read proportion vector, :math:`\pi`.
    #[getter]
    pub fn proportions<'py>(mut slf: PyRefMut<'py, Self>) -> PyResult<Py<PyAny>> {
        let py = slf.py();
        if let Some(prop) = &slf.proportions {
            return Ok(prop.clone_ref(py));
        }
        slf.proportions = None;
        let result = &slf.result;
        let a = {
            let p = py.detach(|| result.proportions());
            let l = PyList::new(py, p)?;
            py.import(intern!(py, "array"))?
                .call_method1(intern!(py, "array"), (intern!(py, "f"), l))
                .map(|a| a.into_pyobject(py))??
                .unbind()
        };
        slf.proportions = Some(a.clone_ref(py));
        Ok(a)
    }

    /// `tuple` of `int`: The number of reads assigned to each region.
    ///
    /// A read is assigned to a region when the primer for this region had
    /// the highest score for the read. It may still fail to pass quality
    /// control (based on the `Mapper` parameters).
    #[getter]
    pub fn assigned_by_region<'py>(slf: PyRef<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let assigned = slf.result.assigned_by_region();
        PyTuple::new(slf.py(), assigned)
    }

    /// `tuple` of `int`: The number of reads mapped to each region.
    ///     
    /// A read is mapped to a region when it was mapped to any database
    /// k-mer of the region it was assigned to by primer-matching, after
    /// passing quality control.
    #[getter]
    pub fn mapped_by_region<'py>(slf: PyRef<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let mapped = slf.result.mapped_by_region();
        PyTuple::new(slf.py(), mapped)
    }

    /// `array` of `int`: The number of reads mapped to each bacterium.
    #[getter]
    pub fn mapped_by_bacterium<'py>(slf: PyRef<'py, Self>) -> PyResult<Bound<'py, PyAny>> {
        let py = slf.py();
        let result = &slf.result;
        let m = py.detach(|| result.mapped_by_bacterium());
        let l = PyList::new(py, m)?;
        py.import(intern!(py, "array"))?
            .call_method1(intern!(py, "array"), (intern!(py, "f"), l))
            .map(|a| a.into_pyobject(py))?
            .map_err(|e| PyErr::from(e))
    }

    /// Run one or more iterations of the read proportion estimation procedure.
    #[pyo3(signature = (n = 1))]
    pub fn refine(&mut self, n: usize) -> PyResult<()> {
        self.frequencies = None;
        self.proportions = None;
        for _i in 0..n {
            self.result.refine();
        }
        Ok(())
    }
}

/// PyO3 bindings to ``papasmurf``, a library for 16S multiple region analysis.
#[pymodule]
#[pyo3(name = "lib")]
pub fn init<'py>(py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {
    m.add("__package__", "papasmurf")?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS").replace(':', "\n"))?;
    m.add("__build__", pyo3_built!(py, build))?;

    m.add_class::<Database>()?;
    m.add_class::<Builder>()?;
    m.add_class::<Mapper>()?;
    m.add_class::<MapperResult>()?;

    Ok(())
}
