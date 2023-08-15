mod builder;
mod kmers;

use serde::Deserialize;
use serde::Serialize;

use crate::matrix::CscMatrix;
use crate::matrix::DenseMatrix;
use crate::primer::Primer;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;
use crate::utils::Interner;

pub use self::builder::Builder;
pub use self::kmers::Kmers;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Entry {
    pub j: usize,
    pub h: usize,
    // pub primer: Paired<Rc<str>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Region {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// A dense, aligned matrix storing unique forward and backward k-mers.
    pub unique_kmers: Paired<Kmers>,
    /// The individual reference entries in this region.
    pub entries: Vec<Entry>,
    /// A sparse matrix storing the k-mer to reference correspondance.
    pub matrix: CscMatrix<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Database {
    /// The size of the k-mers to extract from the reference sequences.
    pub k: usize,
    /// The regions this database contains.
    pub regions: Vec<Region>,
    /// The identifiers of the individual references in the database.
    pub names: OrderedSet<Rc<str>>,
    /// The number of k-mers extracted from each database reference (R vector).
    pub amplified: Vec<u8>,
}
