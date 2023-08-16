mod builder;
mod kmertrie;
mod kmers;

use serde::Deserialize;
use serde::Serialize;

use crate::matrix::CscMatrix;
use crate::matrix::DenseMatrix;
use crate::primer::Primer;
use crate::utils::Interner;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

pub use self::builder::Builder;
pub use self::kmers::Kmers;
pub use self::kmertrie::KmerTrie;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Region {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// The set of forward and backward k-mers in this region.
    pub unique_kmers: Paired<OrderedSet<Rc<str>>>,
    /// A sparse matrix storing the k-mer pair for each database reference.
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
