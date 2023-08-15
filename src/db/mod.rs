mod builder;

use serde::Deserialize;
use serde::Serialize;

use crate::matrix::CscMatrix;
use crate::matrix::DenseMatrix;
use crate::primer::Primer;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

pub use self::builder::Builder;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseEntry {
    pub id: usize,
    pub kmer_index: Paired<usize>,
    // pub primer: Paired<Rc<str>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseRegion {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    // /// The pair of PSSMs for matching this region.
    // pub profile: Paired<lightmotif::ScoringMatrix<lightmotif::Dna>>,
    // /// The set of unique k-mers in this region.
    // pub unique_kmers: Paired<OrderedSet<Rc<str>>>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// The individual reference entries in this region.
    pub entries: Vec<DatabaseEntry>,
    /// A dense, aligned matrix storing unique forward and backward k-mers.
    pub kmers: Paired<DenseMatrix<u8>>,
    /// A sparse matrix storing the k-mer to reference correspondance.
    pub matrix: CscMatrix<f32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Database {
    /// The size of the k-mers to extract from the reference sequences.
    pub k: usize,
    /// The regions this database contains.
    pub regions: Vec<DatabaseRegion>,
    /// The identifiers of the individual references in the database.
    pub names: OrderedSet<Rc<str>>,
    /// The number of k-mers extracted from each database reference (R vector).
    pub amplified: Vec<u8>,
}
