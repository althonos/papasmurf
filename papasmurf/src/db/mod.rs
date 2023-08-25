mod builder;
mod kmers;
mod mapper;

use serde::Deserialize;
use serde::Serialize;

use crate::matrix::CsrMatrix;
use crate::matrix::DenseMatrix;
use crate::matrix::MatrixDimensions;
use crate::primer::Primer;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

pub use self::builder::Builder;
pub use self::kmers::Kmers;
pub use self::mapper::Mapper;
pub use self::mapper::MapperResult;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct UnindexedRegion {
    /// The pair of primers defining this region in the database.
    primer: Paired<Primer>,
    /// The set of unique k-mer pairs in this region.
    unique_pairs: OrderedSet<Paired<usize>>,
    /// The set of forward and backward k-mers in this region.
    unique_kmers: Paired<OrderedSet<Rc<str>>>,
    /// A sparse matrix storing the k-mer pair for each database reference.
    matrix: CsrMatrix<f32>,
}

impl From<UnindexedRegion> for Region {
    fn from(region: UnindexedRegion) -> Self {
        let k = region
            .unique_kmers
            .forward
            .iter()
            .next()
            .map(|kmer| kmer.len())
            .unwrap_or_default();
        let block = region.unique_kmers.as_ref().map(|kmers| {
            let mut block = DenseMatrix::new(k, kmers.len());
            for (j, kmer) in kmers.iter().enumerate() {
                for (i, x) in kmer.as_bytes().iter().enumerate() {
                    block[i][j] = *x;
                }
            }
            block
        });
        Self {
            primer: region.primer,
            unique_pairs: region.unique_pairs,
            matrix: region.matrix,
            block: block.map(Kmers::from),
        }
    }
}

impl From<Region> for UnindexedRegion {
    fn from(region: Region) -> Self {
        let unique_kmers = region.block.map(|kmers| {
            let mut unique_kmers = Vec::with_capacity(kmers.columns());
            let mut s = String::with_capacity(kmers.rows());
            for j in 0..kmers.columns() {
                s.clear();
                for i in 0..kmers.rows() {
                    s.push(kmers[i][j] as char);
                }
                unique_kmers.push(Rc::from(s.as_str()));
            }
            unique_kmers.into()
        });
        Self {
            primer: region.primer,
            unique_pairs: region.unique_pairs,
            unique_kmers,
            matrix: region.matrix,
        }
    }
}

/// A single 16S region from a database, defined by a pair of primers.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(from = "UnindexedRegion", into = "UnindexedRegion")]
pub struct Region {
    /// The pair of primers defining this region in the database.
    primer: Paired<Primer>,
    /// The set of unique k-mer pairs in this region.
    unique_pairs: OrderedSet<Paired<usize>>,
    /// A pair of blocks storing the unique kmers for the forward and backward region.
    block: Paired<Kmers>,
    /// A sparse matrix storing the k-mer pair for each database reference.
    matrix: CsrMatrix<f32>,
}

impl Region {
    /// Get a reference to the primer pair used to define this region.
    #[inline]
    pub fn primer(&self) -> &Paired<Primer> {
        &self.primer
    }

    /// Get a reference to the reference matrix for this region, `M`.
    #[inline]
    pub fn matrix(&self) -> &CsrMatrix<f32> {
        &self.matrix
    }
}

/// A database storing forward and backward k-mers for each region.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Database {
    /// The size of the k-mers to extract from the reference sequences.
    k: usize,
    /// The regions this database contains.
    regions: Vec<Region>,
    /// The identifiers of the individual references in the database.
    names: OrderedSet<Rc<str>>,
    /// The number of k-mers extracted from each database reference (R vector).
    amplified: Vec<u8>,
}

impl Database {
    /// Get the length of the k-mers stored in the database.
    #[inline]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get a reference to the names of the database members.
    #[inline]
    pub fn names(&self) -> &[Rc<str>] {
        self.names.as_slice()
    }

    /// Get a reference to the regions stored in the database.
    #[inline]
    pub fn regions(&self) -> &[Region] {
        self.regions.as_slice()
    }
}

impl AsRef<Database> for &Database {
    #[inline]
    fn as_ref(&self) -> &Database {
        *self
    }
}

impl AsRef<[Region]> for &Database {
    #[inline]
    fn as_ref(&self) -> &[Region] {
        self.regions()
    }
}
