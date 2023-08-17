mod builder;
mod kmers;
mod kmertrie;

use serde::Deserialize;
use serde::Serialize;

use crate::matrix::CsrMatrix;
use crate::matrix::DenseMatrix;
use crate::primer::Primer;

use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

pub use self::builder::Builder;
pub use self::kmers::Kmers;
pub use self::kmertrie::KmerTrie;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct UnindexedRegion {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// The set of forward and backward k-mers in this region.
    pub unique_kmers: Paired<OrderedSet<Rc<str>>>,
    /// A sparse matrix storing the k-mer pair for each database reference.
    pub matrix: CsrMatrix<f32>,
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
        // let mut trie = region.unique_kmers.as_ref().map(|kmers| {
        //     let mut trie = KmerTrie::new(k);
        //     for kmer in kmers.iter() {
        //         trie.insert(kmer);
        //     }
        //     trie
        // });
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
            unique_kmers: region.unique_kmers,
            matrix: region.matrix,
            // trie,
            block,
        }
    }
}

impl From<Region> for UnindexedRegion {
    fn from(region: Region) -> Self {
        Self {
            primer: region.primer,
            unique_pairs: region.unique_pairs,
            unique_kmers: region.unique_kmers,
            matrix: region.matrix,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(from = "UnindexedRegion", into = "UnindexedRegion")]
pub struct Region {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// The set of forward and backward k-mers in this region.
    pub unique_kmers: Paired<OrderedSet<Rc<str>>>,

    // /// A pair of tries storing the unique kmers for the forward and backward region.
    // pub trie: Paired<KmerTrie>,
    /// A pair of blocks storing the unique kmers for the forward and backward region.
    pub block: Paired<DenseMatrix<u8>>,

    /// A sparse matrix storing the k-mer pair for each database reference.
    pub matrix: CsrMatrix<f32>,
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
