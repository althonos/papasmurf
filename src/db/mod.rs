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
pub struct Entry {
    pub id: usize,
    pub kmer_index: Paired<usize>,
    // pub primer: Paired<Rc<str>>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Region {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    // /// The pair of PSSMs for matching this region.
    // pub profile: Paired<lightmotif::ScoringMatrix<lightmotif::Dna>>,
    // /// The set of unique k-mers in this region.
    // pub unique_kmers: Paired<OrderedSet<Rc<str>>>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// The individual reference entries in this region.
    pub entries: Vec<Entry>,
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
    pub regions: Vec<Region>,
    /// The identifiers of the individual references in the database.
    pub names: OrderedSet<Rc<str>>,
    /// The number of k-mers extracted from each database reference (R vector).
    pub amplified: Vec<u8>,
}

#[cfg(feature = "serde")]
mod ser {

    use super::*;
    use serde::ser::Serializer;
    use serde::ser::SerializeStruct;
    use serde::Serialize;

    impl Serialize for Region {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            // let mut kmers = self.kmers
            //     .as_ref()
            //     .map(|matrix|  (0..matrix)   )

            let mut state = serializer.serialize_struct("Region", 5)?;
            state.serialize_field("primer", &self.primer)?;
            state.serialize_field("unique_pairs", &self.unique_pairs)?;
            state.serialize_field("entries", &self.entries)?;
            state.serialize_field("kmers", &self.kmers)?;
            state.serialize_field("matrix", &self.matrix)?;
            state.end()
        }
    }
}

// #[cfg(feature = "serde")]
// mod de {

//     use super::*;
//     use serde::de::Deserializer;
//     use serde::de::Error as DeError;
//     use serde::de::Visitor;

//     use serde::Deserialize;

//     struct PrimerVisitor;

//     impl<'de> Visitor<'de> for PrimerVisitor {
//         type Value = Primer;

//         fn expecting(&self, formatter: &mut Formatter) -> FmtResult {
//             write!(formatter, "a string")
//         }

//         fn visit_str<E: DeError>(self, s: &str) -> Result<Self::Value, E> {
//             Ok(Primer::new(s))
//         }
//     }

//     impl<'de> Deserialize<'de> for Primer {
//         fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//         where
//             D: Deserializer<'de>,
//         {
//             deserializer.deserialize_str(PrimerVisitor)
//         }
//     }
// }