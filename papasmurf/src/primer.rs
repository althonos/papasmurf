use std::fmt::Formatter;
use std::fmt::Result as FmtResult;

use lightmotif::abc::Alphabet;
use lightmotif::abc::Dna;
use lightmotif::dense::DenseMatrix;
use lightmotif::pwm::ScoringMatrix;

use crate::errors::Error;
use crate::seq::mismatches;
use crate::seq::reverse_complement;
use crate::seq::DesambiguationIterator;

// fn _minscore(
//     data: &DenseMatrix<f32, <Dna as Alphabet>::K>,
//     startrow: usize,
//     mismatches: usize,
// ) -> f32 {
//     if mismatches == 0 {
//         (startrow..data.rows())
//             .map(|i| {
//                 *data[i]
//                     .iter()
//                     .filter(|&&x| x > 0.0)
//                     .min_by(|x, y| x.partial_cmp(y).unwrap())
//                     .unwrap()
//             })
//             .sum()
//     } else if startrow == data.rows() - 1 {
//         *data[startrow]
//             .iter()
//             .min_by(|x, y| x.partial_cmp(y).unwrap())
//             .unwrap()
//     } else {
//         let x = *data[startrow]
//             .iter()
//             .filter(|&&x| x > 0.0)
//             .min_by(|x, y| x.partial_cmp(y).unwrap())
//             .unwrap()
//             + _minscore(data, startrow + 1, mismatches);
//         let y = *data[startrow]
//             .iter()
//             .min_by(|x, y| x.partial_cmp(y).unwrap())
//             .unwrap()
//             + _minscore(data, startrow + 1, mismatches - 1);
//         x.min(y)
//     }
// }

#[derive(Debug, Clone)]
pub struct Primer {
    pub template: String,
    pub profile: ScoringMatrix<Dna>,
}

impl Primer {
    /// Create a new primer from the given template DNA.
    ///
    /// # Error
    /// The constructor will fail if given a template string that does not
    /// contain extended DNA symbols.
    ///
    pub fn new<S: Into<String>>(template: S) -> Result<Self, Error> {
        let t = template.into();
        let encoded = DesambiguationIterator::new(&t)?.map(|s| {
            lightmotif::EncodedSequence::encode(&s)
                .expect("DesambiguationIterator only produces valid DNA strings")
        });
        let pssm = lightmotif::CountMatrix::from_sequences(encoded)
            .expect("DesambiguationIterator only produces sequences of the same length")
            .to_freq(0.1)
            .to_scoring(None);
        Ok(Self {
            template: t,
            profile: pssm,
        })
    }

    /// Return the number of symbols in the primer.
    ///
    /// # Example
    /// ``` rust
    /// # extern crate papasmurf;
    /// # use papasmurf::primer::Primer;
    /// let primer = Primer::new("AGGAAGGTGGGGATGACG").unwrap();
    /// assert_eq!(primer.len(), 18);
    /// ```
    pub fn len(&self) -> usize {
        self.template.len()
    }

    /// Compute the number of mismatches between the primer and a sequence.
    pub fn mismatches(&self, seq: &str) -> usize {
        assert_eq!(self.len(), seq.len());
        mismatches(&self.template, seq)
    }

    /// Get the reverse complement of this primer.
    pub fn reverse_complement(&self) -> Primer {
        reverse_complement(&self.template)
            .and_then(Self::new)
            .expect("Primer.reverse_complement always produces a valid Primer")
    }
}

// --- Serde -------------------------------------------------------------------

#[cfg(feature = "serde")]
mod ser {

    use super::*;
    use serde::ser::Serializer;
    use serde::Serialize;

    impl Serialize for Primer {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            serializer.serialize_str(self.template.as_str())
        }
    }
}

#[cfg(feature = "serde")]
mod de {

    use super::*;
    use serde::de::Deserializer;
    use serde::de::Error as DeError;
    use serde::de::Unexpected;
    use serde::de::Visitor;
    use serde::Deserialize;

    struct PrimerVisitor;

    impl<'de> Visitor<'de> for PrimerVisitor {
        type Value = Primer;

        fn expecting(&self, formatter: &mut Formatter) -> FmtResult {
            write!(formatter, "a string")
        }

        fn visit_str<E: DeError>(self, s: &str) -> Result<Self::Value, E> {
            Primer::new(s).map_err(|_| DeError::invalid_value(Unexpected::Str(s), &self))
        }
    }

    impl<'de> Deserialize<'de> for Primer {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            deserializer.deserialize_str(PrimerVisitor)
        }
    }
}
