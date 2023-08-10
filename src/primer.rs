use lightmotif::abc::Alphabet;
use lightmotif::abc::Dna;
use lightmotif::dense::DenseMatrix;
use lightmotif::pwm::ScoringMatrix;

use crate::seq::mismatches;
use crate::seq::reverse_complement;
use crate::seq::DesambiguationIterator;

fn _minscore(
    data: &DenseMatrix<f32, <Dna as Alphabet>::K>,
    startrow: usize,
    mismatches: usize,
) -> f32 {
    if mismatches == 0 {
        (startrow..data.rows())
            .map(|i| {
                *data[i]
                    .iter()
                    .filter(|&&x| x > 0.0)
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap()
            })
            .sum()
    } else if startrow == data.rows() - 1 {
        *data[startrow]
            .iter()
            .min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
    } else {
        let x = *data[startrow]
            .iter()
            .filter(|&&x| x > 0.0)
            .min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
            + _minscore(data, startrow + 1, mismatches);
        let y = *data[startrow]
            .iter()
            .min_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap()
            + _minscore(data, startrow + 1, mismatches - 1);
        x.min(y)
    }
}

#[derive(Debug, Clone)]
pub struct Primer {
    pub template: String,
    pub profile: ScoringMatrix<Dna>,
}

impl Primer {
    pub fn new<S: Into<String>>(template: S) -> Self {
        let t = template.into();
        let mut encoded = DesambiguationIterator::new(&t)
            .map(|s| lightmotif::EncodedSequence::encode(&s).unwrap());
        let pssm = lightmotif::CountMatrix::from_sequences(encoded)
            .unwrap()
            .to_freq(0.1)
            .to_scoring(None);
        Self {
            template: t,
            profile: pssm,
        }
    }

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
        Self::new(reverse_complement(&self.template))
    }
}
