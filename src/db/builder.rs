use std::collections::HashSet;

use lightmotif::abc::Dna;
use lightmotif::pli::Encode;
use lightmotif::pli::Score;
use lightmotif::pli::Threshold;
use lightmotif::pwm::CountMatrix;
use lightmotif::seq::EncodedSequence;

use crate::matrix::DenseMatrix;
use crate::matrix::DokMatrix;
use crate::primer::Primer;
use crate::seq::reverse_complement;
use crate::utils::Interner;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

use super::Database;
use super::Entry;
use super::Region;

#[derive(Debug, Clone)]
struct Sketch {
    pub id: Rc<str>,
    // pub primer: Paired<Rc<str>>,
    pub kmer: Paired<Rc<str>>,
}

#[derive(Debug, Clone)]
pub struct Builder {
    /// The size of the k-mers to extract from the reference sequences.
    k: usize,
    /// The maximum number of mismatches allowed between the primers and a sequence.
    primer_mismatches: usize,
    /// The list of primers used to identify the 16S regions.
    primers: Vec<Paired<Primer>>,
    /// The list of sketches extracted so far, grouped by region.
    sketches: Vec<Vec<Sketch>>,
    /// A string interner, to avoid re-allocating identitical k-mers.
    interner: Interner<str>,
    /// The number of sequences that have been added to the database.
    n: usize,
}

impl Builder {
    /// Create a new database builder using the given primers.
    pub fn new(primers: Vec<Paired<Primer>>) -> Self {
        let mut sketches = Vec::with_capacity(primers.len());
        for _ in 0..primers.len() {
            sketches.push(Vec::new());
        }

        Builder {
            primers,
            sketches,
            interner: Default::default(),
            k: 100,
            primer_mismatches: 2,
            n: 0,
        }
    }

    /// Add a new reference sequence to the database.
    ///
    /// Returns the number of region k-mers successfully extracted from the
    /// sequence.
    ///
    pub fn add<I>(&mut self, id: I, sequence: &str) -> usize
    where
        I: AsRef<str>,
    {
        macro_rules! find_best {
            (
                $pli:ident,
                $primer:expr,
                $striped:ident,
                $scores:ident,
                $seq:expr,
                $pos:ident,
                $mm:ident
            ) => {{
                $pli.score_into(&$striped, &$primer.profile, &mut $scores);
                $pos = 0;
                $mm = usize::MAX;
                // $pos = $pli.argmax(&$scores).unwrap();
                // $mm = $primer.mismatches(&$seq[$pos..$pos + $primer.len()]);
                let indices = $pli.threshold(&$scores, 0.0);
                for i in indices {
                    let mm_i = $primer.mismatches(&$seq[i..i + $primer.len()]);
                    if mm_i < $mm || ((mm_i == $mm) && (i < $pos)) {
                        $mm = mm_i;
                        $pos = i;
                    }
                }
                if $mm > self.primer_mismatches {
                    continue;
                }
            }};
        }

        let pli = lightmotif::Pipeline::avx2().unwrap();
        let mut scores = lightmotif::pli::StripedScores::empty();

        let mut striped = match pli.encode(&sequence[..sequence.len() - self.k]) {
            Ok(encoded) => lightmotif::seq::EncodedSequence::from(encoded).to_striped(),
            Err(_) => return 0,
        };
        if let Some(profile) = self
            .primers
            .iter()
            .flat_map(|pair| [&pair.forward.profile, &pair.backward.profile])
            .max_by_key(|prof| prof.len())
        {
            striped.configure(&profile);
        }

        let mut id_rc: Option<Rc<str>> = None;
        let mut amplified = 0;

        for (region, primer) in self.primers.iter().enumerate() {
            // Find the best position for the forward primer
            let mut fwd_pos;
            let mut fwd_mm;
            find_best!(
                pli,
                primer.forward,
                striped,
                scores,
                sequence,
                fwd_pos,
                fwd_mm
            );
            // Find the best position for the backward primer
            let mut bwd_pos;
            let mut bwd_mm;
            find_best!(
                pli,
                primer.backward,
                striped,
                scores,
                sequence,
                bwd_pos,
                bwd_mm
            );

            // // Extract and intern the sequence of the forward primer
            // let fwd_seq = &sequence[fwd_pos..fwd_pos + primer.forward.len()];
            // let fwd_rc = self.interner.intern(fwd_seq);
            // // Extract and intern the sequence of the backward primer
            // let bwd_seq = &sequence[bwd_pos..sequence.len().min(bwd_pos + primer.backward.len())];
            // let bwd_rc = self.interner.intern(bwd_seq);

            if fwd_pos >= bwd_pos {
                continue;
            }
            if fwd_pos + primer.forward.len() + self.k > sequence.len() {
                continue;
            }
            if bwd_pos < self.k {
                continue;
            }

            // Extract and intern the k-mer
            let fwd_kmer = self.interner.intern(
                &sequence[fwd_pos + primer.forward.len()..fwd_pos + primer.forward.len() + self.k],
            );
            let bwd_kmer = self
                .interner
                .intern(&reverse_complement(&sequence[bwd_pos - self.k..bwd_pos]));

            // Add the amplified k-mer to the current region.
            amplified += 1;
            self.sketches[region].push(Sketch {
                // primer: Paired::new(fwd_rc, bwd_rc),
                kmer: Paired::new(fwd_kmer, bwd_kmer),
                id: id_rc.get_or_insert_with(|| id.as_ref().into()).clone(),
            });
        }

        if amplified > 0 {
            self.n += 1;
        }

        amplified
    }

    /// Build the final database.
    pub fn to_database(&self) -> Database {
        // Extract the unique names of all the references stored so far.
        let names = self
            .sketches
            .iter()
            .flat_map(|entries| entries.iter().map(|kmer| &kmer.id))
            .cloned()
            .collect::<OrderedSet<_>>();

        // Count how many regions were amplified for each reference.
        let mut amplified = vec![0; names.len()];
        for kmer in self.sketches.iter().map(|v| v.iter()).flatten() {
            amplified[names[&kmer.id]] += 1;
        }

        // Group kmers for individual regions
        let mut regions = Vec::with_capacity(self.primers.len());
        for (primer, sketches) in self.primers.iter().zip(self.sketches.iter()) {
            // Extract unique kmers
            let unique = sketches
                .iter()
                .map(|sketch| &sketch.kmer)
                .cloned()
                .collect::<Paired<HashSet<_>>>()
                .map(OrderedSet::from);

            // Encode reference kmers with indices.
            let entries = sketches
                .iter()
                .map(|sketch| Entry {
                    id: names[&sketch.id],
                    // primer: kmer.primer.clone(),
                    kmer_index: Paired::new(
                        unique.forward[&sketch.kmer.forward],
                        unique.backward[&sketch.kmer.backward],
                    ),
                })
                .collect::<Vec<_>>();

            // Extract unique kmer pairs.
            let unique_pairs = entries
                .iter()
                .map(|sketch| &sketch.kmer_index)
                .cloned()
                .collect::<OrderedSet<_>>();

            // Build PSSMs from the kmer block.
            let _profile = unique.as_ref().map(|x| {
                x.iter()
                    .map(AsRef::as_ref)
                    .map(EncodedSequence::<Dna>::encode)
                    .map(Result::unwrap)
                    .collect::<Result<CountMatrix<Dna>, _>>()
                    .expect("sequences stored in the builder should be unambiguous DNA")
                    .to_freq(0.1)
                    .to_scoring(None)
            });

            // Build dense storage for the kmers
            let kmer_block = unique.as_ref().map(|kmers| {
                let mut matrix = DenseMatrix::<u8>::new(self.k, kmers.len());
                for (i, kmer) in kmers.iter().enumerate() {
                    for (j, x) in kmer.as_bytes().iter().enumerate() {
                        matrix[j][i] = *x;
                    }
                }
                matrix
            });

            // Build M_hj matrix
            let mut matrix = DokMatrix::new(unique_pairs.len(), names.len());
            for sketch in entries.iter() {
                let h = unique_pairs[&sketch.kmer_index];
                let j = sketch.id;
                matrix.insert(h, j, 1.0 / amplified[j] as f32);
            }

            // Record region
            regions.push(Region {
                primer: primer.clone(),
                // profile,
                entries,
                unique_pairs,
                matrix: matrix.to_csc(),
                kmers: kmer_block,
            })
        }

        Database {
            k: self.k,
            regions,
            names,
            amplified,
        }
    }
}
