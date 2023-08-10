use std::collections::HashMap;
use std::collections::HashSet;

use lightmotif::abc::Dna;
use lightmotif::pli::BestPosition;
use lightmotif::pli::Encode;
use lightmotif::pli::Score;
use lightmotif::pli::Threshold;
use lightmotif::pwm::CountMatrix;
use lightmotif::seq::EncodedSequence;

use crate::matrix::Matrix;
use crate::primer::Primer;
use crate::seq::mismatches;
use crate::seq::reverse_complement;
use crate::seq::DesambiguationIterator;
use crate::utils::Interner;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

type ScoringMatrix = lightmotif::pwm::ScoringMatrix<lightmotif::Dna>;

#[derive(Debug, Clone)]
pub struct BuilderKmer {
    pub id: Rc<str>,
    // pub primer: Paired<Rc<str>>,
    pub kmer: Paired<Rc<str>>,
}

#[derive(Debug, Clone)]
pub struct DatabaseBuilder {
    /// The size of the k-mers to extract from the reference sequences.
    pub k: usize,
    /// The maximum number of mismatches allowed between the primers and a sequence.
    pub primer_mismatches: usize,

    /// The list of primers used to identify the 16S regions.
    pub primers: Vec<Paired<Primer>>,

    /// The list of k-mers extracted from the sequences so far.
    pub kmers: Vec<Vec<BuilderKmer>>,

    /// A string interner, to avoid re-allocating identitical k-mers.
    pub interner: Interner<str>,

    /// The number of sequences that have been added to the database.
    n: usize,
}

impl DatabaseBuilder {
    /// Create a new database builder using the given primers.
    pub fn new(primers: Vec<Paired<Primer>>) -> Self {
        let mut kmers = Vec::with_capacity(primers.len());
        for _ in 0..primers.len() {
            kmers.push(Vec::new());
        }

        DatabaseBuilder {
            primers,
            kmers,
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
                // $pos = $pli.best_position(&$scores).unwrap();
                // $mm = $primer.mismatches(&$seq[$pos..$pos + $primer.len()]);
                let mut indices = $pli.threshold(&$scores, 0.0);
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
            self.kmers[region].push(BuilderKmer {
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
        let k = self.k;

        // Extract the unique names of all the references stored in the builder.
        let names = self
            .kmers
            .iter()
            .flat_map(|kmers| kmers.iter().map(|kmer| &kmer.id))
            .cloned()
            .collect::<OrderedSet<_>>();

        // Prepare the amplified vector
        let mut amplified = vec![0; names.len()];
        let mut regions = Vec::with_capacity(self.primers.len());
        for (primer, kmers) in self.primers.iter().zip(self.kmers.iter()) {
            // Extract unique kmers
            let unique = kmers
                .iter()
                .map(|kmer| &kmer.kmer)
                .cloned()
                .collect::<Paired<HashSet<_>>>()
                .map(OrderedSet::from);

            // Count amplified per region
            for kmer in kmers.iter() {
                amplified[names[&kmer.id]] += 1;
            }

            // Encode reference kmers with indices.
            let entries = kmers
                .iter()
                .map(|kmer| DatabaseEntry {
                    id: names[&kmer.id],
                    // primer: kmer.primer.clone(),
                    kmer_index: Paired::new(
                        unique.forward[&kmer.kmer.forward],
                        unique.backward[&kmer.kmer.backward],
                    ),
                })
                .collect::<Vec<_>>();

            // Extract unique kmer pairs.
            let unique_pairs = entries
                .iter()
                .map(|entry| &entry.kmer_index)
                .cloned()
                .collect();

            // Build PSSMs from the kmer block.
            let profile = unique.as_ref().map(|x| {
                x.iter()
                    .map(AsRef::as_ref)
                    .map(EncodedSequence::<Dna>::encode)
                    .map(Result::unwrap)
                    .collect::<Result<CountMatrix<Dna>, _>>()
                    .unwrap()
                    .to_freq(0.1)
                    .to_scoring(None)
            });

            // Build dense storage for the kmers
            let kmer_block = unique.as_ref().map(|kmers| {
                let mut matrix = Matrix::<u8>::new(kmers.len(), self.k);
                for (i, kmer) in kmers.iter().enumerate() {
                    matrix[i].copy_from_slice(kmer.as_bytes());
                }
                matrix.transpose()
            });

            // Record region
            regions.push(DatabaseRegion {
                primer: primer.clone(),
                // unique_kmers: unique,
                profile,
                entries,
                unique_pairs,
                kmers: kmer_block,
            })
        }

        Database {
            k,
            regions,
            names,
            amplified,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DatabaseEntry {
    pub id: usize,
    pub kmer_index: Paired<usize>,
    // pub primer: Paired<Rc<str>>,
}

#[derive(Debug, Clone)]
pub struct DatabaseRegion {
    /// The pair of primers defining this region in the database.
    pub primer: Paired<Primer>,
    /// The pair of PSSMs for matching this region.
    pub profile: Paired<lightmotif::ScoringMatrix<lightmotif::Dna>>,
    // /// The set of unique k-mers in this region.
    // pub unique_kmers: Paired<OrderedSet<Rc<str>>>,
    /// The set of unique k-mer pairs in this region.
    pub unique_pairs: OrderedSet<Paired<usize>>,
    /// The individual reference entries in this region.
    pub entries: Vec<DatabaseEntry>,
    /// A dense, aligned matrix storing unique forward and backward k-mers.
    pub kmers: Paired<Matrix<u8>>,
}

#[derive(Debug, Clone)]
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
