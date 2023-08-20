use std::collections::HashSet;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use lightmotif::pli::Encode;
use lightmotif::pli::Score;
use lightmotif::pli::Threshold;

use crate::errors::Error;
use crate::matrix::DokMatrix;
use crate::primer::Primer;
use crate::seq::reverse_complement;
use crate::utils::Interner;
use crate::utils::OrderedSet;
use crate::utils::Paired;
use crate::utils::Rc;

use super::Database;
use super::UnindexedRegion;

/// The k-mer sketched from a single sequence.
///
/// The k-mer are stored as raw strings at this stage, and will be later
/// deduplicated when the builder is transformed into a fully-indexed database.
#[derive(Debug, Clone)]
struct Sketch {
    /// The idenfitier of the reference this sketch originates from.
    pub id: Rc<str>,
    /// The forward and backward k-mers extracted from the reference sequence.
    pub kmer: Paired<Rc<str>>,
}

/// A builder for incremental construction of a new database.
#[derive(Debug)]
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
    /// The number of references added to the database.
    references: AtomicUsize,
}

impl Builder {
    /// Create a new database builder using the given primers.
    pub fn new(mut primers: Vec<Paired<Primer>>) -> Self {
        assert!(primers.len() < u8::MAX as usize);

        // Store sketches independently for each region.
        let mut sketches = Vec::with_capacity(primers.len());
        for _ in 0..primers.len() {
            sketches.push(Vec::new());
        }

        // Reverse-complement the backward primer.
        for pair in primers.iter_mut() {
            pair.backward = pair.backward.reverse_complement();
        }

        Builder {
            primers,
            sketches,
            interner: Default::default(),
            k: 100,
            primer_mismatches: 2,
            references: AtomicUsize::new(0),
        }
    }

    /// Set the number of allowed primer mismatches.
    pub fn with_primer_mismatches(mut self, primer_mismatches: usize) -> Self {
        self.primer_mismatches = primer_mismatches;
        self
    }

    /// Add a new reference sequence to the database.
    ///
    /// Returns the number of region k-mers successfully extracted from the
    /// sequence.
    ///
    /// # Error
    /// The method will return an error when `sequence` does not contain valid
    /// DNA symbols (*A*, *T*, *G*, *C* or *N*).
    ///
    pub fn add<I>(&mut self, id: I, sequence: &str) -> Result<usize, Error>
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

        // Create a lightmotif pipeline to search for the primer.
        let pli = lightmotif::Pipeline::avx2().unwrap();
        let mut scores = lightmotif::pli::StripedScores::<lightmotif::num::U32>::empty();

        // Encode the input sequence
        let mut striped = match pli.encode(&sequence[..sequence.len() - self.k]) {
            Ok(encoded) => lightmotif::seq::EncodedSequence::from(encoded).to_striped(),
            Err(_) => return Err(Error::InvalidDna),
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

            // Extract and intern the sequence of the forward primer
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
                .intern(&reverse_complement(&sequence[bwd_pos - self.k..bwd_pos])?);

            // Add the amplified k-mer to the current region.
            amplified += 1;
            self.sketches[region].push(Sketch {
                // primer: Paired::new(fwd_rc, bwd_rc),
                kmer: Paired::new(fwd_kmer, bwd_kmer),
                id: id_rc.get_or_insert_with(|| id.as_ref().into()).clone(),
            });
        }

        if amplified > 0 {
            self.references.fetch_add(1, Ordering::Relaxed);
        }

        Ok(amplified)
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

            // Extract unique kmer pairs
            let unique_pairs: OrderedSet<Paired<usize>> = sketches
                .iter()
                .map(|sketch| &sketch.kmer)
                .map(|kmer| {
                    Paired::new(
                        unique.forward[&kmer.forward],
                        unique.backward[&kmer.backward],
                    )
                })
                .collect::<HashSet<Paired<_>>>()
                .into();

            // Build M_hj matrix
            let mut matrix = DokMatrix::new(unique_pairs.len(), names.len());
            for sketch in sketches.iter() {
                let j = names[&sketch.id];
                let h = unique_pairs[&Paired::new(
                    unique.forward[&sketch.kmer.forward],
                    unique.backward[&sketch.kmer.backward],
                )];
                if amplified[j] > 0 {
                    matrix.insert(h, j, 1.0 / amplified[j] as f32);
                }
            }

            let mut region_primer = primer.clone();
            region_primer.backward = region_primer.backward.reverse_complement();

            // Record region
            regions.push(
                UnindexedRegion {
                    primer: region_primer,
                    unique_pairs,
                    matrix: matrix.to_csr(),
                    unique_kmers: unique,
                }
                .into(),
            )
        }

        Database {
            k: self.k,
            regions,
            names,
            amplified,
        }
    }
}

impl From<Builder> for Database {
    fn from(builder: Builder) -> Self {
        builder.to_database()
    }
}

impl<'b> From<&'b Builder> for Database {
    fn from(builder: &'b Builder) -> Self {
        builder.to_database()
    }
}
