use std::collections::HashSet;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::RwLock;

use lightmotif::abc::Dna;
use lightmotif::pli::dispatch::Dispatch;
use lightmotif::pli::Encode;
use lightmotif::pli::Pipeline;
use lightmotif::pli::Score;
use lightmotif::pli::Stripe;
use lightmotif::pli::Threshold;

use crate::errors::Error;
use crate::matrix::DokMatrix;
use crate::primer::Primer;
use crate::seq::count_ambiguous;
use crate::seq::reverse_complement;
use crate::seq::DisambiguationIterator;
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
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Sketch {
    /// The name of the reference bacterium this sketch originates from.
    pub name: Rc<str>,
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
    sketches: Vec<RwLock<HashSet<Sketch>>>,
    /// A string interner, to avoid re-allocating identitical k-mers.
    interner: Interner<str>,
    /// The number of references added to the database.
    references: AtomicUsize,
    /// The lightmotif pipeline to run PSSM-related operations.
    pipeline: Pipeline<Dna, Dispatch>,
    /// The length of the largest primer, or `None`.
    largest_primer: Option<usize>,
}

impl Builder {
    /// Create a new database builder using the given primers.
    pub fn new(mut primers: Vec<Paired<Primer>>) -> Self {
        assert!(primers.len() < u8::MAX as usize);

        // Store sketches independently for each region.
        let mut sketches = Vec::with_capacity(primers.len());
        for _ in 0..primers.len() {
            sketches.push(RwLock::new(HashSet::new()));
        }

        // Reverse-complement the backward primer.
        for pair in primers.iter_mut() {
            pair.backward = pair.backward.reverse_complement();
        }

        // Compute length of largest primer
        let largest_primer = primers
            .iter()
            .flat_map(|pair| [pair.forward.profile(), pair.backward.profile()])
            .map(|prof| prof.len())
            .max();

        Builder {
            primers,
            sketches,
            interner: Default::default(),
            k: 100,
            primer_mismatches: 2,
            references: AtomicUsize::new(0),
            pipeline: Pipeline::dispatch(),
            largest_primer,
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
    /// sequence. If the sequence contains ambiguous IUPAC DNA symbols, then
    /// the combination of all possible sequences is generated and used to
    /// extract k-mers, up to 3 ambiguous positions.
    ///
    /// # Error
    /// The method will return an error when `sequence` does not contain valid
    /// IUPAC DNA symbols.
    ///
    pub fn add<I>(&self, name: I, sequence: &str) -> Result<usize, Error>
    where
        I: AsRef<str>,
    {
        let name_ = name.as_ref();
        let mut n = 0;
        if count_ambiguous(&sequence)? <= 3 {
            for dna in DisambiguationIterator::new(&sequence).unwrap() {
                n += self.add_unambiguous(name_, &dna)?;
            }
        }
        Ok(n)
    }

    // Add a single unambiguous sequence to the builder.
    fn add_unambiguous<I>(&self, name: I, sequence: &str) -> Result<usize, Error>
    where
        I: AsRef<str>,
    {
        macro_rules! find_best {
            (
                $primer:expr,
                $striped:ident,
                $scores:ident,
                $seq:expr,
                $pos:ident,
                $mm:ident
            ) => {{
                self.pipeline
                    .score_into(&$primer.profile(), &$striped, &mut $scores);
                $pos = 0;
                $mm = usize::MAX;
                let coordinates = self.pipeline.threshold(&$scores, 0.0);
                for c in coordinates {
                    let i = $scores.offset(c);
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
        let mut scores = lightmotif::scores::StripedScores::<f32, _>::empty();
        // Encode the input sequence
        let mut striped = match self.pipeline.encode(&sequence[..sequence.len() - self.k]) {
            Ok(encoded) => self.pipeline.stripe(encoded),
            Err(_) => return Err(Error::InvalidDna),
        };
        if let Some(&n) = self.largest_primer.as_ref() {
            striped.configure_wrap(n);
        }

        let mut name_rc: Option<Rc<str>> = None;
        let mut amplified = 0;

        for (region, primer) in self.primers.iter().enumerate() {
            // Find the best position for the forward primer
            let mut fwd_pos;
            let mut fwd_mm;
            find_best!(primer.forward, striped, scores, sequence, fwd_pos, fwd_mm);
            // Find the best position for the backward primer
            let mut bwd_pos;
            let mut bwd_mm;
            find_best!(primer.backward, striped, scores, sequence, bwd_pos, bwd_mm);

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
            if self.sketches[region]
                .write()
                .expect("lock was poisoned")
                .insert(Sketch {
                    // primer: Paired::new(fwd_rc, bwd_rc),
                    kmer: Paired::new(fwd_kmer, bwd_kmer),
                    name: name_rc.get_or_insert_with(|| name.as_ref().into()).clone(),
                })
            {
                amplified += 1;
            }
        }

        if amplified > 0 {
            self.references.fetch_add(1, Ordering::Relaxed);
        }

        Ok(amplified)
    }

    /// Build the final database.
    pub fn to_database(&self) -> Database {
        let sketches_ref = self
            .sketches
            .iter()
            .map(|s| s.read().expect("lock was poisoned"))
            .collect::<Vec<_>>();

        // Extract the unique names of all the references stored so far.
        let names = sketches_ref
            .iter()
            .flat_map(|entries| entries.iter().map(|kmer| &kmer.name))
            .cloned()
            .collect::<OrderedSet<_>>();

        // Count how many regions were amplified for each reference.
        let mut amplified = vec![0; names.len()];
        for kmer in sketches_ref.iter().map(|v| v.iter()).flatten() {
            amplified[names[&kmer.name]] += 1;
        }

        // Group kmers for individual regions
        let mut regions = Vec::with_capacity(self.primers.len());
        for (primer, sketches) in self.primers.iter().zip(sketches_ref.iter()) {
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
                let j = names[&sketch.name];
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
