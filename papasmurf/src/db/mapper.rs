use std::assert_eq;
use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::RwLock;

use crate::db::Database;
use crate::errors::Error;
use crate::matrix::CooMatrix;
use crate::matrix::CsrMatrix;
use crate::matrix::DokMatrix;
use crate::matrix::Dot;
use crate::matrix::MatrixDimensions;
use crate::matrix::NonZeroElements;
use crate::matrix::NonZeroElementsMut;
use crate::matrix::VerticalStack;
use crate::primer::Primer;
use crate::utils::Paired;

/// Dead simple counter.
type Counter<K> = HashMap<K, usize>;

/// A helper for mapping 16S reads from a sample to a k-mer database.
#[derive(Debug)]
pub struct Mapper<D: AsRef<Database>> {
    /// The database referenced to by the mapper.
    db: D,
    /// The requested k-mer length for matching.
    kmer_length: usize,
    /// The number of allowed mismatches in the primer region.
    primer_mismatches: usize,
    /// The number of allowed mismatches in the database k-mers region.
    kmer_mismatches: usize,
    /// The constant error probability per nucleotide.
    error_probability: f64,
    /// The length of the region where to look for a primer in the reads.
    primer_region: usize,
    /// Whether or not reads shorter than the database k-mers can be mapped.
    partial_hits: bool,
    /// The number of reads given by the mapper so far.
    added_reads: AtomicUsize,
    /// Minimum read frequency to include in mapping.
    min_read_frequency: f64,
    /// Minimum read count to include in mapping.
    min_read_count: usize,
    /// Maximum number of ambiguous bases per read.
    max_ambiguous: usize,
    /// The read counters recording kmers and count per region.
    reads: Vec<RwLock<Counter<Paired<String>>>>,
}

impl<D: AsRef<Database>> Mapper<D> {
    /// Create a new mapper for the given database and default parameters.
    pub fn new(db: D) -> Self {
        let r = db.as_ref().regions.len();
        Self {
            kmer_length: db.as_ref().k,
            primer_mismatches: 2,
            kmer_mismatches: 2,
            error_probability: 0.005,
            primer_region: 20,
            partial_hits: false,
            min_read_frequency: 1e-4,
            min_read_count: 2,
            max_ambiguous: 0,
            added_reads: AtomicUsize::new(0),
            reads: (0..r).map(|_| RwLock::new(Counter::new())).collect(),
            db,
        }
    }

    /// Get a reference to the database used by this mapper.
    pub fn as_database(&self) -> &Database {
        self.db.as_ref()
    }

    /// Set the number of allowed mismatches in the primer.
    ///
    /// The database references the primer sequences used to define each
    /// region of the 16S gene. In the original SMURF implementation, a
    /// read is discarded when there is not perfect match to any primer of
    /// the database. To allow for reads of worse quality to be processed,
    /// PAPASMURF allows modifying the maximum number of mismatches between
    /// the read and the primers.
    pub fn with_primer_mismatches(mut self, primer_mismatches: usize) -> Self {
        self.primer_mismatches = primer_mismatches;
        self
    }

    /// Set the length of k-mer regions to match.
    pub fn with_kmer_mismatches(mut self, kmer_mismatches: usize) -> Self {
        self.kmer_mismatches = kmer_mismatches;
        self
    }

    /// Set the number of allowed mismatches in the k-mer region.
    pub fn with_kmer_length(mut self, kmer_length: usize) -> Result<Self, Error> {
        if kmer_length > self.db.as_ref().k {
            Err(Error::InvalidKmerLength)
        } else {
            self.kmer_length = kmer_length;
            Ok(self)
        }
    }

    /// Set the error probability used for computing the probability of origin.
    pub fn with_error_probability(mut self, error_probability: f64) -> Self {
        self.error_probability = error_probability;
        self
    }

    /// Toggle whether partial hits are enabled.
    ///
    /// Once the primer sequence removed, a read may be shorter than the
    /// k-mers in the database. If partial hits are disabled, then the read
    /// will be discarded. Otherwise, the partial sequence will be used to
    /// count for mismatches and compute the probability of origin.
    pub fn with_partial_hits(mut self, partial_hits: bool) -> Self {
        self.partial_hits = partial_hits;
        self
    }

    /// Scan a sequence with a primer to find the minimum number of mismatches.
    fn scan_primer(&self, primer: &Primer, sequence: &str) -> (isize, usize) {
        let min_offset = -(primer.len() as isize) + 1;
        let max_offset = self.primer_region.min(sequence.len() - primer.len()) as isize;

        let mut min_i = isize::MAX;
        let mut min_mm = usize::MAX;

        for i in min_offset..max_offset {
            let mm = if i < 0 {
                let q = &primer.template()[(-i) as usize..];
                let t = &sequence[..(primer.len() as isize + i) as usize];
                crate::seq::mismatches(q, t) + (-i) as usize
            } else {
                let q = &primer.template();
                let t = &sequence[i as usize..(i + primer.len() as isize) as usize];
                crate::seq::mismatches(q, t)
            };
            if mm == 0 {
                return (i, mm);
            }
            if mm < min_mm {
                min_i = i;
                min_mm = mm;
            }
        }

        (min_i, min_mm)
    }

    /// Add a read to the mapper.
    pub fn add(&self, read: Paired<&str>) -> Result<bool, Error> {
        let db = self.db.as_ref();
        self.added_reads.fetch_add(1, Ordering::Relaxed);

        // Exclude reads with ambiguous bases
        let ambig = read.as_ref().map(|s| {
            s.matches(|c| c != 'A' && c != 'C' && c != 'G' && c != 'T')
                .count()
        });
        if ambig.forward > self.max_ambiguous || ambig.backward > self.max_ambiguous {
            // println!("discarding: read contains too many ambiguous bases (fwd={} bwd={})", ambig.forward, ambig.backward);
            return Ok(false);
        }

        let mut mapped = false;

        for (r, region) in db.regions.iter().enumerate() {
            let hit_fwd = self.scan_primer(&region.primer.forward, &read.forward);
            let hit_bwd = self.scan_primer(&region.primer.backward, &read.backward);

            let primer_mismatches = Paired::new(hit_fwd.1, hit_bwd.1);
            let pos = Paired::new(hit_fwd.0, hit_bwd.0);

            // Skip if primers mismatch the reads
            if primer_mismatches.forward > self.primer_mismatches
                || primer_mismatches.backward > self.primer_mismatches
            {
                // println!(
                //    "discarding: primer mismatch fwd={} bwd={} (max={})",
                //    primer_mismatches.forward, primer_mismatches.backward, self.primer_mismatches
                // );
                // return Ok(false);
                continue;
            }

            // Create the kmer pair
            let mut kmer = Paired::new(
                &read.forward[(pos.forward + region.primer.forward.len() as isize) as usize..],
                &read.backward[(pos.backward + region.primer.backward.len() as isize) as usize..],
            );

            // Check that the kmer is long enough for the database regions or that
            // partial mapping is enabled in the mapper.
            if kmer.forward.len() > self.kmer_length {
                kmer.forward = &kmer.forward[..self.kmer_length];
            } else if kmer.forward.len() < self.kmer_length && !self.partial_hits {
                // println!("discarding: partial forward kmer");
                continue;
            }
            if kmer.backward.len() > self.kmer_length {
                kmer.backward = &kmer.backward[..self.kmer_length];
            } else if kmer.backward.len() < self.kmer_length && !self.partial_hits {
                // println!("discarding: partial reverse kmer");
                continue;
            }

            // Length is correct, record the read to compute frequencies
            *self.reads[r]
                .write()
                .expect("lock was poisoned")
                .entry(kmer.map(String::from))
                .or_default() += 1;

            // Read is mapped in (at least) one region
            mapped = true;
        }

        // All done for now, the heavy lifting will happen in `Mapper::finish`
        Ok(mapped)
    }

    /// Finish mapping and return the partial results.
    ///
    /// Once all the reads have been processed by the mapper, the final
    /// probability of origin for each read is computed and aggregated for
    /// all regions.
    pub fn finish(self) -> MapperResult<D> {
        // Keep a reference to the database
        let db = self.db.as_ref();

        // Count mapped reads per region
        let mut mapped_reads = vec![0; db.regions.len()];
        let mut assigned_reads = vec![0; db.regions.len()];

        // Count frequency of individual reads per region
        let mut frequencies = vec![Vec::new(); db.regions.len()];

        // Compute E_i,h matrix for each region separately
        let mut expected = vec![HashMap::new(); db.regions.len()];

        // Build region-specific Q matrix while filtering low abundance reads
        for (r, (region, reads)) in db.regions.iter().zip(&self.reads).enumerate() {
            // Retrieve the shared buffer
            let reads = reads.read().expect("lock was poisoned");

            // Compute the minimum number of required reads
            let total_count = reads.values().sum::<usize>();
            let threshold =
                (self.min_read_frequency * total_count as f64).max(self.min_read_count as f64);
            println!("region={} threshold={:?}", r, threshold);

            // Map to region pairs and build the E matrix
            // NOTE(@althonos): embarassingly parallelizable with a lock
            //                  wrapping the result buffers.
            for (kmer, count) in reads.iter().filter(|(_, v)| **v as f64 >= threshold) {
                // Compute mismatches between the read kmer and all the database kmers
                let mismatch = Paired::new(
                    region
                        .block
                        .forward
                        .mismatches(&kmer.forward)
                        .expect("sequence must be valid here"),
                    region
                        .block
                        .backward
                        .mismatches(&kmer.backward)
                        .expect("sequence must be valid here"),
                );
                // Compute mismatch probability
                let mut mapped = false;
                for (h, pair) in region.unique_pairs.iter().enumerate() {
                    let ne = mismatch.forward[pair.forward] as usize
                        + mismatch.backward[pair.backward] as usize;
                    if ne <= self.kmer_mismatches {
                        let l = kmer.forward.len() + kmer.backward.len();
                        let proba = (self.error_probability / 3.0).powi(ne as i32)
                            * (1.0 - self.error_probability).powi((l - ne) as i32);
                        // println!("h={} ne={} e={}", h, ne, e);
                        if proba > 0.0 {
                            let i = frequencies[r].len();
                            expected[r].insert((i, h), proba as f64);
                            mapped = true;
                        }
                    }
                }
                // Record read frequency if it matched any reference k-mer
                if mapped {
                    frequencies[r].push(*count);
                }
            }

            // Count reads assigned and mapped to this region
            assigned_reads[r] = reads.values().sum();
            mapped_reads[r] = frequencies[r].iter().sum();
        }

        // Merge regions Q_i,j into the same matrix
        let mut q = CsrMatrix::<f64>::new(0, db.names.len());
        for (r, (region, expected)) in db.regions.iter().zip(expected.into_iter()).enumerate() {
            let e = DokMatrix::with_data(frequencies[r].len(), region.unique_pairs.len(), expected);
            q.vstack(&e.to_csr().dot(&region.matrix));
        }

        // Merge regions read frequency vectors F_i
        let mapped_count = mapped_reads.iter().sum::<usize>();
        let freq_vec = frequencies
            .iter()
            .flatten()
            .map(|f| *f as f64 / mapped_count as f64)
            .collect::<Vec<f64>>();
        assert_eq!(freq_vec.len(), q.rows());

        // Normalize read frequencies to regions probability
        for (_, j, x) in q.non_zero_elements_mut() {
            *x = *x / db.amplified[j] as f64;
        }

        // Compute initial proportion vector pi = y @ Q
        let mut pi = vec![0.0; q.columns()];
        for (i, j, x) in q.non_zero_elements() {
            pi[j] += *x * freq_vec[i];
        }
        let total = pi.iter().sum::<f64>();
        for x in pi.iter_mut() {
            *x /= total;
        }

        // Return the initial approximation
        // (needs .refine() to get accurate).
        MapperResult {
            pi,
            q: q.into(),
            db: self.db,
            y: freq_vec,
            assigned_reads,
            mapped_reads,
        }
    }
}

impl<D: AsRef<Database>> AsRef<D> for Mapper<D> {
    fn as_ref(&self) -> &D {
        &self.db
    }
}

impl<D: AsRef<Database>> AsRef<Database> for Mapper<D> {
    fn as_ref(&self) -> &Database {
        self.db.as_ref()
    }
}

/// The results of a database mapping.
///
/// Once all reads have been mapped against the database k-mers, the final
/// `Q` probability matrix is computed by aggregating all regions.
#[derive(Debug, Clone)]
pub struct MapperResult<D: AsRef<Database>> {
    db: D,
    q: CooMatrix<f64>, // dim[i, j]
    y: Vec<f64>,       // dim[i]
    pi: Vec<f64>,      // dim[j]
    assigned_reads: Vec<usize>,
    mapped_reads: Vec<usize>,
}

impl<D: AsRef<Database>> MapperResult<D> {
    /// Get a reference to the database used by the mapper.
    #[inline]
    pub fn as_database(&self) -> &Database {
        self.db.as_ref()
    }

    /// Get a reference to the number of assigned reads per region.
    #[inline]
    pub fn assigned_by_region(&self) -> &[usize] {
        &self.assigned_reads
    }

    /// Get a reference to the number of mapped reads per region.
    #[inline]
    pub fn mapped_by_region(&self) -> &[usize] {
        &self.mapped_reads
    }

    /// Get a reference to the read probability matrix, `Q`.
    #[inline]
    pub fn probabilities(&self) -> &CooMatrix<f64> {
        &self.q
    }

    /// Get a reference to the read proportion vector, `π`.
    #[inline]
    pub fn proportions(&self) -> &[f64] {
        &self.pi
    }

    /// Compute the bacterium frequency vector, `X`.
    pub fn frequencies(&self) -> Vec<f64> {
        let db = self.db.as_ref();

        // Compute frequency normalizing proportions by number of
        // regions and hard-thresholding frequencies
        let mut freq = self
            .pi
            .iter()
            .zip(&db.amplified)
            .map(|(&pi, &r)| if pi > 1e-10 { pi / r as f64 } else { 0.0 })
            .collect::<Vec<_>>();

        // Renormalize
        let tot = freq.iter().sum::<f64>();
        for x in freq.iter_mut() {
            *x /= tot;
        }

        freq
    }

    /// Compute the number of reads mapped to each bacterium.
    pub fn mapped_by_bacterium(&self) -> Vec<usize> {
        let mut mapped = vec![0; self.q.columns()];
        for (_, j, _) in self.q.non_zero_elements() {
            mapped[j] += 1;
        }
        mapped
    }

    /// Run one iteration of the read proportion estimation procedure.
    ///
    /// Returns the L1 error.
    pub fn refine(&mut self) -> f64 {
        let _db = self.db.as_ref();

        // Estimate theta
        let mut theta = vec![0.0; self.q.rows()];
        for (i, j, x) in self.q.non_zero_elements() {
            theta[i] += *x * self.pi[j];
        }

        // Reweight the counts (reuse theta buffer)
        for i in 0..self.q.rows() {
            theta[i] = self.y[i] / (theta[i] + f64::EPSILON);
        }

        // Update the refinement factor
        let mut factor = vec![0.0; self.q.columns()];
        for (i, j, x) in self.q.non_zero_elements() {
            factor[j] += *x * theta[i];
        }

        // Compute L1 error
        let l1 = factor
            .iter()
            .zip(&self.pi)
            .map(|(fact, prop)| (1.0f64 - fact).abs() * prop)
            .sum::<f64>();

        // Update unknown read proportion vector
        for j in 0..self.q.columns() {
            self.pi[j] *= factor[j];
        }

        // Remove bacteria of low frequency
        // for x in self.pi.iter_mut() {
        //     if *x < 1e-10 {
        //         *x = 0.0;
        //     }
        // }

        // Return L1 error
        l1
    }
}

impl<D: AsRef<Database>> AsRef<D> for MapperResult<D> {
    fn as_ref(&self) -> &D {
        &self.db
    }
}

impl<D: AsRef<Database>> AsRef<Database> for MapperResult<D> {
    fn as_ref(&self) -> &Database {
        self.db.as_ref()
    }
}
