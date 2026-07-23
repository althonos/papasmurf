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
use crate::matrix::VerticalStack;
use crate::primer::Primer;
use crate::utils::Paired;

/// A helper for mapping 16S reads from a sample to a k-mer database.
#[derive(Debug)]
pub struct Mapper<D: AsRef<Database>> {
    /// The database referenced to by the mapper.
    db: D,
    /// The requested k-mer length for matching.
    kmer_length: usize,
    /// The `E` read matching probability matrix for each region.
    expected: Vec<RwLock<HashMap<(usize, usize), f64>>>,
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
    reads: AtomicUsize,
    /// The number of reads assigned to each region after primer scanning.
    assigned_reads: Vec<AtomicUsize>,
    /// The number of reads mapped to each region after quality filtering.
    mapped_reads: Vec<AtomicUsize>,
}

impl<D: AsRef<Database>> Mapper<D> {
    /// Create a new mapper for the given database.
    pub fn new(db: D) -> Self {
        let r = db.as_ref().regions.len();
        let expected = db
            .as_ref()
            .regions
            .iter()
            .map(|_| RwLock::from(HashMap::new()))
            .collect();
        Self {
            kmer_length: db.as_ref().k,
            primer_mismatches: 2,
            kmer_mismatches: 2,
            error_probability: 0.005,
            primer_region: 20,
            partial_hits: false,
            reads: AtomicUsize::new(0),
            assigned_reads: (0..r).map(|_| AtomicUsize::new(0)).collect(),
            mapped_reads: (0..r).map(|_| AtomicUsize::new(0)).collect(),
            expected,
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

        // Find the best matching primer and primer position
        let (r, region, pos, primer_mismatches) = db
            .regions
            .iter()
            .enumerate()
            .map(|(r, region)| {
                (
                    r,
                    region,
                    self.scan_primer(&region.primer.forward, &read.forward),
                    self.scan_primer(&region.primer.backward, &read.backward),
                )
            })
            .min_by(|x, y| (x.2 .1 + x.3 .1).partial_cmp(&(y.2 .1 + y.3 .1)).unwrap())
            .map(|x| {
                (
                    x.0,
                    x.1,
                    Paired::new(x.2 .0, x.3 .0),
                    Paired::new(x.2 .1, x.3 .1),
                )
            })
            .unwrap();

        // Keep count of the number of reads assigned to each region
        self.assigned_reads[r].fetch_add(1, Ordering::Relaxed);

        // Skip if primers mismatch the reads
        if primer_mismatches.forward > self.primer_mismatches
            || primer_mismatches.backward > self.primer_mismatches
        {
            //println!(
            //    "discarding: primer mismatch fwd={} bwd={} (max={})",
            //    primer_mismatches.forward, primer_mismatches.backward, self.primer_mismatches
            //);
            return Ok(false);
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
            //println!("discarding: partial forward kmer");
            return Ok(false);
        }
        if kmer.backward.len() > self.kmer_length {
            kmer.backward = &kmer.backward[..self.kmer_length];
        } else if kmer.backward.len() < self.kmer_length && !self.partial_hits {
            //println!("discarding: partial reverse kmer");
            return Ok(false);
        }

        // Compute mismatches between the read kmer and all the database kmers
        let mismatch = Paired::new(
            region.block.forward.mismatches(kmer.forward)?,
            region.block.backward.mismatches(kmer.backward)?,
        );

        // Record the read if it matches any database kmer
        let mut mapped = false;
        let mut index = None;

        for (h, pair) in region.unique_pairs.iter().enumerate() {
            let ne =
                mismatch.forward[pair.forward] as usize + mismatch.backward[pair.backward] as usize;
            if ne <= self.kmer_mismatches {
                let l = kmer.forward.len() + kmer.backward.len();
                let e = (self.error_probability / 3.0).powi(ne as i32)
                    * (1.0 - self.error_probability).powi((l - ne) as i32);
                //println!("h={} ne={} e={}", h, ne, e);
                if e > 0.0 {
                    let i = match index {
                        Some(i) => i,
                        None => {
                            let i = self.reads.fetch_add(1, Ordering::Relaxed);
                            index = Some(i);
                            i
                        }
                    };

                    self.expected[r]
                        .write()
                        .expect("lock was poisoned")
                        .insert((i, h), e as f64);
                    mapped = true;
                }
            }
        }

        // Keep count of the number of reads mapped to each region
        if mapped {
            self.mapped_reads[r].fetch_add(1, Ordering::Relaxed);
        }

        // Return whether the read was mapped
        Ok(mapped)
    }

    /// Finish mapping and return the partial results.
    ///
    /// Once all the reads have been processed by the mapper, the final
    /// probability of origin for each read is computed and aggregated for
    /// all regions.
    pub fn finish(self) -> MapperResult<D> {
        let db = self.db.as_ref();
        let reads = self.reads.load(Ordering::Relaxed);

        // Compute the Q_i,j matrix
        let mut q_matrix = CsrMatrix::<f64>::new(0, db.names.len());
        for (region, expected) in db.regions.iter().zip(self.expected) {
            let e = DokMatrix::with_data(
                reads,
                region.unique_pairs.len(),
                expected.into_inner().unwrap(),
            );
            let q = e.to_csr().dot(&region.matrix);
            q_matrix.vstack(q);
        }

        // Compute initial read proportion vector pi as the sum of mapped reads
        let n = self
            .mapped_reads
            .iter()
            .map(|x| x.load(Ordering::Relaxed))
            .sum::<usize>();
        let pi = vec![1.0 / n as f64; q_matrix.columns()];

        // Recover counts by region
        let assigned_reads = self
            .assigned_reads
            .into_iter()
            .map(|count| count.load(Ordering::Relaxed))
            .collect();
        let mapped_reads = self
            .mapped_reads
            .into_iter()
            .map(|count| count.load(Ordering::Relaxed))
            .collect();

        MapperResult {
            db: self.db,
            q: q_matrix.to_coo(),
            pi,
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
    q: CooMatrix<f64>,
    pi: Vec<f64>,
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

        // Compute frequency by normalizing by number of regions
        let mut freq = self
            .pi
            .iter()
            .zip(&db.amplified)
            .map(|(pi, r)| if *r > 0 { pi / *r as f64 } else { 0.0 })
            .map(|f| if f < 1e-10 { 0.0 } else { f })
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
        let n = self.mapped_reads.iter().sum::<usize>();

        // Compute d_i = sum_j Q_i,j p_j for each i
        let mut dens = vec![0.0; self.q.rows()];
        for (i, j, x) in self.q.non_zero_elements() {
            dens[i] += x * self.pi[j];
        }

        // Compute u_j = sum_i Q_i,j / d_i for each j
        let mut up = vec![0.0; self.q.columns()];
        for (i, j, x) in self.q.non_zero_elements() {
            up[j] += *x / (dens[i] + f64::EPSILON);
        }

        // Compute update factor
        let mut factor = vec![0.0; self.q.columns()];
        for j in 0..self.q.columns() {
            factor[j] = up[j] / n as f64;
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
        for x in self.pi.iter_mut() {
            if *x < 1e-10 {
                *x = 0.0;
            }
        }

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
