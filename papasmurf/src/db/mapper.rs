use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::RwLock;

use crate::db::Database;
use crate::errors::Error;
use crate::matrix::CooMatrix;
use crate::matrix::DokMatrix;
use crate::matrix::Dot;
use crate::matrix::MatrixDimensions;
use crate::matrix::NonZeroElements;
use crate::primer::Primer;
use crate::utils::Paired;

/// A helper for mapping 16S reads from a sample to a k-mer database.
#[derive(Debug)]
pub struct Mapper<D: AsRef<Database>> {
    pub db: D,
    pub expected: Vec<RwLock<HashMap<(usize, usize), f32>>>,
    primer_mismatches: usize,
    kmer_mismatches: usize,
    error_probability: f32,
    primer_region: usize,
    partial_hits: bool,
    pub reads: AtomicUsize,
}

impl<D: AsRef<Database>> Mapper<D> {
    /// Create a new mapper for the given database.
    pub fn new(db: D) -> Self {
        let expected = db
            .as_ref()
            .regions
            .iter()
            .map(|_| RwLock::from(HashMap::new()))
            .collect();
        Self {
            expected,
            db,
            primer_mismatches: 2,
            kmer_mismatches: 2,
            error_probability: 0.005,
            primer_region: 20,
            partial_hits: false,
            reads: AtomicUsize::new(0),
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

    /// Set the number of allowed mismatches in the k-mer region.
    pub fn with_kmer_mismatches(mut self, kmer_mismatches: usize) -> Self {
        self.kmer_mismatches = kmer_mismatches;
        self
    }

    /// Set the error probability used for computing the probability of origin.
    pub fn with_error_probability(mut self, error_probability: f32) -> Self {
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

        // Add a new row to the E_i,h matrices
        let i = self.reads.fetch_add(1, Ordering::Relaxed);

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

        // Skip if primers mismatch the reads
        if primer_mismatches.forward > self.primer_mismatches
            || primer_mismatches.backward > self.primer_mismatches
        {
            return Ok(false);
        }

        // Create the kmer pair
        let mut kmer = Paired::new(
            &read.forward[(pos.forward + region.primer.forward.len() as isize) as usize..],
            &read.backward[(pos.backward + region.primer.backward.len() as isize) as usize..],
        );

        // Check that the kmer is long enough for the database regions or that
        // partial mapping is enabled in the mapper.
        if kmer.forward.len() > db.k {
            kmer.forward = &kmer.forward[..db.k];
        } else if kmer.forward.len() < db.k && !self.partial_hits {
            return Ok(false);
        }
        if kmer.backward.len() > db.k {
            kmer.backward = &kmer.backward[..db.k];
        } else if kmer.backward.len() < db.k && !self.partial_hits {
            return Ok(false);
        }

        // Compute mismatches between the read kmer and all the database kmers
        let mismatch = Paired::new(
            region.block.forward.mismatches(kmer.forward)?,
            region.block.backward.mismatches(kmer.backward)?,
        );

        // Record the read if it matches any database kmer
        let mut mapped = false;
        for (h, pair) in region.unique_pairs.iter().enumerate() {
            if mismatch.forward[pair.forward] as usize <= self.kmer_mismatches
                && mismatch.backward[pair.backward] as usize <= self.kmer_mismatches
            {
                let l = kmer.forward.len() + kmer.backward.len();
                let ne =
                    (mismatch.forward[pair.forward] + mismatch.backward[pair.backward]) as usize;
                let e = (self.error_probability / 3.0).powf(ne as f32)
                    * (1.0 - self.error_probability).powf((l - ne) as f32);
                if e > 0.0 {
                    self.expected[r]
                        .write()
                        .expect("lock was poisoned")
                        .insert((i, h), e);
                    mapped = true;
                }
            }
        }

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
        let mut q_matrix = CooMatrix::<f32>::new(reads, db.names.len());
        for (region, expected) in db.regions.iter().zip(self.expected) {
            let e = DokMatrix::with_data(
                reads,
                region.unique_pairs.len(),
                expected.into_inner().unwrap(),
            );
            let q = e.to_csr().dot(&region.matrix);
            q_matrix = q_matrix + q.into_coo();
        }

        MapperResult {
            db: self.db,
            pi: vec![1.0 / q_matrix.columns() as f32; q_matrix.columns()],
            q: q_matrix,
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
    q: CooMatrix<f32>,
    pi: Vec<f32>,
}

impl<D: AsRef<Database>> MapperResult<D> {
    /// Get a reference to the database used by the mapper.
    #[inline]
    pub fn as_database(&self) -> &Database {
        self.db.as_ref()
    }

    /// Get a reference to read probability matrix, `Q`.
    #[inline]
    pub fn probabilities(&self) -> &CooMatrix<f32> {
        &self.q
    }

    /// Get a reference to the read proportion vector, `π`.
    #[inline]
    pub fn proportions(&self) -> &[f32] {
        &self.pi
    }

    /// Compute the bacterium frequency vector, `X`.
    pub fn frequencies(&self) -> Vec<f32> {
        let db = self.db.as_ref();
        let mut x = Vec::with_capacity(self.q.columns());
        for j in 0..self.q.columns() {
            if db.amplified[j] > 0 {
                x.push(self.pi[j] / db.amplified[j] as f32);
            } else {
                x.push(0.0);
            }
        }
        let tot = x.iter().sum::<f32>();
        if tot > 0.0 {
            for j in 0..self.q.columns() {
                x[j] /= tot;
            }
        }
        x
    }

    /// Compute the number of reads mapped to each reference bacterium.
    pub fn mapped(&self) -> Vec<usize> {
        let mut mapped = vec![0; self.q.columns()];
        for (_, j, _) in self.q.non_zero_elements() {
            mapped[j] += 1;
        }
        mapped
    }

    /// Run one iteration of the read proportion estimation procedure.
    pub fn refine(&mut self) {
        let _db = self.db.as_ref();
        let mut up = vec![0.0; self.q.columns()];
        let mut dens = vec![0.0; self.q.rows()];
        dens.fill(0.0);
        for (i, j, x) in self.q.non_zero_elements() {
            dens[i] += x * self.pi[j];
        }
        up.fill(0.0);
        for (i, j, x) in self.q.non_zero_elements() {
            if dens[i] > 0.0 {
                up[j] += *x / dens[i]
            }
        }
        for j in 0..self.q.columns() {
            self.pi[j] *= up[j] / self.q.rows() as f32;
        }
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
