use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::RwLock;

use super::db::Database;
use super::errors::Error;
use super::matrix::CooMatrix;

use super::matrix::DokMatrix;
use super::matrix::Dot;
use super::matrix::MatrixDimensions;
use super::matrix::NonZeroElements;
use super::primer::Primer;
use super::utils::Paired;

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

    pub fn with_primer_mismatches(mut self, primer_mismatches: usize) -> Self {
        self.primer_mismatches = primer_mismatches;
        self
    }

    pub fn with_kmer_mismatches(mut self, kmer_mismatches: usize) -> Self {
        self.kmer_mismatches = kmer_mismatches;
        self
    }

    pub fn with_error_probability(mut self, error_probability: f32) -> Self {
        self.error_probability = error_probability;
        self
    }

    pub fn with_partial_hits(mut self, partial_hits: bool) -> Self {
        self.partial_hits = partial_hits;
        self
    }

    fn scan_primer(&self, primer: &Primer, sequence: &str) -> (isize, usize) {
        let min_offset = -(primer.len() as isize) + 1;
        let max_offset = self.primer_region.min(sequence.len() - primer.len()) as isize;

        let mut min_i = isize::MAX;
        let mut min_mm = usize::MAX;

        for i in min_offset..max_offset {
            let mm = if i < 0 {
                let q = &primer.template[(-i) as usize..];
                let t = &sequence[..(primer.len() as isize + i) as usize];
                super::seq::mismatches(q, t) + (-i) as usize
            } else {
                let q = &primer.template;
                let t = &sequence[i as usize..(i + primer.len() as isize) as usize];
                super::seq::mismatches(q, t)
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
            q_matrix = q_matrix + q.to_coo();
        }

        let mut mapped = vec![0; q_matrix.columns()];
        for (_, j, _) in q_matrix.non_zero_elements() {
            mapped[j] += 1;
        }

        MapperResult {
            db: self.db,
            pi: vec![1.0 / q_matrix.columns() as f32; q_matrix.columns()],
            x: vec![0.0; q_matrix.columns()],
            q: q_matrix,
            mapped,
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

#[derive(Debug, Clone)]
pub struct MapperResult<D: AsRef<Database>> {
    db: D,
    pub q: CooMatrix<f32>,
    pub mapped: Vec<usize>,
    pub pi: Vec<f32>,
    pub x: Vec<f32>,
}

impl<D: AsRef<Database>> MapperResult<D> {
    pub fn refine(&mut self) {
        let db = self.db.as_ref();

        // Compute the pi_j vector
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

        // Compute the X_j matrix
        self.x.fill(0.0);
        for j in 0..self.q.columns() {
            if db.amplified[j] > 0 {
                self.x[j] = self.pi[j] / db.amplified[j] as f32;
            }
        }
        let tot = self.x.iter().sum::<f32>();
        if tot > 0.0 {
            for j in 0..self.q.columns() {
                self.x[j] /= tot;
            }
        }
    }
}
