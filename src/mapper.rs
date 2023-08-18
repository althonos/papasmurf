use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::sync::RwLock;

use super::db::Database;
use super::matrix::CooMatrix;
use super::matrix::DenseMatrix;
use super::matrix::DokMatrix;
use super::matrix::Dot;
use super::matrix::MatrixDimensions;
use super::matrix::NonZeroElements;
use super::primer::Primer;
use super::utils::Paired;

fn simd_mismatches(query: &[u8], db: &DenseMatrix<u8>, out: &mut [u8]) {
    use std::arch::x86_64::*;
    unsafe {
        let _k = db.rows();

        let mut c = 0;

        let ones = _mm256_set1_epi8(1);

        while c + std::mem::size_of::<__m256i>() * 4 < db.columns() {
            let mut m1 = _mm256_setzero_si256();
            let mut m2 = _mm256_setzero_si256();
            let mut m3 = _mm256_setzero_si256();
            let mut m4 = _mm256_setzero_si256();

            for i in 0..query.len() {
                if query[i] != b'N' {
                    let q = _mm256_set1_epi8(query[i] as i8);
                    let r1 = _mm256_load_si256(db[i].as_ptr().add(c) as *const _);
                    let r2 = _mm256_load_si256(db[i].as_ptr().add(c + 32) as *const _);
                    let r3 = _mm256_load_si256(db[i].as_ptr().add(c + 64) as *const _);
                    let r4 = _mm256_load_si256(db[i].as_ptr().add(c + 96) as *const _);
                    m1 = _mm256_add_epi8(m1, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r1), ones));
                    m2 = _mm256_add_epi8(m2, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r2), ones));
                    m3 = _mm256_add_epi8(m3, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r3), ones));
                    m4 = _mm256_add_epi8(m4, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r4), ones));
                }
            }

            _mm256_storeu_si256(out.as_mut_ptr().add(c) as *mut _, m1);
            _mm256_storeu_si256(out.as_mut_ptr().add(c + 32) as *mut _, m2);
            _mm256_storeu_si256(out.as_mut_ptr().add(c + 64) as *mut _, m3);
            _mm256_storeu_si256(out.as_mut_ptr().add(c + 96) as *mut _, m4);
            c += std::mem::size_of::<__m256i>() * 4;
        }

        while c + std::mem::size_of::<__m256i>() < db.columns() {
            let mut m1 = _mm256_setzero_si256();

            for i in 0..query.len() {
                if query[i] != b'N' {
                    let q = _mm256_set1_epi8(query[i] as i8);
                    let r1 = _mm256_load_si256(db[i][c..].as_ptr() as *const _);
                    m1 = _mm256_add_epi8(m1, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r1), ones));
                }
            }

            _mm256_storeu_si256(out[c..].as_mut_ptr() as *mut _, m1);
            c += std::mem::size_of::<__m256i>();
        }

        while c < db.columns() {
            let mut m = 0;
            for i in 0..query.len() {
                if query[i] != b'N' && query[i] != db[i][c] {
                    m += 1;
                }
            }
            out[c] = m;
            c += 1;
        }
    }
}

#[derive(Debug)]
pub struct CooBuilder<T> {
    data: RwLock<Vec<(usize, usize, T)>>,
}

impl<T> CooBuilder<T> {
    pub fn new() -> Self {
        Self {
            data: RwLock::new(Vec::new()),
        }
    }

    pub fn insert(&self, i: usize, j: usize, x: T) {
        let mut w = self.data.write().expect("lock was poisoned");
        w.push((i, j, x));
    }

    pub fn len(&self) -> usize {
        let mut r = self.data.read().expect("lock was poisoned");
        r.len()
    }
}

impl<T: Clone> CooBuilder<T> {
    pub fn to_coo_with_dimensions(&self, rows: usize, cols: usize) -> CooMatrix<T> {
        let mut w = self.data.write().expect("lock was poisoned");
        w.sort_by_key(|&(i, j, _)| (i, j));

        let mut coo = CooMatrix::new(rows, cols);
        for (i, j, x) in w.iter() {
            coo.insert(*i, *j, x.clone());
        }

        coo
    }
}

#[derive(Debug)]
pub struct Mapper<'db> {
    pub db: &'db Database,
    pub expected: Vec<CooBuilder<f32>>,
    primer_mismatches: usize,
    kmer_mismatches: usize,
    error_probability: f32,
    primer_region: usize,
    partial_hits: bool,
    pub reads: AtomicUsize,
}

impl<'db> Mapper<'db> {
    pub fn new(db: &'db Database) -> Self {
        let expected = db.regions.iter().map(|_| CooBuilder::new()).collect();
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

    pub fn add(&self, read: Paired<&str>) -> bool {
        // Add a new row to the E_i,h matrices
        let i = self.reads.fetch_add(1, Ordering::Relaxed);

        // Find the best matching primer and primer position
        let (r, pos, primer_mismatches) = self
            .db
            .regions
            .iter()
            .enumerate()
            .map(|(r, region)| {
                (
                    r,
                    self.scan_primer(&region.primer.forward, &read.forward),
                    self.scan_primer(&region.primer.backward.reverse_complement(), &read.backward),
                )
            })
            .min_by(|x, y| (x.1 .1 + x.2 .1).partial_cmp(&(y.1 .1 + y.2 .1)).unwrap())
            .map(|x| {
                (
                    x.0,
                    Paired::new(x.1 .0, x.2 .0),
                    Paired::new(x.1 .1, x.2 .1),
                )
            })
            .unwrap();
        let region = &self.db.regions[r];

        // Skip if primers mismatch the reads
        if primer_mismatches.forward > self.primer_mismatches
            || primer_mismatches.backward > self.primer_mismatches
        {
            return false;
        }

        // Create the kmer pair
        let mut kmer = Paired::new(
            &read.forward[(pos.forward + region.primer.forward.len() as isize) as usize..],
            &read.backward[(pos.backward + region.primer.backward.len() as isize) as usize..],
        );

        // Check that the kmer is long enough for the database regions or that
        // partial mapping is enabled in the mapper.
        if kmer.forward.len() > self.db.k {
            kmer.forward = &kmer.forward[..self.db.k];
        } else if kmer.forward.len() < self.db.k && !self.partial_hits {
            return false;
        }
        if kmer.backward.len() > self.db.k {
            kmer.backward = &kmer.backward[..self.db.k];
        } else if kmer.backward.len() < self.db.k && !self.partial_hits {
            return false;
        }

        // Compute mismatches between the read kmer and all the database kmers
        // let mut mismatch = Paired::<HashMap<usize, u8>>::default();
        // for (x, mm) in region
        //     .trie
        //     .forward
        //     .fuzzy_search(kmer.forward, self.kmer_mismatches)
        // {
        //     let h = region.unique_kmers.forward[x.as_str()];
        //     mismatch.forward.insert(h, mm as u8);
        // }
        // for (x, mm) in region
        //     .trie
        //     .backward
        //     .fuzzy_search(kmer.backward, self.kmer_mismatches)
        // {
        //     let h = region.unique_kmers.backward[x.as_str()];
        //     mismatch.backward.insert(h, mm as u8);
        // }
        let mut mismatch = region
            .block
            .as_ref()
            .map(|matrix| vec![0u8; matrix.columns()]);
        simd_mismatches(
            kmer.forward.as_bytes(),
            &region.block.forward,
            &mut mismatch.forward,
        );
        simd_mismatches(
            kmer.backward.as_bytes(),
            &region.block.backward,
            &mut mismatch.backward,
        );

        // Record the read if it matches any database kmer
        let mut mapped = false;
        for (h, pair) in region.unique_pairs.iter().enumerate() {
            // if let Some(mm_fwd) = mismatch.forward.get(&pair.forward) {
            // if let Some(mm_bwd) = mismatch.backward.get(&pair.backward) {
            // let ne = (mm_fwd + mm_bwd) as usize;
            if mismatch.forward[pair.forward] as usize <= self.kmer_mismatches
                && mismatch.backward[pair.backward] as usize <= self.kmer_mismatches
            {
                let l = kmer.forward.len() + kmer.backward.len();
                let ne =
                    (mismatch.forward[pair.forward] + mismatch.backward[pair.backward]) as usize;
                let e = (self.error_probability / 3.0).powf(ne as f32)
                    * (1.0 - self.error_probability).powf((l - ne) as f32);
                if e > 0.0 {
                    self.expected[r].insert(i, h, e);
                    mapped = true;
                }
            }
            // }
            // }
        }

        mapped
    }

    pub fn finish(self) -> MapperResult {
        // Compute the Q_i,j matrix
        println!("Computing Q matrix");
        let mut q_matrix =
            CooMatrix::<f32>::new(self.reads.load(Ordering::Relaxed), self.db.names.len());
        for (region, expected) in self.db.regions.iter().zip(self.expected) {
            let e = expected.to_coo_with_dimensions(q_matrix.rows(), region.unique_pairs.len());
            let q = e.to_csr().dot(&region.matrix);
            q_matrix = q_matrix + q.to_coo();
        }

        // Compute the pi_j vector
        println!("Computing Pi vector");
        let mut pi = vec![1.0; q_matrix.columns()];
        let mut up = vec![0.0; q_matrix.columns()];
        let mut dens = vec![0.0; q_matrix.rows()];
        for _it in 0..10 {
            // println!("iteration {}", it);
            dens.fill(0.0);
            for (i, j, x) in q_matrix.non_zero_elements() {
                dens[i] += x * pi[j];
            }
            up.fill(0.0);
            for (i, j, x) in q_matrix.non_zero_elements() {
                if dens[i] > 0.0 {
                    up[j] += *x / dens[i]
                }
            }
            for j in 0..q_matrix.columns() {
                pi[j] *= up[j] / q_matrix.rows() as f32;
            }
        }

        println!("Computing X_j vector");
        let mut xj = vec![0.0; q_matrix.columns()];
        for j in 0..q_matrix.columns() {
            if self.db.amplified[j] > 0 {
                xj[j] = pi[j] / self.db.amplified[j] as f32;
            }
        }
        let tot = xj.iter().sum::<f32>();
        if tot > 0.0 {
            for j in 0..q_matrix.columns() {
                xj[j] /= tot;
            }
        }

        MapperResult {
            q: q_matrix,
            pi,
            x: xj,
        }
    }
}

pub struct MapperResult {
    pub q: CooMatrix<f32>,
    pub pi: Vec<f32>,
    pub x: Vec<f32>,
}
