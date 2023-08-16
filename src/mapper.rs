use std::collections::HashMap;

use super::db::Database;
use super::db::KmerTrie;
use super::matrix::DenseMatrix;
use super::matrix::DokMatrix;
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

#[derive(Debug, Clone)]
pub struct Mapper<'db> {
    pub db: &'db Database,
    pub tries: Vec<Paired<KmerTrie>>,    
    pub expected: Vec<DokMatrix<f32>>,
    primer_mismatches: usize,
    kmer_mismatches: usize,
    error_probability: f32,
    primer_region: usize,
}

impl<'db> Mapper<'db> {
    pub fn new(db: &'db Database) -> Self {
        let expected = db
            .regions
            .iter()
            .map(|region| DokMatrix::new(0, region.unique_pairs.len()))
            .collect();
        let tries = db.regions.iter()
            .map(|region| region.unique_kmers.as_ref().map(|kmers| {
                let mut trie = KmerTrie::new(db.k);
                for kmer in kmers.iter() {
                    trie.insert(kmer);
                }
                trie
            }))
            .collect();
        Self {
            expected,
            db,
            tries,
            primer_mismatches: 2,
            kmer_mismatches: 2,
            error_probability: 0.005,
            primer_region: 30,
        }
    }

    fn scan(&self, primer: &Primer, sequence: &str) -> (usize, usize) {
        let offset = self.primer_region.min(sequence.len() - primer.len());
        let mut min_i = usize::MAX;
        let mut min_mm = usize::MAX;
        for i in 0..offset {
            let mm = primer.mismatches(&sequence[i..i + primer.len()]);
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

    pub fn add(&mut self, read: Paired<&str>) -> bool {
        // Add a new row to the E_i,h matrices
        let i = self.expected[0].rows();
        for e in self.expected.iter_mut() {
            e.grow(1, 0);
        }

        // Find the best matching primer and primer position
        let (r, pos, primer_mismatches) = self
            .db
            .regions
            .iter()
            .enumerate()
            .map(|(r, region)| {
                (
                    r,
                    self.scan(&region.primer.forward, &read.forward),
                    self.scan(&region.primer.backward.reverse_complement(), &read.backward),
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
        if primer_mismatches.forward > 2 || primer_mismatches.backward > 2 {
            return false;
        }

        // Create the kmer pair
        let mut kmer = Paired::new(
            &read.forward[pos.forward + region.primer.forward.len()..],
            &read.backward[pos.backward + region.primer.backward.len()..],
        );
        if kmer.forward.len() > self.db.k {
            kmer.forward = &kmer.forward[..self.db.k];
        }
        if kmer.backward.len() > self.db.k {
            kmer.backward = &kmer.backward[..self.db.k];
        }

        // Find the k-mers with few mismatches
        // let mm_bwd = self.tries[r].backward.fuzzy_search(kmer.backward, self.kmer_mismatches);

        // // Compute mismatches between the read kmer and all the database kmers
        let mut mismatch = Paired::<HashMap<usize, u8>>::default();
        for (x, mm) in self.tries[r].forward.fuzzy_search(kmer.forward, self.kmer_mismatches) {
            let h = region.unique_kmers.forward[&crate::utils::Rc::from(x)];
            mismatch.forward.insert(h, mm as u8);
        }
        for (x, mm) in self.tries[r].backward.fuzzy_search(kmer.backward, self.kmer_mismatches) {
            let h = region.unique_kmers.backward[&crate::utils::Rc::from(x)];
            mismatch.backward.insert(h, mm as u8);
        }

        // simd_mismatches(
        //     kmer.forward.as_bytes(),
        //     &region.unique_kmers.forward,
        //     &mut mismatch.forward,
        // );
        // simd_mismatches(
        //     kmer.backward.as_bytes(),
        //     &region.unique_kmers.backward,
        //     &mut mismatch.backward,
        // );

        // Record the read if it matches any database kmer
        let mut mapped = false;
        for (h, pair) in region.unique_pairs.iter().enumerate() {
            if let (Some(mm_fwd), Some(mm_bwd)) = (mismatch.forward.get(&pair.forward), mismatch.backward.get(&pair.backward)) {
                let ne = (mm_fwd + mm_bwd) as usize;
                let l = kmer.forward.len() + kmer.backward.len();
                let e = (self.error_probability / 3.0).powf(ne as f32)
                    * (1.0 - self.error_probability).powf((l - ne) as f32);
                if e > 0.0 && ne <= self.kmer_mismatches {
                    self.expected[r].insert(i, h, e);
                    mapped = true;
                }
            }
        }

        mapped
    }
}
