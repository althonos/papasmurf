use super::db::Database;

use super::matrix::DenseMatrix;
use super::matrix::DokMatrix;
use super::matrix::MatrixDimensions;
use super::matrix::NonZeroElements;
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
    pub expected: Vec<DokMatrix<f32>>,
}

impl<'db> Mapper<'db> {
    pub fn new(db: &'db Database) -> Self {
        let expected = db
            .regions
            .iter()
            .map(|region| DokMatrix::new(0, region.unique_pairs.len()))
            .collect();
        Self { expected, db }
    }

    pub fn add(&mut self, read: Paired<&str>) -> bool {
        let i = self.expected[0].rows();

        // Add a new row to the E_i,h matrices
        for e in self.expected.iter_mut() {
            e.grow(1, 0);
        }

        //
        let (r, pos, primer_mismatches) =
            self.db
                // let (r, pos) = db
                .regions
                .iter()
                .enumerate()
                // .map(|(r, region)| {
                //     pli.score_into(&striped.forward, &region.profile.forward, &mut scores);
                //     let fwd_pos = pli.argmax(&scores).unwrap();
                //     let fwd_score = scores[fwd_pos];
                //     pli.score_into(&striped.backward, &region.profile.backward, &mut scores);
                //     let bwd_pos = pli.argmax(&scores).unwrap();
                //     let bwd_score = scores[bwd_pos];
                //     (r, Paired::new((fwd_pos, fwd_score), (bwd_pos, bwd_score)))
                // })
                // .max_by(|(_, p1), (_, p2)| {
                //     (p1.forward.1 + p2.backward.1)
                //         .partial_cmp(&(p2.forward.1 + p2.backward.1))
                //         .unwrap()
                // })
                // .map(|(r, p)| (r, p.map(|x| x.0)))
                // .unwrap();
                .map(|(r, region)| {
                    (
                        r,
                        (0..read.forward.len() - region.primer.forward.len())
                            .map(|i| {
                                (
                                    i,
                                    region.primer.forward.mismatches(
                                        &read.forward[i..i + region.primer.forward.len()],
                                    ),
                                )
                            })
                            .min_by_key(|(_, s)| *s)
                            .unwrap(),
                        (0..read.backward.len() - region.primer.backward.len())
                            .map(|i| {
                                (
                                    i,
                                    region.primer.backward.reverse_complement().mismatches(
                                        &read.backward[i..i + region.primer.backward.len()],
                                    ),
                                )
                            })
                            .min_by_key(|(_, s)| *s)
                            .unwrap(),
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

        if primer_mismatches.forward > 2 || primer_mismatches.backward > 2 {
            return false;
        }
        let mut kmer = Paired::new(
            &read.forward[pos.forward + self.db.regions[r].primer.forward.len()..],
            &read.backward[pos.backward + self.db.regions[r].primer.backward.len()..],
        );
        // let mut kmer = Paired::new(
        //     &seq.forward[pos.forward..],
        //     &seq.backward[pos.backward..]
        // );
        if kmer.forward.len() > self.db.k {
            kmer.forward = &kmer.forward[..self.db.k];
        }
        if kmer.backward.len() > self.db.k {
            kmer.backward = &kmer.backward[..self.db.k];
        }

        // if pos.forward + builder.k > read.forward.seq().len() {
        //     println!("{:?}", &pos.backward);
        //     continue;
        // }
        // if pos.backward + builder.k > read.backward.seq().len() {
        //     println!("{:?}", &pos.backward);
        //     continue;
        // }

        // let mut kmer = Paired::new(
        //     if pos.forward + builder.k > read.forward.seq().len() {
        //         &read.forward.seq()[pos.forward..]
        //     } else {
        //         &read.forward.seq()[pos.forward..pos.forward + builder.k]
        //     },
        //     if pos.backward + builder.k > read.backward.seq().len() {
        //         &read.backward.seq()[pos.backward..]
        //     } else {
        //         &read.backward.seq()[pos.backward..pos.backward + builder.k]
        //     }
        // );

        let mut mismatch = self.db.regions[r]
            .unique_kmers
            .as_ref()
            .map(|block| vec![0u8; block.columns()]);
        simd_mismatches(
            kmer.forward.as_bytes(),
            &self.db.regions[r].unique_kmers.forward,
            &mut mismatch.forward,
        );
        simd_mismatches(
            kmer.backward.as_bytes(),
            &self.db.regions[r].unique_kmers.backward,
            &mut mismatch.backward,
        );

        let mut mapped = false;
        for (h, pair) in self.db.regions[r].unique_pairs.iter().enumerate() {
            let mm = Paired::new(
                mismatch.forward[pair.forward],
                mismatch.backward[pair.backward],
            );
            const PE: f32 = 0.005;
            let ne = (mm.forward + mm.backward) as f32;
            let l = kmer.forward.len() + kmer.backward.len();
            let e = (PE / 3.0).powf(ne) * (1.0 - PE).powf(l as f32 - ne);
            if e > 0.0 && ne <= 2.0 {
                self.expected[r].insert(i, h, e);
                mapped = true;
            }
        }

        mapped
    }
}
