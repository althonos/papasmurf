use std::ops::Deref;

use crate::errors::Error;
use crate::matrix::DenseMatrix;
use crate::matrix::MatrixDimensions;

/// A generic data structure to store all the kmers of a database region.
#[derive(Debug, Clone)]
pub struct Kmers {
    pub block: DenseMatrix<u8>,
}

impl Kmers {
    pub fn mismatches(&self, query: &str) -> Result<Vec<u8>, Error> {
        crate::seq::validate(query)?;
        
        let q = query.as_bytes();
        let b = &self.block;
        let mut out = vec![0u8; b.columns() ];

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if std::is_x86_feature_detected!("avx2") {
            unsafe { self::avx2::mismatches(q, b, out.as_mut()); }
            return Ok(out);
        }

        unimplemented!()
    }
}

impl From<DenseMatrix<u8>> for Kmers {
    fn from(block: DenseMatrix<u8>) -> Self {
        Kmers { block }
    }
}

impl Deref for Kmers {
    type Target = DenseMatrix<u8>;
    fn deref(&self) -> &Self::Target {
        &self.block
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod avx2 {

    use super::*;

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    #[target_feature(enable = "avx2")]
    pub unsafe fn mismatches(query: &[u8], kmers: &DenseMatrix<u8>, out: &mut [u8]) {
        let ones = _mm256_set1_epi8(1);
        
        let mut c = 0;
        while c + std::mem::size_of::<__m256i>() * 4 < kmers.columns() {
            let mut m1 = _mm256_setzero_si256();
            let mut m2 = _mm256_setzero_si256();
            let mut m3 = _mm256_setzero_si256();
            let mut m4 = _mm256_setzero_si256();

            for i in 0..query.len() {
                let kmerptr = kmers[i].as_ptr();
                if query[i] != b'N' {
                    let q = _mm256_set1_epi8(query[i] as i8);
                    let r1 = _mm256_load_si256(kmerptr.add(c) as *const _);
                    let r2 = _mm256_load_si256(kmerptr.add(c + 32) as *const _);
                    let r3 = _mm256_load_si256(kmerptr.add(c + 64) as *const _);
                    let r4 = _mm256_load_si256(kmerptr.add(c + 96) as *const _);
                    m1 = _mm256_add_epi8(m1, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r1), ones));
                    m2 = _mm256_add_epi8(m2, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r2), ones));
                    m3 = _mm256_add_epi8(m3, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r3), ones));
                    m4 = _mm256_add_epi8(m4, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r4), ones));
                }
            }

            let outptr = out.as_mut_ptr();
            _mm256_storeu_si256(outptr.add(c) as *mut _, m1);
            _mm256_storeu_si256(outptr.add(c + 32) as *mut _, m2);
            _mm256_storeu_si256(outptr.add(c + 64) as *mut _, m3);
            _mm256_storeu_si256(outptr.add(c + 96) as *mut _, m4);
            c += std::mem::size_of::<__m256i>() * 4;
        }
        while c + std::mem::size_of::<__m256i>() < kmers.columns() {
            let mut m1 = _mm256_setzero_si256();

            for i in 0..query.len() {
                if query[i] != b'N' {
                    let q = _mm256_set1_epi8(query[i] as i8);
                    let r1 = _mm256_load_si256(kmers[i][c..].as_ptr() as *const _);
                    m1 = _mm256_add_epi8(m1, _mm256_andnot_si256(_mm256_cmpeq_epi8(q, r1), ones));
                }
            }

            _mm256_storeu_si256(out[c..].as_mut_ptr() as *mut _, m1);
            c += std::mem::size_of::<__m256i>();
        }
        while c < kmers.columns() {
            let mut m = 0;
            for i in 0..query.len() {
                if query[i] != b'N' && query[i] != kmers[i][c] {
                    m += 1;
                }
            }
            out[c] = m;
            c += 1;
        }
    }
}

#[cfg(feature = "serde")]
mod ser {

    use super::*;
    use serde::ser::SerializeSeq;
    use serde::ser::Serializer;
    use serde::Serialize;

    impl Serialize for Kmers {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            let mut seq = serializer.serialize_seq(Some(self.block.rows()))?;
            for i in 0..self.block.rows() {
                let row = std::str::from_utf8(&self.block[i])
                    .expect("kmer should always be an ASCII string");
                seq.serialize_element(row)?
            }
            seq.end()
        }
    }
}

#[cfg(feature = "serde")]
mod de {

    use super::*;

    use std::fmt::Formatter;
    use std::fmt::Result as FmtResult;

    use serde::de::Deserializer;

    use serde::de::SeqAccess;
    use serde::de::Visitor;
    use serde::Deserialize;

    struct KmersVisitor;

    impl<'de> Visitor<'de> for KmersVisitor {
        type Value = Kmers;

        fn expecting(&self, formatter: &mut Formatter) -> FmtResult {
            write!(formatter, "a sequence of strings")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let mut buffer = Vec::<String>::with_capacity(seq.size_hint().unwrap_or(0));
            while let Some(x) = seq.next_element()? {
                buffer.push(x);
            }

            let rows = buffer.len();
            let cols = buffer.first().map(|row| row.len()).unwrap_or(0);
            let mut matrix = DenseMatrix::new(rows, cols);
            for (i, row) in buffer.into_iter().enumerate() {
                matrix[i].copy_from_slice(row.as_bytes());
            }

            Ok(matrix.into())
        }
    }

    impl<'de> Deserialize<'de> for Kmers {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            Ok(deserializer.deserialize_seq(KmersVisitor)?)
        }
    }
}
