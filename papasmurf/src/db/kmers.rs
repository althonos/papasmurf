use std::ops::Deref;

use crate::errors::Error;
use crate::matrix::DenseMatrix;
use crate::matrix::MatrixDimensions;

/// A generic data structure to store all the kmers of a database region.
#[derive(Debug, Clone)]
pub struct Kmers {
    block: DenseMatrix<u8>,
}

impl Kmers {
    /// Create a new k-mer storage for the given sequences.
    pub fn new<I>(kmers: I) -> Result<Self, Error>
    where
        I: IntoIterator,
        <I as IntoIterator>::Item: AsRef<[u8]>,
        <I as IntoIterator>::IntoIter: ExactSizeIterator,
    {
        let mut it = kmers.into_iter().peekable();
        let rows = it.peek().map(|kmer| kmer.as_ref().len()).unwrap_or(0);
        let cols = it.len();
        let mut block = DenseMatrix::new(rows, cols);

        for (i, item) in it.enumerate() {
            if i > cols {
                return Err(Error::InvalidDimensions);
            }
            let seq = item.as_ref();
            for (j, b) in seq.iter().enumerate() {
                if j > rows {
                    return Err(Error::InvalidDimensions);
                }
                match *b {
                    b'A' | b'C' | b'G' | b'T' => block[j][i] = *b,
                    _ => return Err(Error::InvalidDna),
                }
            }
        }

        Ok(Self::from(block))
    }

    /// Compute the number of mismatches between all k-mers and the query.
    #[allow(unreachable_code)]
    pub fn mismatches(&self, query: &str) -> Result<Vec<u8>, Error> {
        crate::seq::validate(query)?;

        let q = query.as_bytes();
        let b = &self.block;
        let mut out = vec![0u8; b.columns()];

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if std::is_x86_feature_detected!("avx2") {
            unsafe {
                self::avx2::mismatches(q, b, out.as_mut());
            }
            return Ok(out);
        }

        #[cfg(target_arch = "x86")]
        if std::is_x86_feature_detected!("sse2") {
            unsafe {
                self::sse2::mismatches(q, b, out.as_mut());
            }
            return Ok(out);
        }

        #[cfg(target_arch = "x86_64")]
        {
            unsafe {
                self::sse2::mismatches(q, b, out.as_mut());
            }
            return Ok(out);
        }

        for c in 0..self.block.columns() {
            let mut m = 0;
            for i in 0..q.len() {
                if q[i] != b'N' && q[i] != self.block[i][c] {
                    m += 1;
                }
            }
            out[c] = m;
        }
        Ok(out)
    }
}

impl MatrixDimensions for Kmers {
    fn rows(&self) -> usize {
        self.block.rows()
    }

    fn columns(&self) -> usize {
        self.block.columns()
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

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod sse2 {

    use super::*;

    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    #[target_feature(enable = "sse2")]
    pub unsafe fn mismatches(query: &[u8], kmers: &DenseMatrix<u8>, out: &mut [u8]) {
        let ones = _mm_set1_epi8(1);

        let mut c = 0;
        while c + std::mem::size_of::<__m128i>() * 4 < kmers.columns() {
            let mut m1 = _mm_setzero_si128();
            let mut m2 = _mm_setzero_si128();
            let mut m3 = _mm_setzero_si128();
            let mut m4 = _mm_setzero_si128();

            for i in 0..query.len() {
                let kmerptr = kmers[i].as_ptr();
                if query[i] != b'N' {
                    let q = _mm_set1_epi8(query[i] as i8);
                    let r1 = _mm_load_si128(kmerptr.add(c) as *const _);
                    let r2 = _mm_load_si128(kmerptr.add(c + 16) as *const _);
                    let r3 = _mm_load_si128(kmerptr.add(c + 32) as *const _);
                    let r4 = _mm_load_si128(kmerptr.add(c + 48) as *const _);
                    m1 = _mm_add_epi8(m1, _mm_andnot_si128(_mm_cmpeq_epi8(q, r1), ones));
                    m2 = _mm_add_epi8(m2, _mm_andnot_si128(_mm_cmpeq_epi8(q, r2), ones));
                    m3 = _mm_add_epi8(m3, _mm_andnot_si128(_mm_cmpeq_epi8(q, r3), ones));
                    m4 = _mm_add_epi8(m4, _mm_andnot_si128(_mm_cmpeq_epi8(q, r4), ones));
                }
            }

            let outptr = out.as_mut_ptr();
            _mm_storeu_si128(outptr.add(c) as *mut _, m1);
            _mm_storeu_si128(outptr.add(c + 16) as *mut _, m2);
            _mm_storeu_si128(outptr.add(c + 32) as *mut _, m3);
            _mm_storeu_si128(outptr.add(c + 48) as *mut _, m4);
            c += std::mem::size_of::<__m128i>() * 4;
        }
        while c + std::mem::size_of::<__m128i>() < kmers.columns() {
            let mut m1 = _mm_setzero_si128();

            for i in 0..query.len() {
                if query[i] != b'N' {
                    let q = _mm_set1_epi8(query[i] as i8);
                    let r1 = _mm_load_si128(kmers[i][c..].as_ptr() as *const _);
                    m1 = _mm_add_epi8(m1, _mm_andnot_si128(_mm_cmpeq_epi8(q, r1), ones));
                }
            }

            _mm_storeu_si128(out[c..].as_mut_ptr() as *mut _, m1);
            c += std::mem::size_of::<__m128i>();
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

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_new() {
        let kmers = Kmers::new(&["ATTAT", "ATGCA"]).unwrap();
        assert_eq!(kmers.block.rows(), 5);
        assert_eq!(kmers.block.columns(), 2);
        assert_eq!(&kmers.block[0], b"AA");
        assert_eq!(&kmers.block[1], b"TT");
        assert_eq!(&kmers.block[2], b"TG");
        assert_eq!(&kmers.block[3], b"AC");
        assert_eq!(&kmers.block[4], b"TA");
    }

    #[test]
    fn test_mismatches() {
        let kmers = Kmers::new(&[
            "AAACA", "AACAC", "AACGG", "AAGCA", "AATCC", "ACAAG", "ACAGG", "ACCCA", "ACTCA",
            "ACTGC", "AGACG", "AGTAA", "AGTTC", "ATCAC", "ATTAG", "CACAC", "CACAG", "CACTA",
            "CAGAC", "CATTC", "CCCAT", "CCCTA", "CGTGC", "CTAAT", "CTACT", "CTAGT", "CTCCG",
            "CTGAA", "CTGGA", "CTTCT", "CTTTC", "GAAGA", "GAATC", "GAGGT", "GGACC", "GGGGA",
            "GTTCA", "TAGCG", "TCCTA", "TCTCA", "TCTTG", "TGTTA", "TGTTC", "TTAAC", "TTCAA",
            "TTCAC", "TTCGA", "TTTCT",
        ])
        .unwrap();

        let mm = kmers.mismatches("CGTGC").unwrap();
        assert_eq!(mm[22], 0);

        let mm = kmers.mismatches("AAGCA").unwrap();
        assert_eq!(mm[0], 1);
        assert_eq!(mm[1], 3);
        assert_eq!(mm[2], 3);
        assert_eq!(mm[3], 0);
    }
}
