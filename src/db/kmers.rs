use std::ops::Deref;

use serde::Deserialize;
use serde::Serialize;

use crate::matrix::DenseMatrix;
use crate::matrix::MatrixDimensions;

/// A generic data structure to store all the kmers of a database region.
#[derive(Debug, Clone)]
pub struct Kmers {
    block: DenseMatrix<u8>,
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
    use serde::de::Error;
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
