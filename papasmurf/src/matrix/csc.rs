use std::iter::FusedIterator;

use serde::Deserialize;
use serde::Serialize;

use super::MatrixDimensions;
use super::NonZeroElements;
use super::Unique;

// --- CscMatrix ---------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CscMatrix<T> {
    pub(super) rows: usize,
    pub(super) data: Vec<T>,
    pub(super) col_index: Vec<usize>,
    pub(super) row_index: Vec<usize>,
}

impl<T> CscMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        CscMatrix {
            rows,
            data: Vec::new(),
            row_index: Vec::new(),
            col_index: vec![0; cols + 1],
        }
    }
}

impl<T: PartialEq> CscMatrix<T> {
    /// Find the indices of unique columns in the CSC matrix.
    pub fn unique_columns(&self) -> Unique {
        // Compute number of nonzero elements per columns
        let mut lengths: Vec<usize> = vec![0; self.columns()];
        for i in 1..=self.columns() {
            lengths[i - 1] = self.col_index[i] - self.col_index[i - 1];
        }

        // Sort columns by number of nonzero elements
        let mut colptr = (0..self.columns()).collect::<Vec<_>>();
        colptr.sort_by_key(|&j| (lengths[j], j));

        // Prepare the result vector
        let mut indices = Vec::new();
        let mut reverse = vec![None; self.columns()];

        // Iterate by groups
        let mut i = 0;
        while i < colptr.len() {
            // Advance until a larger column is reached
            let mut i2 = i;
            while i2 < colptr.len() && lengths[colptr[i2]] == lengths[colptr[i]] {
                i2 += 1;
            }

            // Get group of indices with equal length
            let group = &colptr[i..i2];
            debug_assert!(group.iter().all(|&j| lengths[j] == lengths[group[0]]));
            i = i2;

            // Zero-sized columns are duplicates
            if lengths[group[0]] == 0 {
                let j1 = group[0];
                for &j2 in group.iter().skip(1) {
                    reverse[j2] = Some(indices.len());
                }
                reverse[j1] = Some(indices.len());
                indices.push(j1);
                continue;
            }

            // Compare columns pair by pair
            for (n, &j1) in group.iter().enumerate() {
                if reverse[j1].is_some() {
                    continue;
                }
                for &j2 in group[n + 1..].iter() {
                    if reverse[j2].is_some() {
                        continue;
                    }
                    let r1 = &self.row_index[self.col_index[j1]..self.col_index[j1 + 1]];
                    let r2 = &self.row_index[self.col_index[j2]..self.col_index[j2 + 1]];
                    let d1 = &self.data[self.col_index[j1]..self.col_index[j1 + 1]];
                    let d2 = &self.data[self.col_index[j2]..self.col_index[j2 + 1]];
                    if r1 == r2 && d1 == d2 {
                        reverse[j2] = Some(indices.len());
                    }
                }
                reverse[j1] = Some(indices.len());
                indices.push(j1);
            }
        }

        Unique {
            indices,
            reverse: reverse.into_iter().map(Option::unwrap).collect(),
        }
    }

    /// Return a boolean vector denoting duplicate columns.
    pub fn duplicate_columns(&self) -> Vec<bool> {
        let uniq = self.unique_columns();
        let mut dup = vec![true; self.columns()];
        for &j in uniq.indices.iter() {
            dup[j] = false;
        }
        dup
    }
}

impl<T> Default for CscMatrix<T> {
    fn default() -> Self {
        Self {
            rows: 0,
            data: Vec::new(),
            col_index: Vec::new(),
            row_index: Vec::new(),
        }
    }
}

impl<T> MatrixDimensions for CscMatrix<T> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn columns(&self) -> usize {
        self.col_index.len() - 1
    }
}

// --- NonZeroIter -------------------------------------------------------------

pub struct NonZeroIter<'m, T> {
    matrix: &'m CscMatrix<T>,
    col: usize,
    ptr: usize,
}

impl<'mx, T> Iterator for NonZeroIter<'mx, T> {
    type Item = (usize, usize, &'mx T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr >= self.matrix.data.len() {
            return None;
        }
        while self.ptr >= self.matrix.col_index[self.col + 1] {
            if self.col + 1 > self.matrix.col_index.len() {
                return None;
            }
            self.col += 1;
        }
        self.ptr += 1;
        Some((
            self.matrix.row_index[self.ptr - 1],
            self.col,
            &self.matrix.data[self.ptr - 1],
        ))
    }
}

impl<'mx, T> ExactSizeIterator for NonZeroIter<'mx, T> {
    fn len(&self) -> usize {
        self.matrix.data.len() - self.ptr
    }
}

impl<'mx, T> FusedIterator for NonZeroIter<'mx, T> {}

impl<'m, T: 'm> NonZeroElements<'m, T> for CscMatrix<T> {
    type Iter = NonZeroIter<'m, T>;
    fn non_zero_elements(&'m self) -> Self::Iter {
        NonZeroIter {
            col: 0,
            ptr: 0,
            matrix: self,
        }
    }
}

#[cfg(test)]
mod test {

    use super::super::dok::DokMatrix;
    use super::*;

    #[test]
    fn non_zero_elements() {
        let m1 = CscMatrix::<u8>::new(2, 2);
        let mut it = m1.non_zero_elements();
        assert_eq!(it.next(), None);

        let mut a = DokMatrix::<u8>::new(2, 2);
        a.insert(0, 0, 1);
        a.insert(0, 1, 2);
        a.insert(1, 0, 3);
        let m2 = a.to_csc();

        let mut it = m2.non_zero_elements();
        assert_eq!(it.next(), Some((0, 0, &1)));
        assert_eq!(it.next(), Some((1, 0, &3)));
        assert_eq!(it.next(), Some((0, 1, &2)));
        assert_eq!(it.next(), None);
    }

    #[test]
    fn unique_columns() {
        let m1 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(1, 2, 1);
            a.insert(0, 3, 1);
            a.insert(1, 3, 1);
            a.to_csc()
        };
        let u1 = m1.unique_columns();
        assert_eq!(u1.indices, &[0, 1, 2, 3]);
        assert_eq!(u1.reverse, &[0, 1, 2, 3]);

        let m2 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(0, 3, 1);
            a.to_csc()
        };
        let u2 = m2.unique_columns();
        assert_eq!(u2.indices, &[0, 1]);
        assert_eq!(u2.reverse, &[0, 1, 0, 1]);

        let m3 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(1, 1, 1);
            a.insert(1, 2, 1);
            a.insert(0, 3, 1);
            a.insert(1, 3, 1);
            a.to_csc()
        };
        let u3 = m3.unique_columns();
        assert_eq!(u3.indices, &[0, 2, 1]);
        assert_eq!(u3.reverse, &[0, 2, 1, 2]);

        let m4 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(0, 3, 2);
            a.to_csc()
        };
        let u4 = m4.unique_columns();
        assert_eq!(u4.indices, &[0, 1, 3]);
        assert_eq!(u4.reverse, &[0, 1, 0, 2]);
    }

    #[test]
    fn duplicate_columns() {
        let m1 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(1, 2, 1);
            a.insert(0, 3, 1);
            a.insert(1, 3, 1);
            a.to_csc()
        };
        assert_eq!(m1.duplicate_columns(), vec![false, false, false, false]);

        let m2 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(0, 3, 1);
            a.to_csc()
        };
        assert_eq!(m2.duplicate_columns(), vec![false, false, true, true]);

        let m3 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(1, 1, 1);
            a.insert(1, 2, 1);
            a.insert(0, 3, 1);
            a.insert(1, 3, 1);
            a.to_csc()
        };
        assert_eq!(m3.duplicate_columns(), vec![false, false, false, true]);

        let m4 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(0, 3, 2);
            a.to_csc()
        };
        assert_eq!(m4.duplicate_columns(), vec![false, false, true, false]);
    }
}
