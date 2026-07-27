use std::iter::FusedIterator;

use serde::Deserialize;
use serde::Serialize;

use super::coo::CooMatrix;
use super::csr::CsrMatrix;
use super::MatrixDimensions;
use super::NonZeroElements;
use super::Unique;

// --- CscMatrix ---------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
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

    /// Convert the matrix into COO format without cloning data.
    pub fn into_coo(self) -> CooMatrix<T> {
        let mut coo = CooMatrix::new(self.rows(), self.columns());
        coo.reserve(self.nnz());

        let mut ptr = 0;
        let mut col = 0;
        let mut it = self.data.into_iter();

        while let Some(x) = it.next() {
            while ptr >= self.col_index[col + 1] {
                col += 1;
            }
            ptr += 1;
            coo.i.push(self.row_index[ptr - 1]);
            coo.j.push(col);
            coo.data.push(x);
        }

        coo
    }
}

impl<T: Clone> CscMatrix<T> {
    /// Build a CSR matrix by cloning the data.
    pub fn to_csr(&self) -> CsrMatrix<T> {
        let nnz = self.data.len();

        let mut count = vec![0; self.rows()];
        for j in self.row_index.iter() {
            count[*j] += 1;
        }

        let mut data = self.data.clone();
        let mut col_index = vec![0; nnz];
        let mut row_index = vec![0; self.rows() + 1];
        for j in 1..=self.rows {
            row_index[j] = row_index[j - 1] + count[j - 1];
        }

        for i in 0..self.columns() {
            for src in self.col_index[i]..self.col_index[i + 1] {
                let j = self.row_index[src];
                let dest = row_index[j];
                col_index[dest] = i;
                data[dest] = self.data[src].clone();
                row_index[j] += 1;
            }
        }

        row_index.insert(0, 0);
        row_index.pop();

        CsrMatrix {
            cols: self.columns(),
            data,
            col_index,
            row_index,
        }
    }

    /// Build a COO matrix by cloning data.
    pub fn to_coo(&self) -> CooMatrix<T> {
        self.to_csr().into_coo()
    }

    pub fn select_columns(&self, columns: &[usize]) -> Self {
        let mut data = Vec::new();
        let mut row_index = Vec::new();
        let mut col_index = Vec::with_capacity(columns.len() + 1);

        col_index.push(0);
        for &j in columns.iter() {
            debug_assert!(j < self.columns());
            let col_elem = self.col_index[j]..self.col_index[j + 1];
            row_index.extend_from_slice(&self.row_index[col_elem.clone()]);
            data.extend_from_slice(&self.data[col_elem]);
            col_index.push(row_index.len());
        }

        Self {
            rows: self.rows,
            data,
            row_index,
            col_index,
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

impl<T: Clone> From<&CsrMatrix<T>> for CscMatrix<T> {
    fn from(value: &CsrMatrix<T>) -> Self {
        value.to_csc()
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

    use std::assert_eq;

    use super::super::dok::DokMatrix;
    use super::*;

    #[test]
    fn to_csr() {
        let mut a = DokMatrix::<u8>::new(2, 3);
        a.insert(0, 0, 1);
        a.insert(0, 1, 2);
        a.insert(1, 0, 3);
        a.insert(0, 2, 4);

        let c = a.to_csc().to_csr();
        let mut it = c.non_zero_elements();
        assert_eq!(it.next(), Some((0, 0, &1)));
        assert_eq!(it.next(), Some((0, 1, &2)));
        assert_eq!(it.next(), Some((0, 2, &4)));
        assert_eq!(it.next(), Some((1, 0, &3)));
        assert_eq!(it.next(), None);
    }

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

    #[test]
    fn select_columns() {
        let m1 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(1, 2, 1);
            a.insert(0, 3, 1);
            a.insert(1, 3, 1);
            a.to_csc()
        };
        let b1 = m1.select_columns(&[0, 2]);
        assert_eq!(b1.columns(), 2);
        assert_eq!(b1.rows(), m1.rows());
        let mut it = b1.non_zero_elements();
        assert_eq!(it.next(), Some((1, 1, &1)));
        assert_eq!(it.next(), None);

        let m2 = {
            let mut a = DokMatrix::<u8>::new(2, 4);
            a.insert(0, 1, 1);
            a.insert(1, 2, 1);
            a.insert(0, 3, 1);
            a.insert(1, 3, 1);
            a.to_csc()
        };
        let b2 = m2.select_columns(&[0, 1, 2, 3]);
        assert_eq!(b2.columns(), m2.columns());
        assert_eq!(b2.rows(), m2.rows());
        assert_eq!(b2, m2);
    }
}
