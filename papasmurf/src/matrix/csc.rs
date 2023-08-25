use std::iter::FusedIterator;

use serde::Deserialize;
use serde::Serialize;

use super::MatrixDimensions;
use super::NonZeroElements;

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
}
