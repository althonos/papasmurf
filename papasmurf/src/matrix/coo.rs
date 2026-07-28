use std::assert_eq;
use std::cmp::Ordering;
use std::iter::FusedIterator;
use std::ops::Add;
use std::ops::Range;

use serde::Deserialize;
use serde::Serialize;

use super::csc::CscMatrix;
use super::csr::CsrMatrix;
use super::MatrixDimensions;
use super::NonZeroElements;
use super::NonZeroElementsMut;
use super::VerticalStack;

// --- CooMatrix ---------------------------------------------------------------

/// A sparse matrix in coordinate (COO) format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CooMatrix<T> {
    pub(super) rows: usize,
    pub(super) cols: usize,
    pub(super) i: Vec<usize>,
    pub(super) j: Vec<usize>,
    pub(super) data: Vec<T>,
}

impl<T> CooMatrix<T> {
    /// Create a new COO matrix with the given dimensions.
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            i: Vec::new(),
            j: Vec::new(),
            data: Vec::new(),
        }
    }

    /// Reserve space for the given number of non-zero elements.
    pub fn reserve(&mut self, nnz: usize) {
        self.i.reserve(nnz);
        self.j.reserve(nnz);
        self.data.reserve(nnz);
    }

    /// Convert the matrix into CSR format without cloning data.
    pub fn into_csr(self) -> CsrMatrix<T> {
        let mut csr = CsrMatrix::new(self.rows, self.cols);
        let mut it = self.i.into_iter().zip(self.j).zip(self.data).peekable();

        for i in 0..self.rows {
            csr.row_index[i] = csr.col_index.len();
            while let Some(((x, _), _)) = it.peek() {
                if *x != i {
                    break;
                }
                let ((_, y), z) = it.next().unwrap();
                csr.col_index.push(y);
                csr.data.push(z);
            }
        }

        csr.row_index[self.rows] = csr.col_index.len();
        csr
    }

    /// Insert a new non-zero element at the end of the matrix.
    pub(super) fn insert(&mut self, i: usize, j: usize, data: T) {
        assert!(i < self.rows, "row index out of range");
        assert!(j < self.cols, "column index out of range");

        if self.nnz() > 0 {
            let last_i = *self.i.last().unwrap();
            let last_j = *self.j.last().unwrap();
            assert!(
                (i, j) > (last_i, last_j),
                "{:?} > {:?}",
                (i, j),
                (last_i, last_j)
            );
        }

        self.i.push(i);
        self.j.push(j);
        self.data.push(data);
    }
}

impl<T: Clone> CooMatrix<T> {
    /// Build a CSR matrix by cloning data.
    pub fn to_csr(&self) -> CsrMatrix<T> {
        let mut csr = CsrMatrix::new(self.rows, self.cols);
        let mut it = self.non_zero_elements().peekable();

        for i in 0..self.rows {
            csr.row_index[i] = csr.col_index.len();
            while let Some((x, _, _)) = it.peek() {
                if *x != i {
                    break;
                }
                let (_, y, z) = it.next().unwrap();
                csr.col_index.push(y);
                csr.data.push(z.clone());
            }
        }

        csr.row_index[self.rows] = csr.col_index.len();
        csr
    }
}

impl<T> MatrixDimensions for CooMatrix<T> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn columns(&self) -> usize {
        self.cols
    }
}

impl<T, U> Add<&CooMatrix<U>> for CooMatrix<T>
where
    T: Add<Output = T> + Clone + Default + PartialEq,
    U: Into<T> + Clone,
{
    type Output = CooMatrix<T>;
    fn add(self, rhs: &CooMatrix<U>) -> Self {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);

        let zero = T::default();

        let mut out = Self::new(self.rows, rhs.cols);
        let mut x = 0;
        let mut y = 0;

        while x < self.data.len() && y < rhs.data.len() {
            let i1 = self.i[x];
            let j1 = self.j[x];
            let d1 = &self.data[x];
            let i2 = rhs.i[y];
            let j2 = rhs.j[y];
            let d2 = &rhs.data[y];
            match (i1, j1).cmp(&(i2, j2)) {
                Ordering::Less => {
                    out.insert(i1, j1, d1.clone());
                    x += 1;
                }
                Ordering::Greater => {
                    out.insert(i2, j2, d2.clone().into());
                    y += 1;
                }
                Ordering::Equal => {
                    let d = d1.clone() + d2.clone().into();
                    if d != zero {
                        out.insert(i1, j1, d1.clone() + d2.clone().into());
                    }
                    x += 1;
                    y += 1;
                }
            }
        }

        while x < self.data.len() {
            out.i.push(self.i[x]);
            out.j.push(self.j[x]);
            out.data.push(self.data[x].clone());
            x += 1;
        }
        while y < rhs.data.len() {
            out.i.push(rhs.i[y]);
            out.j.push(rhs.j[y]);
            out.data.push(rhs.data[y].clone().into());
            y += 1;
        }

        out
    }
}

impl<T, U> Add<CooMatrix<U>> for CooMatrix<T>
where
    T: Add<Output = T> + Clone + Default + PartialEq,
    U: Into<T> + Clone,
{
    type Output = CooMatrix<T>;
    fn add(self, rhs: CooMatrix<U>) -> Self::Output {
        self.add(&rhs)
    }
}

impl<T: Clone> From<&CsrMatrix<T>> for CooMatrix<T> {
    fn from(csr: &CsrMatrix<T>) -> CooMatrix<T> {
        csr.to_coo()
    }
}

impl<T> From<CsrMatrix<T>> for CooMatrix<T> {
    fn from(csr: CsrMatrix<T>) -> CooMatrix<T> {
        csr.into_coo()
    }
}

impl<T: Clone> From<&CscMatrix<T>> for CooMatrix<T> {
    fn from(csc: &CscMatrix<T>) -> CooMatrix<T> {
        csc.to_coo()
    }
}

impl<T> From<CscMatrix<T>> for CooMatrix<T> {
    fn from(csc: CscMatrix<T>) -> CooMatrix<T> {
        csc.into_coo()
    }
}

impl<T> VerticalStack<&CooMatrix<T>> for CooMatrix<T>
where
    T: Clone,
{
    fn vstack(&mut self, other: &CooMatrix<T>) {
        assert_eq!(self.columns(), other.columns());
        let offset = self.rows();
        for (i, j, x) in other.non_zero_elements() {
            self.insert(i + offset, j, x.clone());
        }
    }
}

impl<T> VerticalStack<CooMatrix<T>> for CooMatrix<T>
where
    T: Clone,
{
    fn vstack(&mut self, other: CooMatrix<T>) {
        self.vstack(&other);
    }
}

// --- NonZeroIter -------------------------------------------------------------

pub struct NonZeroIter<'m, T> {
    matrix: &'m CooMatrix<T>,
    pos: Range<usize>,
}

impl<'mx, T> Iterator for NonZeroIter<'mx, T> {
    type Item = (usize, usize, &'mx T);
    fn next(&mut self) -> Option<Self::Item> {
        let pos = self.pos.next()?;
        Some((
            self.matrix.i[pos],
            self.matrix.j[pos],
            &self.matrix.data[pos],
        ))
    }
}

impl<'mx, T> ExactSizeIterator for NonZeroIter<'mx, T> {
    fn len(&self) -> usize {
        self.pos.len()
    }
}

impl<'mx, T> FusedIterator for NonZeroIter<'mx, T> {}

impl<'m, T: 'm> NonZeroElements<'m> for CooMatrix<T> {
    type Elem = T;
    type Iter = NonZeroIter<'m, T>;
    fn nnz(&'m self) -> usize {
        self.data.len()
    }
    fn non_zero_elements(&'m self) -> Self::Iter {
        NonZeroIter {
            pos: 0..self.data.len(),
            matrix: self,
        }
    }
}

// --- NonZeroIterMut ----------------------------------------------------------

pub struct NonZeroIterMut<'m, T> {
    matrix: &'m mut CooMatrix<T>,
    pos: Range<usize>,
}

impl<'mx, T> Iterator for NonZeroIterMut<'mx, T> {
    type Item = (usize, usize, &'mx mut T);
    fn next(&mut self) -> Option<Self::Item> {
        let pos = self.pos.next()?;
        Some((self.matrix.i[pos], self.matrix.j[pos], unsafe {
            std::mem::transmute(&mut self.matrix.data[pos])
        }))
    }
}

impl<'mx, T> ExactSizeIterator for NonZeroIterMut<'mx, T> {
    fn len(&self) -> usize {
        self.pos.len()
    }
}

impl<'mx, T> FusedIterator for NonZeroIterMut<'mx, T> {}

impl<'m, T: 'm> NonZeroElementsMut<'m> for CooMatrix<T> {
    type IterMut = NonZeroIterMut<'m, T>;
    fn non_zero_elements_mut(&'m mut self) -> Self::IterMut {
        NonZeroIterMut {
            pos: 0..self.data.len(),
            matrix: self,
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn to_csr() {
        let mut coo_matrix = CooMatrix::<u8>::new(4, 6);
        coo_matrix.insert(0, 0, 10);
        coo_matrix.insert(0, 1, 20);
        coo_matrix.insert(1, 1, 30);
        coo_matrix.insert(1, 3, 40);
        coo_matrix.insert(2, 2, 50);
        coo_matrix.insert(2, 3, 60);
        coo_matrix.insert(2, 4, 70);
        coo_matrix.insert(3, 5, 80);

        let csr = coo_matrix.to_csr();
        println!("{:?}", csr);
        assert_eq!(csr.data, vec![10, 20, 30, 40, 50, 60, 70, 80]);
        assert_eq!(csr.col_index, vec![0, 1, 1, 3, 2, 3, 4, 5]);
        assert_eq!(csr.row_index, vec![0, 2, 4, 7, 8]);
    }

    #[test]
    fn coo_coo_add() {
        let mut coo_matrix = CooMatrix::<u8>::new(4, 6);
        coo_matrix.insert(0, 0, 10);
        coo_matrix.insert(0, 1, 20);
        coo_matrix.insert(1, 1, 30);
        coo_matrix.insert(1, 3, 40);
        coo_matrix.insert(2, 2, 50);
        coo_matrix.insert(2, 3, 60);
        coo_matrix.insert(2, 4, 70);
        coo_matrix.insert(3, 5, 80);

        let c2 = coo_matrix.clone() + &coo_matrix;
        let mut nz = c2.non_zero_elements();

        assert_eq!(nz.next(), Some((0, 0, &20)));
        assert_eq!(nz.next(), Some((0, 1, &40)));
        assert_eq!(nz.next(), Some((1, 1, &60)));
        assert_eq!(nz.next(), Some((1, 3, &80)));
        assert_eq!(nz.next(), Some((2, 2, &100)));
        assert_eq!(nz.next(), Some((2, 3, &120)));
        assert_eq!(nz.next(), Some((2, 4, &140)));
        assert_eq!(nz.next(), Some((3, 5, &160)));
    }
}
