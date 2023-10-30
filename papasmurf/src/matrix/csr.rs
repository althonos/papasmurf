use std::iter::FusedIterator;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;

use serde::Deserialize;
use serde::Serialize;

use super::coo::CooMatrix;
use super::Dot;
use super::MatrixDimensions;
use super::NonZeroElements;

// --- CsrMatrix ---------------------------------------------------------------

/// A sparse matrix in compressed sparse row (CSR) format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CsrMatrix<T> {
    pub(super) cols: usize,
    pub(super) data: Vec<T>,
    pub(super) col_index: Vec<usize>,
    pub(super) row_index: Vec<usize>,
}

impl<T> CsrMatrix<T> {
    /// Create a new CSR matrix with the given dimensions.
    pub fn new(rows: usize, cols: usize) -> Self {
        CsrMatrix {
            cols,
            data: Vec::new(),
            col_index: Vec::new(),
            row_index: vec![0; rows + 1],
        }
    }

    /// Reserve space for the given number of non-zero elements.
    pub fn reserve(&mut self, nnz: usize) {
        self.data.reserve(nnz);
        self.col_index.reserve(nnz);
    }

    /// Convert the matrix into COO format without cloning data.
    pub fn into_coo(self) -> CooMatrix<T> {
        let mut coo = CooMatrix::new(self.rows(), self.columns());
        coo.reserve(self.nnz());

        let mut ptr = 0;
        let mut row = 0;
        let mut it = self.data.into_iter();

        while let Some(x) = it.next() {
            while ptr >= self.row_index[row + 1] {
                row += 1;
            }
            ptr += 1;
            coo.i.push(row);
            coo.j.push(self.col_index[ptr - 1]);
            coo.data.push(x);
        }

        coo
    }
}

impl<T: Clone> CsrMatrix<T> {
    /// Build a COO matrix by cloning data.
    pub fn to_coo(&self) -> CooMatrix<T> {
        let mut coo = CooMatrix::new(self.rows(), self.columns());
        for (i, j, x) in self.non_zero_elements() {
            coo.insert(i, j, x.clone());
        }
        coo
    }
}

impl<T: Add<Output = T> + PartialEq + Clone + Default> Add<&CsrMatrix<T>> for CsrMatrix<T> {
    type Output = CooMatrix<T>;
    fn add(self, rhs: &CsrMatrix<T>) -> Self::Output {
        assert_eq!(self.rows(), rhs.rows());
        assert_eq!(self.columns(), rhs.columns());

        let mut out = CooMatrix::new(self.rows(), self.columns());

        let mut l = self.non_zero_elements();
        let mut r = self.non_zero_elements();
        let mut l_item = l.next();
        let mut r_item = r.next();

        loop {
            match (l_item, r_item) {
                (Some((li, lj, lx)), Some((ri, rj, _))) if (li, lj) < (ri, rj) => {
                    out.insert(li, lj, lx.clone());
                    l_item = l.next();
                }
                (Some((li, lj, lx)), None) => {
                    out.insert(li, lj, lx.clone());
                    l_item = l.next();
                }
                (Some((li, lj, _)), Some((ri, rj, rx))) if (li, lj) > (ri, rj) => {
                    out.insert(ri, rj, rx.clone());
                    r_item = r.next();
                }
                (None, Some((ri, rj, rx))) => {
                    out.insert(ri, rj, rx.clone());
                    r_item = r.next();
                }
                (Some((_, _, lx)), Some((ri, rj, rx))) => {
                    let x = lx.clone() + rx.clone();
                    if x != T::default() {
                        out.insert(ri, rj, x);
                    }
                    l_item = l.next();
                    r_item = r.next();
                }
                (None, None) => break,
            }
        }

        out
    }
}

impl<T: Add<Output = T> + PartialEq + Clone + Default> Add<CsrMatrix<T>> for CsrMatrix<T> {
    type Output = CooMatrix<T>;
    fn add(self, rhs: CsrMatrix<T>) -> Self::Output {
        self.add(&rhs)
    }
}

/// Implement sparse matrix multiplication using the Gustavson algorithm.
///
/// # References
/// - Gustavson, Fred G. ‘Two Fast Algorithms for Sparse Matrices:
///   Multiplication and Permuted Transposition’. ACM Transactions on
///   Mathematical Software 4, no. 3 (September 1978): 250–69.
///   https://doi.org/10.1145/355791.355796.
/// - Alexandrov, Luben. ‘Parallel Sparse Matrix-Matrix Multiplication’.
///   Hochschulschrift, Institut für Theoretische Informatik (ITI),
///   Karlsruher Institut für Technologie (December 2014).
///   https://doi.org/10.5445/IR/1000128898.
impl<T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default> Dot<&CsrMatrix<T>>
    for CsrMatrix<T>
{
    type Output = CsrMatrix<T>;
    fn dot(self, rhs: &CsrMatrix<T>) -> CsrMatrix<T> {
        assert_eq!(self.columns(), rhs.rows());

        let mut out = CsrMatrix::new(self.rows(), rhs.columns());

        let mut ip = 0;
        let mut x = vec![T::default(); rhs.columns()];
        let mut xb = vec![usize::MAX; rhs.columns()];

        for i in 0..self.row_index.len() - 1 {
            out.row_index[i] = ip;

            let start_row_a = self.row_index[i];
            let end_row_a = self.row_index[i + 1];

            for jp in start_row_a..end_row_a {
                let j = self.col_index[jp];

                let start_row_b = rhs.row_index[j];
                let end_row_b = rhs.row_index[j + 1];

                for kp in start_row_b..end_row_b {
                    let k = rhs.col_index[kp];
                    if xb[k] != i {
                        out.col_index.push(k);
                        ip += 1;
                        xb[k] = i;
                        x[k] = self.data[jp].clone() * rhs.data[kp].clone();
                    } else {
                        x[k] += self.data[jp].clone() * rhs.data[kp].clone();
                    }
                }
            }

            for vp in out.row_index[i]..ip {
                let v = out.col_index[vp];
                out.data.push(x[v].clone());
            }
        }

        let n = out.row_index.len();
        out.row_index[n - 1] = ip;

        let mut tmp = Vec::new();
        for p in 0..out.row_index.len() - 1 {
            let row = out.row_index[p]..out.row_index[p + 1];

            tmp.clear();
            for i in row {
                tmp.push((out.col_index[i], out.data[i].clone()));
            }
            tmp.sort_unstable_by_key(|(j, _)| *j);

            let mut k = out.row_index[p];
            for (j, x) in tmp.drain(..) {
                out.col_index[k] = j;
                out.data[k] = x;
                k += 1;
            }
        }

        out
    }
}

impl<T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default> Dot<CsrMatrix<T>>
    for CsrMatrix<T>
{
    type Output = CsrMatrix<T>;
    fn dot(self, rhs: CsrMatrix<T>) -> CsrMatrix<T> {
        self.dot(&rhs)
    }
}

impl<T> Default for CsrMatrix<T> {
    fn default() -> Self {
        Self {
            cols: 0,
            data: Vec::new(),
            col_index: Vec::new(),
            row_index: Vec::new(),
        }
    }
}

impl<T> MatrixDimensions for CsrMatrix<T> {
    #[inline]
    fn rows(&self) -> usize {
        self.row_index.len() - 1
    }

    #[inline]
    fn columns(&self) -> usize {
        self.cols
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

// --- NonZeroIter -------------------------------------------------------------

pub struct NonZeroIter<'m, T> {
    matrix: &'m CsrMatrix<T>,
    row: usize,
    ptr: usize,
}

impl<'mx, T> Iterator for NonZeroIter<'mx, T> {
    type Item = (usize, usize, &'mx T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr >= self.matrix.data.len() {
            return None;
        }
        while self.ptr >= self.matrix.row_index[self.row + 1] {
            if self.row + 1 > self.matrix.row_index.len() {
                return None;
            }
            self.row += 1;
        }
        self.ptr += 1;
        Some((
            self.row,
            self.matrix.col_index[self.ptr - 1],
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

impl<'m, T: 'm> NonZeroElements<'m, T> for CsrMatrix<T> {
    type Iter = NonZeroIter<'m, T>;
    fn non_zero_elements(&'m self) -> Self::Iter {
        NonZeroIter {
            row: 0,
            ptr: 0,
            matrix: self,
        }
    }
}

#[cfg(test)]
mod test {

    use super::super::dok::DokMatrix;
    use super::super::Dot;
    use super::*;

    #[test]
    fn into_coo() {
        let mut a = DokMatrix::<u8>::new(2, 2);
        a.insert(0, 0, 1);
        a.insert(0, 1, 2);
        a.insert(1, 0, 3);

        let c = a.to_csr().into_coo();
        let mut it = c.non_zero_elements();
        assert_eq!(it.next(), Some((0, 0, &1)));
        assert_eq!(it.next(), Some((0, 1, &2)));
        assert_eq!(it.next(), Some((1, 0, &3)));
        assert_eq!(it.next(), None);
    }

    #[test]
    fn csr_csr_dot() {
        let mut a = DokMatrix::<u8>::new(2, 2);
        a.insert(0, 0, 1);
        a.insert(0, 1, 2);
        a.insert(1, 0, 3);

        let c = a.to_csr().dot(a.to_csr()).to_coo();
        let mut it = c.non_zero_elements();
        assert_eq!(it.next(), Some((0, 0, &7)));
        assert_eq!(it.next(), Some((0, 1, &2)));
        assert_eq!(it.next(), Some((1, 0, &3)));
        assert_eq!(it.next(), Some((1, 1, &6)));
        assert_eq!(it.next(), None);
    }

    #[test]
    fn non_zero_elements() {
        let m1 = CsrMatrix::<u8>::new(2, 2);
        let mut it = m1.non_zero_elements();
        assert_eq!(it.next(), None);

        let mut a = DokMatrix::<u8>::new(2, 2);
        a.insert(0, 0, 1);
        a.insert(0, 1, 2);
        a.insert(1, 0, 3);
        let m2 = a.to_csr();

        let mut it = m2.non_zero_elements();
        assert_eq!(it.next(), Some((0, 0, &1)));
        assert_eq!(it.next(), Some((0, 1, &2)));
        assert_eq!(it.next(), Some((1, 0, &3)));
        assert_eq!(it.next(), None);
    }
}
