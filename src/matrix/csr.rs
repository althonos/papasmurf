use std::cmp::Ordering;
use std::ops::Add;
use std::ops::Mul;

use serde::Deserialize;
use serde::Serialize;

use super::coo::CooMatrix;
use super::csc::CscMatrix;
use super::MatrixDimensions;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CsrMatrix<T> {
    pub(super) cols: usize,
    pub(super) data: Vec<T>,
    pub(super) col_index: Vec<usize>,
    pub(super) row_index: Vec<usize>,
}

impl<T> CsrMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        CsrMatrix {
            cols,
            data: Vec::new(),
            col_index: Vec::new(),
            row_index: vec![0; rows + 1],
        }
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
    }
}

impl<T: Add<Output = T> + Mul<Output = T> + PartialEq + Clone> CsrMatrix<T> {
    pub fn dot(&self, rhs: &CscMatrix<T>) -> CooMatrix<T> {
        assert_eq!(self.cols, rhs.rows);

        let mut out = CooMatrix::new(self.rows(), rhs.columns());

        for i in 0..self.rows() {
            if self.row_index[i] == self.row_index[i + 1] {
                continue;
            }

            let row_cols = &self.col_index[self.row_index[i]..self.row_index[i + 1]];
            let row_data = &self.data[self.row_index[i]..self.row_index[i + 1]];

            for j in 0..rhs.columns() {
                if rhs.col_index[j] == rhs.col_index[j + 1] {
                    continue;
                }

                let col_rows = &rhs.row_index[rhs.col_index[j]..rhs.col_index[j + 1]];
                let col_data = &rhs.data[rhs.col_index[j]..rhs.col_index[j + 1]];

                let mut x: Option<T> = None;
                let mut k1 = 0;
                let mut k2 = 0;

                while k1 < row_cols.len() && k2 < col_rows.len() {
                    match row_cols[k1].cmp(&col_rows[k2]) {
                        Ordering::Less => k1 += 1,
                        Ordering::Greater => k2 += 1,
                        Ordering::Equal => {
                            let p = row_data[k1].clone() * col_data[k2].clone();
                            if let Some(n) = x.as_mut() {
                                *n = n.clone() + p;
                            } else {
                                x = Some(p);
                            }
                            k1 += 1;
                            k2 += 1;
                        }
                    }
                }

                if let Some(res) = x {
                    out.i.push(i);
                    out.j.push(j);
                    out.data.push(res);
                }
            }
        }

        out
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

#[cfg(test)]
mod test {

    use super::super::dok::DokMatrix;

    #[test]
    fn csr_csc_dot() {
        let mut a = DokMatrix::<u8>::new(2, 2);
        a.insert(0, 0, 1);
        a.insert(0, 1, 2);
        a.insert(1, 0, 3);

        let c = a.to_csr().dot(&a.to_csc());
        let mut it = c.iter();
        assert_eq!(it.next(), Some((0, 0, &7)));
        assert_eq!(it.next(), Some((0, 1, &2)));
        assert_eq!(it.next(), Some((1, 0, &3)));
        assert_eq!(it.next(), Some((1, 1, &6)));
        assert_eq!(it.next(), None);
    }
}
