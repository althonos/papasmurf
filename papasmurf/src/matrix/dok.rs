use std::collections::HashMap;
use std::iter::FusedIterator;
use std::ops::AddAssign;

use serde::Deserialize;
use serde::Serialize;

use super::csr::CsrMatrix;
use super::MatrixDimensions;
use super::NonZeroElements;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DokMatrix<T> {
    data: HashMap<(usize, usize), T>,
    rows: usize,
    cols: usize,
}

impl<T> DokMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            data: Default::default(),
            rows,
            cols,
        }
    }

    pub fn with_data(rows: usize, cols: usize, data: HashMap<(usize, usize), T>) -> Self {
        assert!( data.keys().map(|(i, j)| i).max().map(|&x| x < rows).unwrap_or(true));
        assert!( data.keys().map(|(i, j)| j).max().map(|&x| x < cols).unwrap_or(true));
        Self {
            rows,
            cols,
            data,
        }
    }

    pub fn insert(&mut self, i: usize, j: usize, data: T) {
        assert!(i < self.rows);
        assert!(j < self.cols);
        self.data.insert((i, j), data);
    }

    pub fn grow(&mut self, rows: usize, cols: usize) {
        self.rows += rows;
        self.cols += cols;
    }
}

impl<T: Default + Clone> DokMatrix<T> {
    pub fn get(&self, i: usize, j: usize) -> T {
        self.data.get(&(i, j)).cloned().unwrap_or_default()
    }
}

impl<T: Clone> DokMatrix<T> {
    pub fn to_csr(&self) -> CsrMatrix<T> {
        let mut indices = self.data.keys().collect::<Vec<_>>();
        indices.sort_unstable();

        let mut csr = CsrMatrix::new(self.rows, self.cols);
        let mut it = indices.into_iter().peekable();

        for i in 0..self.rows {
            csr.row_index[i] = csr.col_index.len();
            while let Some((x, _)) = it.peek() {
                if *x != i {
                    break;
                }
                let (x, y) = it.next().unwrap();
                csr.col_index.push(*y);
                csr.data.push(self.data.get(&(*x, *y)).unwrap().clone());
            }
        }

        csr.row_index[self.rows] = csr.col_index.len();
        csr
    }
}

impl<T> Default for DokMatrix<T> {
    fn default() -> Self {
        Self {
            data: HashMap::new(),
            rows: 0,
            cols: 0,
        }
    }
}

impl<T: AddAssign + Clone> AddAssign<Self> for DokMatrix<T> {
    fn add_assign(&mut self, rhs: Self) {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);
        for (coord, rval) in rhs.data.iter() {
            if let Some(lval) = self.data.get_mut(&coord) {
                lval.add_assign(rval.clone());
            } else {
                self.data.insert(*coord, rval.clone());
            }
        }
    }
}

impl<T> AsRef<HashMap<(usize, usize), T>> for DokMatrix<T> {
    fn as_ref(&self) -> &HashMap<(usize, usize), T> {
        &self.data
    }
}

impl<T> MatrixDimensions for DokMatrix<T> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn columns(&self) -> usize {
        self.cols
    }
}

pub struct NonZeroIter<'m, T> {
    it: std::collections::hash_map::Iter<'m, (usize, usize), T>,
}

impl<'mx, T> Iterator for NonZeroIter<'mx, T> {
    type Item = (usize, usize, &'mx T);
    fn next(&mut self) -> Option<Self::Item> {
        self.it.next().map(|((i, j), x)| (*i, *j, x))
    }
}

impl<'mx, T> ExactSizeIterator for NonZeroIter<'mx, T> {
    fn len(&self) -> usize {
        self.it.len()
    }
}

impl<'mx, T> FusedIterator for NonZeroIter<'mx, T> {}

impl<'m, T: 'm> NonZeroElements<'m, T> for DokMatrix<T> {
    type Iter = NonZeroIter<'m, T>;
    fn non_zero_elements(&'m self) -> Self::Iter {
        NonZeroIter {
            it: self.data.iter(),
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn to_csr() {
        let mut dok_matrix = DokMatrix::<u8>::new(4, 6);
        dok_matrix.insert(0, 0, 10);
        dok_matrix.insert(0, 1, 20);
        dok_matrix.insert(1, 1, 30);
        dok_matrix.insert(1, 3, 40);
        dok_matrix.insert(2, 2, 50);
        dok_matrix.insert(2, 3, 60);
        dok_matrix.insert(2, 4, 70);
        dok_matrix.insert(3, 5, 80);

        let csr = dok_matrix.to_csr();
        println!("{:?}", csr);
        assert_eq!(csr.data, vec![10, 20, 30, 40, 50, 60, 70, 80]);
        assert_eq!(csr.col_index, vec![0, 1, 1, 3, 2, 3, 4, 5]);
        assert_eq!(csr.row_index, vec![0, 2, 4, 7, 8]);
    }

    // #[test]
    // fn to_csc() {
    //     let mut dok_matrix = DokMatrix::<u8>::new(6, 4);
    //     dok_matrix.insert(0, 0, 10);
    //     dok_matrix.insert(1, 0, 20);
    //     dok_matrix.insert(1, 1, 30);
    //     dok_matrix.insert(3, 1, 40);
    //     dok_matrix.insert(2, 2, 50);
    //     dok_matrix.insert(3, 2, 60);
    //     dok_matrix.insert(4, 2, 70);
    //     dok_matrix.insert(5, 3, 80);

    //     let csc = dok_matrix.to_csc();
    //     assert_eq!(csc.data, vec![10, 20, 30, 40, 50, 60, 70, 80]);
    //     assert_eq!(csc.row_index, vec![0, 1, 1, 3, 2, 3, 4, 5]);
    //     assert_eq!(csc.col_index, vec![0, 2, 4, 7, 8]);
    // }

    #[test]
    fn non_zero_elements() {
        let mut dok_matrix = DokMatrix::<u8>::new(4, 6);
        dok_matrix.insert(0, 0, 10);
        dok_matrix.insert(1, 3, 40);

        let mut it = dok_matrix.non_zero_elements();
        assert!(matches!(it.next(), Some((1, 3, &40)) | Some((0, 0, &10))));
        assert!(matches!(it.next(), Some((1, 3, &40)) | Some((0, 0, &10))));
        assert_eq!(it.next(), None);
    }
}
