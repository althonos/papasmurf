use std::cmp::Ordering;
use std::ops::Add;

use super::csr::CsrMatrix;
use super::MatrixDimensions;

#[derive(Debug, Clone)]
pub struct CooMatrix<T> {
    pub(super) rows: usize,
    pub(super) cols: usize,
    pub(super) i: Vec<usize>,
    pub(super) j: Vec<usize>,
    pub(super) data: Vec<T>,
}

impl<T> CooMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            i: Vec::new(),
            j: Vec::new(),
            data: Vec::new(),
        }
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    pub fn insert(&mut self, i: usize, j: usize, data: T) {
        assert!(i < self.rows);
        assert!(j < self.cols);

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

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &T)> {
        self.i
            .iter()
            .zip(self.j.iter())
            .zip(self.data.iter())
            .map(|((i, j), x)| (*i, *j, x))
    }

    pub fn grow(&mut self, rows: usize, cols: usize) {
        self.rows += rows;
        self.cols += cols;
    }
}

impl<T: Clone> CooMatrix<T> {
    pub fn to_csr(&self) -> CsrMatrix<T> {
        let mut csr = CsrMatrix::new(self.rows, self.cols);

        // let mut indices = self.data.keys().collect::<Vec<_>>();
        // indices.sort_unstable();

        let mut it = self.iter().peekable();

        for i in 0..self.rows {
            csr.row_index[i] = csr.col_index.len();
            while let Some((x, _, _)) = it.peek() {
                if *x != i {
                    break;
                }
                let (x, y, z) = it.next().unwrap();
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

impl<T: Add<Output = T> + Clone> Add for CooMatrix<T> {
    type Output = CooMatrix<T>;
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);

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
                    out.insert(i2, j2, d2.clone());
                    y += 1;
                }
                Ordering::Equal => {
                    // FIXME: May be zero
                    out.insert(i1, j1, d1.clone() + d2.clone());
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
            out.data.push(rhs.data[y].clone());
            y += 1;
        }

        out
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
}
