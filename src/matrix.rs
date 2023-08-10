use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;

use lightmotif::num::Unsigned;
use lightmotif::num::U32;

// --- MatrixDimensions --------------------------------------------------------

pub trait MatrixDimensions {
    fn rows(&self) -> usize;
    fn columns(&self) -> usize;
    fn shape(&self) -> (usize, usize) {
        (self.rows(), self.columns())
    }
}

impl<T: MatrixDimensions> MatrixDimensions for &T {
    fn rows(&self) -> usize {
        (*self).rows()
    }
    fn columns(&self) -> usize {
        (*self).columns()
    }
    fn shape(&self) -> (usize, usize) {
        (*self).shape()
    }
}

// --- Matrix ------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct Matrix<T: Default + Copy, A: Unsigned = U32> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
    offset: usize,
    _alignment: std::marker::PhantomData<A>,
}

impl<T: Default + Copy, A: Unsigned> Matrix<T, A> {
    pub fn new(rows: usize, cols: usize) -> Self {
        // Always over-allocate columns to avoid alignment issues.
        let c = cols + (A::USIZE - cols % A::USIZE) * (cols % A::USIZE > 0) as usize;
        let data = vec![T::default(); (rows + 1) * c];

        // compute offset to aligned memory
        let mut offset = 0;
        while data[offset..].as_ptr() as usize % c > 0 {
            offset += 1
        }

        Self {
            offset,
            rows,
            cols,
            data,
            _alignment: std::marker::PhantomData,
        }
    }

    pub fn transpose(&self) -> Self {
        let mut t = Self::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t[j][i] = self[i][j];
            }
        }
        t
    }

    #[inline]
    pub fn stride(&self) -> usize {
        let x = std::mem::size_of::<T>();
        let c = self.cols * x;
        let b =
            self.cols + (A::USIZE - self.cols % A::USIZE) * ((self.cols % A::USIZE) > 0) as usize;
        b / x + ((b % x) > 0) as usize
    }
}

impl<T: Default + Copy, A: Unsigned> Index<usize> for Matrix<T, A> {
    type Output = [T];
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        let c = self.stride();
        let row = self.offset + c * index;
        &self.data[row..row + self.cols]
    }
}

impl<T: Default + Copy, A: Unsigned> IndexMut<usize> for Matrix<T, A> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let c = self.stride();
        let row = self.offset + c * index;
        &mut self.data[row..row + self.cols]
    }
}

impl<T: Default + Copy, A: Unsigned> MatrixDimensions for Matrix<T, A> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn columns(&self) -> usize {
        self.cols
    }
}

// --- DokMatrix ---------------------------------------------------------------

#[derive(Debug, Clone)]
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

    pub fn insert(&mut self, i: usize, j: usize, data: T) {
        assert!(i < self.rows);
        assert!(j < self.cols);
        self.data.insert((i, j), data);
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
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

    pub fn to_csc(&self) -> CscMatrix<T> {
        let mut indices = self.data.keys().map(|(i, j)| (j, i)).collect::<Vec<_>>();
        indices.sort_unstable();

        let mut csc = CscMatrix::new(self.rows, self.cols);
        let mut it = indices.into_iter().peekable();

        for j in 0..self.cols {
            csc.col_index[j] = csc.row_index.len();
            while let Some((y, _)) = it.peek() {
                if **y != j {
                    break;
                }
                let (y, x) = it.next().unwrap();
                csc.row_index.push(*x);
                csc.data.push(self.data.get(&(*x, *y)).unwrap().clone());
            }
        }

        csc.col_index[self.cols] = csc.row_index.len();
        csc
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

// --- CooMatrix ---------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct CooMatrix<T> {
    rows: usize,
    cols: usize,
    i: Vec<usize>,
    j: Vec<usize>,
    data: Vec<T>,
}

impl<T> CooMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            i: Vec::new(),
            j: Vec::new(),
            data: Vec::new()
        }
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    pub fn insert(&mut self, i: usize, j: usize, data: T) {
        assert!(i < self.rows);
        assert!(j < self.cols);

        if self.nnz() > 0 {
            assert!( (&i, &j) > (self.i.last().unwrap(), self.j.last().unwrap()) );
        }

        self.i.push(i);
        self.j.push(j);
        self.data.push(data);
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, &T)> {
        self.i.iter().zip(self.j.iter()).zip(self.data.iter())
            .map(|((i, j), x)| (*i, *j, x) )
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

// --- CscMatrix ---------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct CscMatrix<T> {
    rows: usize,
    data: Vec<T>,
    col_index: Vec<usize>,
    row_index: Vec<usize>,
}

impl<T> CscMatrix<T> {
    fn new(rows: usize, cols: usize) -> Self {
        CscMatrix {
            rows,
            data: Vec::new(),
            row_index: Vec::new(),
            col_index: vec![0; cols + 1],
        }
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
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

// --- CsrMatrix ---------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct CsrMatrix<T> {
    cols: usize,
    data: Vec<T>,
    col_index: Vec<usize>,
    row_index: Vec<usize>,
}

impl<T> CsrMatrix<T> {
    fn new(rows: usize, cols: usize) -> Self {
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

    use super::*;

    fn test_coo_to_csr() {
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
    fn test_dok_to_csr() {
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

    #[test]
    fn test_dok_to_csc() {
        let mut dok_matrix = DokMatrix::<u8>::new(6, 4);
        dok_matrix.insert(0, 0, 10);
        dok_matrix.insert(1, 0, 20);
        dok_matrix.insert(1, 1, 30);
        dok_matrix.insert(3, 1, 40);
        dok_matrix.insert(2, 2, 50);
        dok_matrix.insert(3, 2, 60);
        dok_matrix.insert(4, 2, 70);
        dok_matrix.insert(5, 3, 80);

        let csc = dok_matrix.to_csc();
        assert_eq!(csc.data, vec![10, 20, 30, 40, 50, 60, 70, 80]);
        assert_eq!(csc.row_index, vec![0, 1, 1, 3, 2, 3, 4, 5]);
        assert_eq!(csc.col_index, vec![0, 2, 4, 7, 8]);
    }

    #[test]
    fn test_csr_dot() {
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
