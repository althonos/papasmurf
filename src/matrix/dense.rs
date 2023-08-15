use std::ops::Index;
use std::ops::IndexMut;

use serde::Deserialize;
use serde::Serialize;
use typenum::Unsigned;
use typenum::U32;

use super::MatrixDimensions;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DenseMatrix<T: Default + Copy, A: Unsigned = U32> {
    pub(super) data: Vec<T>,
    pub(super) rows: usize,
    pub(super) cols: usize,
    pub(super) offset: usize,
    pub(super) _alignment: std::marker::PhantomData<A>,
}

impl<T: Default + Copy, A: Unsigned> DenseMatrix<T, A> {
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
        let _c = self.cols * x;
        let b =
            self.cols + (A::USIZE - self.cols % A::USIZE) * ((self.cols % A::USIZE) > 0) as usize;
        b / x + ((b % x) > 0) as usize
    }
}

impl<T: Default + Copy, A: Unsigned> Index<usize> for DenseMatrix<T, A> {
    type Output = [T];
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        let c = self.stride();
        let row = self.offset + c * index;
        &self.data[row..row + self.cols]
    }
}

impl<T: Default + Copy, A: Unsigned> IndexMut<usize> for DenseMatrix<T, A> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        let c = self.stride();
        let row = self.offset + c * index;
        &mut self.data[row..row + self.cols]
    }
}

impl<T: Default + Copy, A: Unsigned> MatrixDimensions for DenseMatrix<T, A> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn columns(&self) -> usize {
        self.cols
    }
}
