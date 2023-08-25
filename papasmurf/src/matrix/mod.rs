//! Various formats of sparse and dense matrices.

mod coo;
mod csr;
mod dense;
mod dok;

pub use self::coo::CooMatrix;
pub use self::csr::CsrMatrix;
pub use self::dense::DenseMatrix;
pub use self::dok::DokMatrix;

/// Used to get the dimensions of a matrix.
pub trait MatrixDimensions {
    /// Get the number of rows of the matrix.
    fn rows(&self) -> usize;
    /// Get the number of columns of the matrix.
    fn columns(&self) -> usize;
    /// Get the shape of the matrix.
    fn shape(&self) -> (usize, usize) {
        (self.rows(), self.columns())
    }
}

/// Used to count and iterate over the non-zero elements of a matrix.
impl<M: MatrixDimensions> MatrixDimensions for &M {
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

/// Used to count and iterate over the non-zero elements of a matrix.
pub trait NonZeroElements<'m, T: 'm> {
    /// An iterator over the the non-zero elements of the matrix.
    type Iter: Iterator<Item = (usize, usize, &'m T)> + ExactSizeIterator;
    /// Get an iterator over the non-zero elements of the matrix.
    fn non_zero_elements(&'m self) -> Self::Iter;
    /// Get the number of non-zero elements of the matrix.
    fn nnz(&'m self) -> usize {
        self.non_zero_elements().len()
    }
}

/// The dot-product operator for matrices.
pub trait Dot<Rhs = Self> {
    type Output;
    /// Compute the dot-product between this matrix and another.
    fn dot(self, rhs: Rhs) -> Self::Output;
}
