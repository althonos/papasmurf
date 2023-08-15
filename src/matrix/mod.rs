mod coo;
mod csc;
mod csr;
mod dense;
mod dok;

pub use self::coo::CooMatrix;
pub use self::csc::CscMatrix;
pub use self::csr::CsrMatrix;
pub use self::dense::DenseMatrix;
pub use self::dok::DokMatrix;

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
