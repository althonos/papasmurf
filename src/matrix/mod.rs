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

pub trait NonZeroElements<'m, T: 'm> {
    type Iter: Iterator<Item = (usize, usize, &'m T)> + ExactSizeIterator;
    fn non_zero_elements(&'m self) -> Self::Iter;
    fn nnz(&'m self) -> usize {
        self.non_zero_elements().len()
    }
}

pub trait Dot<Rhs = Self> {
    type Output;
    fn dot(self, rhs: Rhs) -> Self::Output;
}

// impl<'mx, T, M: NonZeroElements<T>> NonZeroElements<'mx, T> for &M {
//     type IntoIter = <M as NonZeroElements<T>>::IntoIter;
//     fn non_zero_elements(&'mx self) -> Self::IntoIter {
//         (*self).non_zero_elements()
//     }
// }
