//! Various formats of sparse and dense matrices.

mod coo;
mod csc;
mod csr;
mod dense;
mod dok;
mod vector;

use std::ops::AddAssign;

pub use self::coo::CooMatrix;
pub use self::csc::CscMatrix;
pub use self::csr::CsrMatrix;
pub use self::dense::DenseMatrix;
pub use self::dok::DokMatrix;
pub use self::vector::Vector;

// --- Matrix dimensions--------------------------------------------------------

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

#[derive(Debug, Clone, Copy)]
pub struct Rows;

#[derive(Debug, Clone, Copy)]
pub struct Columns;

pub trait Dimension {
    fn as_usize() -> usize;
}

impl Dimension for Rows {
    fn as_usize() -> usize {
        0
    }
}

impl Dimension for Columns {
    fn as_usize() -> usize {
        1
    }
}

// --- Non zero iterator -------------------------------------------------------

/// Used to count and iterate over the non-zero elements of a matrix.
pub trait NonZeroElements<'m> {
    /// The element type.
    type Elem: 'm;
    /// An iterator over the the non-zero elements of the matrix.
    type Iter: Iterator<Item = (usize, usize, &'m Self::Elem)> + ExactSizeIterator;
    /// Get an iterator over the non-zero elements of the matrix.
    fn non_zero_elements(&'m self) -> Self::Iter;
    /// Get the number of non-zero elements of the matrix.
    fn nnz(&'m self) -> usize {
        self.non_zero_elements().len()
    }
}

pub trait NonZeroElementsMut<'m>: NonZeroElements<'m> {
    type IterMut: Iterator<Item = (usize, usize, &'m mut Self::Elem)> + ExactSizeIterator;
    fn non_zero_elements_mut(&'m mut self) -> Self::IterMut;
}

// --- Matrix operations -------------------------------------------------------

/// The dot-product operator for matrices.
pub trait Dot<Rhs = Self> {
    type Output;
    /// Compute the dot-product between this matrix and another.
    fn dot(self, rhs: Rhs) -> Self::Output;
}

/// Vertical concatenation for matrices.
pub trait VerticalStack<Rhs = Self> {
    fn vstack(&mut self, other: Rhs);
}

/// Sum along a given axis.
pub trait AxisSum<D: Dimension> {
    type Output;
    fn sum_axis(&self, dim: D) -> Self::Output;
}

macro_rules! impl_sum {
    ($MAT:ident) => {
        impl<T> AxisSum<Rows> for $MAT<T>
        where
            T: Default + Clone + AddAssign,
        {
            type Output = Vector<T>;
            fn sum_axis(&self, _dim: Rows) -> Self::Output {
                let mut out = Vector::new(self.columns());
                for (_, j, x) in self.non_zero_elements() {
                    out[j] += x.clone();
                }
                out
            }
        }

        impl<T> AxisSum<Columns> for $MAT<T>
        where
            T: Default + Clone + AddAssign,
        {
            type Output = Vector<T>;
            fn sum_axis(&self, _dim: Columns) -> Self::Output {
                let mut out = Vector::new(self.rows());
                for (i, _, x) in self.non_zero_elements() {
                    out[i] += x.clone();
                }
                out
            }
        }
    };
}

impl_sum!(CsrMatrix);
impl_sum!(CscMatrix);
impl_sum!(CooMatrix);
impl_sum!(DokMatrix);

// --- Unique ------------------------------------------------------------------

/// Indices for unique rows / columns of a matrix.
#[derive(Debug, Clone)]
pub struct Unique {
    pub indices: Vec<usize>,
    pub reverse: Vec<usize>,
}
