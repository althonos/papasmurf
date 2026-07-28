use std::ops::AddAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;
use std::process::Output;

use crate::matrix::MatrixDimensions;
use crate::matrix::NonZeroElements;

use super::Dot;

#[derive(Debug, Clone)]
pub struct Vector<T> {
    data: Vec<T>,
}

impl<T> Vector<T> {
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl<T: Default> Vector<T> {
    pub fn new(length: usize) -> Self {
        Self {
            data: (0..length).map(|_| T::default()).collect(),
        }
    }
}

impl<T> AsRef<[T]> for Vector<T> {
    fn as_ref(&self) -> &[T] {
        self.data.as_slice()
    }
}

impl<T> AsMut<[T]> for Vector<T> {
    fn as_mut(&mut self) -> &mut [T] {
        self.data.as_mut_slice()
    }
}

impl<T> From<Vec<T>> for Vector<T> {
    fn from(data: Vec<T>) -> Self {
        Self { data }
    }
}

impl<T> Index<usize> for Vector<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T> IndexMut<usize> for Vector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<T, M> Dot<M> for Vector<T>
where
    T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default,
    M: MatrixDimensions + for<'m> NonZeroElements<'m, Elem = T>,
{
    type Output = Vector<T>;
    fn dot(self, mat: M) -> Self::Output {
        assert_eq!(self.len(), mat.rows());
        let mut out = Vector::new(mat.columns());
        for (i, j, x) in mat.non_zero_elements() {
            out[j] += self[i].clone() * x.clone();
        }
        out
    }
}

impl<T, M> Dot<Vector<T>> for M
where
    T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default,
    M: MatrixDimensions + for<'m> NonZeroElements<'m, Elem = T>,
{
    type Output = Vector<T>;
    fn dot(self, v: Vector<T>) -> Self::Output {
        assert_eq!(self.columns(), v.len());
        let mut out = Vector::new(self.rows());
        for (i, j, x) in self.non_zero_elements() {
            out[i] += x.clone() * v[j].clone();
        }
        out
    }
}
