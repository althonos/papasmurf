use std::iter::Sum;
use std::ops::AddAssign;
use std::ops::DivAssign;
use std::ops::Index;
use std::ops::IndexMut;
use std::ops::Mul;

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

    pub fn iter(&self) -> <&Vec<T> as IntoIterator>::IntoIter {
        self.data.iter()
    }

    pub fn iter_mut(&mut self) -> <&mut Vec<T> as IntoIterator>::IntoIter {
        self.data.iter_mut()
    }
}

impl<T: Default> Vector<T> {
    pub fn new(length: usize) -> Self {
        Self {
            data: (0..length).map(|_| T::default()).collect(),
        }
    }
}

impl<T> Vector<T>
where
    T: for<'a> Sum<&'a T> + DivAssign + Clone,
{
    pub fn normalize(&mut self) {
        let total = self.data.iter().sum::<T>();
        for x in self.data.iter_mut() {
            *x /= total.clone();
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

impl<T, M> Dot<&M> for Vector<T>
where
    T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default,
    M: MatrixDimensions + for<'m> NonZeroElements<'m, Elem = T>,
{
    type Output = Vector<T>;
    fn dot(self, mat: &M) -> Self::Output {
        (&self).dot(mat)
    }
}

impl<T, M> Dot<&M> for &Vector<T>
where
    T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default,
    M: MatrixDimensions + for<'m> NonZeroElements<'m, Elem = T>,
{
    type Output = Vector<T>;
    fn dot(self, mat: &M) -> Self::Output {
        assert_eq!(self.len(), mat.rows());
        let mut out = Vector::new(mat.columns());
        for (i, j, x) in mat.non_zero_elements() {
            out[j] += self[i].clone() * x.clone();
        }
        out
    }
}

impl<T, M> Dot<Vector<T>> for &M
where
    T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default,
    M: MatrixDimensions + for<'m> NonZeroElements<'m, Elem = T>,
{
    type Output = Vector<T>;
    fn dot(self, v: Vector<T>) -> Self::Output {
        self.dot(&v)
    }
}

impl<T, M> Dot<&Vector<T>> for &M
where
    T: AddAssign + Mul<Output = T> + PartialEq + Clone + Default,
    M: MatrixDimensions + for<'m> NonZeroElements<'m, Elem = T>,
{
    type Output = Vector<T>;
    fn dot(self, v: &Vector<T>) -> Self::Output {
        assert_eq!(self.columns(), v.len());
        let mut out = Vector::new(self.rows());
        for (i, j, x) in self.non_zero_elements() {
            out[i] += x.clone() * v[j].clone();
        }
        out
    }
}

impl<T> IntoIterator for Vector<T> {
    type Item = <Vec<T> as IntoIterator>::Item;
    type IntoIter = <Vec<T> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Vector<T> {
    type Item = <&'a Vec<T> as IntoIterator>::Item;
    type IntoIter = <&'a Vec<T> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}
