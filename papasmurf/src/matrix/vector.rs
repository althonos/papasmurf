use std::ops::Index;
use std::ops::IndexMut;

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
