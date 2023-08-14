use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;
use std::ops::Index;
use std::sync::Arc;
use std::sync::RwLock;

/// A pair of values for paired-end reads.
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Paired<T> {
    pub forward: T,
    pub backward: T,
}

impl<T> Paired<T> {
    #[inline]
    pub fn new(forward: T, backward: T) -> Self {
        Self { forward, backward }
    }

    #[inline]
    pub fn as_ref(&self) -> Paired<&T> {
        Paired::new(&self.forward, &self.backward)
    }

    #[inline]
    pub fn as_mut(&mut self) -> Paired<&mut T> {
        Paired::new(&mut self.forward, &mut self.backward)
    }

    #[inline]
    pub fn map<U, F>(self, f: F) -> Paired<U>
    where
        F: FnOnce(T) -> U + Copy,
    {
        Paired::new(f(self.forward), f(self.backward))
    }
}

impl<T> From<(T, T)> for Paired<T> {
    fn from(t: (T, T)) -> Self {
        Self::new(t.0, t.1)
    }
}

impl<T> FromIterator<(T, T)> for Paired<HashSet<T>>
where
    T: PartialEq + Eq + Hash,
{
    fn from_iter<I: IntoIterator<Item = (T, T)>>(it: I) -> Self {
        let mut p = Paired::<HashSet<T>>::default();
        for (x, y) in it {
            p.forward.insert(x);
            p.backward.insert(y);
        }
        p
    }
}

impl<T> FromIterator<Paired<T>> for Paired<HashSet<T>>
where
    T: PartialEq + Eq + Hash,
{
    fn from_iter<I: IntoIterator<Item = Paired<T>>>(it: I) -> Self {
        Self::from_iter(it.into_iter().map(|p| (p.forward, p.backward)))
    }
}
