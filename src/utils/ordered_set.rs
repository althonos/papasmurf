use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;
use std::ops::Index;
use std::sync::Arc;
use std::sync::RwLock;

/// An immutable ordered set, mapping unique elements to a fixed index.
#[derive(Debug, Clone)]
pub struct OrderedSet<T> {
    data: Vec<T>,
    indices: HashMap<T, usize>,
}

impl<T> OrderedSet<T> {
    pub fn len(&self) -> usize {
        self.data.len()
    }
    pub fn iter<'a>(&'a self) -> <&'a Self as IntoIterator>::IntoIter {
        (&self).into_iter()
    }
}

impl<T> FromIterator<T> for OrderedSet<T>
where
    T: Eq + Hash + Clone + Ord,
{
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut dedup = iter
            .into_iter()
            .enumerate()
            .map(|(i, x)| (x, i))
            .collect::<HashMap<_, _>>();
        let mut keys = dedup.keys().cloned().collect::<Vec<_>>();
        keys.sort_unstable();
        for (i, k) in keys.iter().cloned().enumerate() {
            dedup.insert(k, i);
        }

        Self {
            data: keys,
            indices: dedup,
        }
    }
}

impl<T> From<HashSet<T>> for OrderedSet<T>
where
    T: Eq + Hash + Clone + Ord,
{
    fn from(set: HashSet<T>) -> Self {
        let mut data = set.into_iter().collect::<Vec<_>>();
        data.sort_unstable();
        let mut indices = HashMap::with_capacity(data.len());
        for (i, k) in data.iter().cloned().enumerate() {
            indices.insert(k, i);
        }

        Self { data, indices }
    }
}

impl<'a, T> IntoIterator for &'a OrderedSet<T> {
    type Item = &'a T;
    type IntoIter = <&'a Vec<T> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter {
        (&self.data).into_iter()
    }
}

impl<T> Index<usize> for OrderedSet<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T, K> Index<&K> for OrderedSet<T>
where
    T: Borrow<K> + Eq + Hash,
    K: Hash + Eq + ?Sized,
{
    type Output = usize;
    fn index(&self, index: &K) -> &Self::Output {
        &self.indices[index]
    }
}
