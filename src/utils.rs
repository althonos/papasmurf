use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;
use std::ops::Index;
use std::sync::Arc;
use std::sync::RwLock;

/// The reference counter type used by the entire crate.
pub type Rc<T> = Arc<T>;

/// A generic interner for borrowed data.
#[derive(Debug)]
pub struct Interner<T>
where
    T: ?Sized,
{
    data: Rc<RwLock<HashMap<Rc<T>, Rc<T>>>>,
}

impl<T: ?Sized> Clone for Interner<T> {
    fn clone(&self) -> Self {
        Self {
            data: self.data.clone(),
        }
    }
}

impl<T: ?Sized> Default for Interner<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: ?Sized> Interner<T> {
    pub fn new() -> Self {
        Self {
            data: Default::default(),
        }
    }

    /// Intern a DNA sequence or get a reference to a previously interned string.
    pub fn intern(&self, obj: &T) -> Rc<T>
    where
        T: Hash + Eq,
        for<'a> &'a T: Into<Rc<T>>,
    {
        {
            let r = self.data.read().expect("failed to acquired lock");
            if let Some(rc) = r.get(obj) {
                return rc.clone();
            }
        }

        let mut w = self.data.write().expect("failed to acquired lock");
        let rc: Rc<T> = obj.into();
        w.insert(rc.clone(), rc.clone());
        rc
    }
}

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
