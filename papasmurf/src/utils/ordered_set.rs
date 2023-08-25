use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Formatter;
use std::fmt::Result as FmtResult;
use std::hash::Hash;
use std::ops::Index;

/// An immutable ordered set, mapping unique elements to a fixed index.
#[derive(Debug, Clone)]
pub struct OrderedSet<T> {
    data: Vec<T>,
    indices: HashMap<T, usize>,
}

impl<T> OrderedSet<T> {
    /// Return the number of items in the ordered set.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Get an iterator over references to the items stored in the set.
    pub fn iter<'a>(&'a self) -> <&'a Self as IntoIterator>::IntoIter {
        (&self).into_iter()
    }

    /// Get a reference to the memory where the items are stored.
    pub fn as_slice(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl<T> AsRef<[T]> for &OrderedSet<T> {
    fn as_ref(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl<T> FromIterator<T> for OrderedSet<T>
where
    T: Eq + Hash + Clone + Ord,
{
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut indices = iter
            .into_iter()
            .enumerate()
            .map(|(i, x)| (x, i))
            .collect::<HashMap<_, _>>();
        let mut data = indices.keys().cloned().collect::<Vec<_>>();
        data.sort_unstable();
        for (i, k) in data.iter().cloned().enumerate() {
            indices.insert(k, i);
        }

        Self { data, indices }
    }
}

impl<T> From<Vec<T>> for OrderedSet<T>
where
    T: Eq + Hash + Clone + Ord,
{
    fn from(mut data: Vec<T>) -> Self {
        data.sort_unstable();
        data.dedup();
        let mut indices = HashMap::new();
        for (i, k) in data.iter().cloned().enumerate() {
            indices.insert(k, i);
        }
        Self { data, indices }
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

// --- Serde -------------------------------------------------------------------

#[cfg(feature = "serde")]
mod ser {

    use super::*;
    use serde::ser::SerializeSeq;
    use serde::ser::Serializer;
    use serde::Serialize;

    impl<T: Serialize> Serialize for OrderedSet<T> {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            let mut seq = serializer.serialize_seq(Some(self.len()))?;
            for element in self {
                seq.serialize_element(element)?;
            }
            seq.end()
        }
    }
}

#[cfg(feature = "serde")]
mod de {

    use super::*;
    use serde::de::Deserializer;
    use serde::de::SeqAccess;
    use serde::de::Visitor;
    use serde::Deserialize;

    struct OrderedSetVisitor<T> {
        _marker: std::marker::PhantomData<T>,
    }

    impl<T> Default for OrderedSetVisitor<T> {
        fn default() -> Self {
            Self {
                _marker: std::marker::PhantomData,
            }
        }
    }

    impl<'de, T> Visitor<'de> for OrderedSetVisitor<T>
    where
        T: Clone + Eq + Hash + Deserialize<'de>,
    {
        type Value = OrderedSet<T>;

        fn expecting(&self, formatter: &mut Formatter) -> FmtResult {
            write!(formatter, "a sequence of values")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let mut data = Vec::with_capacity(seq.size_hint().unwrap_or(0));
            let mut indices = HashMap::with_capacity(seq.size_hint().unwrap_or(0));

            while let Some(item) = seq.next_element::<T>()? {
                indices.insert(item.clone(), data.len());
                data.push(item);
            }

            Ok(OrderedSet { data, indices })
        }
    }

    impl<'de, T> Deserialize<'de> for OrderedSet<T>
    where
        T: Clone + Eq + Hash + Deserialize<'de>,
    {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            Ok(deserializer.deserialize_seq(OrderedSetVisitor::default())?)
        }
    }
}
