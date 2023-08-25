use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;

/// A pair of values for paired-end reads.
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Paired<T> {
    /// The value for the forward read.
    pub forward: T,
    /// The value for the backward read.
    pub backward: T,
}

impl<T> Paired<T> {
    /// Create a new pair of values from the given values.
    #[inline]
    pub fn new(forward: T, backward: T) -> Self {
        Self { forward, backward }
    }

    /// Get a reference over the pair of values.
    #[inline]
    pub fn as_ref(&self) -> Paired<&T> {
        Paired::new(&self.forward, &self.backward)
    }

    /// Get a mutable reference over the pair of values.
    #[inline]
    pub fn as_mut(&mut self) -> Paired<&mut T> {
        Paired::new(&mut self.forward, &mut self.backward)
    }

    /// Apply the same function to the forward and the backward values.
    #[inline]
    pub fn map<U, F>(self, f: F) -> Paired<U>
    where
        F: FnOnce(T) -> U + Copy,
    {
        Paired::new(f(self.forward), f(self.backward))
    }

    /// Apply a function to both the forward and the backward values.
    #[inline]
    pub fn merge<U, F>(self, f: F) -> U
    where
        F: FnOnce(T, T) -> U,
    {
        f(self.forward, self.backward)
    }
}

impl<T, E: Debug> Paired<Result<T, E>> {
    /// Unwrap a pair of results.
    #[inline]
    pub fn unwrap(self) -> Paired<T> {
        self.map(Result::unwrap)
    }

    /// Transpose a pair of results into a result of pair.
    pub fn transpose(self) -> Result<Paired<T>, E> {
        Ok(Paired::new(self.forward?, self.backward?))
    }
}

impl<T> Paired<Option<T>> {
    /// Unwrap a pair of options.
    #[inline]
    pub fn unwrap(self) -> Paired<T> {
        self.map(Option::unwrap)
    }

    /// Transpose a pair of options into an option of pair.
    pub fn transpose(self) -> Option<Paired<T>> {
        Some(Paired::new(self.forward?, self.backward?))
    }
}

impl<T> From<(T, T)> for Paired<T> {
    fn from(t: (T, T)) -> Self {
        Self::new(t.0, t.1)
    }
}

impl<T> From<Paired<T>> for (T, T) {
    fn from(pair: Paired<T>) -> Self {
        (pair.forward, pair.backward)
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

#[cfg(feature = "serde")]
mod ser {

    use super::*;
    use serde::ser::SerializeTuple;
    use serde::ser::Serializer;
    use serde::Serialize;

    impl<T: Serialize> Serialize for Paired<T> {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            let mut tup = serializer.serialize_tuple(2)?;
            tup.serialize_element(&self.forward)?;
            tup.serialize_element(&self.backward)?;
            tup.end()
        }
    }
}

#[cfg(feature = "serde")]
mod de {

    use super::*;

    use std::fmt::Formatter;
    use std::fmt::Result as FmtResult;

    use serde::de::Deserializer;
    use serde::de::Error;
    use serde::de::SeqAccess;
    use serde::de::Visitor;
    use serde::Deserialize;

    struct PairedVisitor<T> {
        _marker: std::marker::PhantomData<T>,
    }

    impl<T> Default for PairedVisitor<T> {
        fn default() -> Self {
            Self {
                _marker: std::marker::PhantomData,
            }
        }
    }

    impl<'de, T> Visitor<'de> for PairedVisitor<T>
    where
        T: Deserialize<'de>,
    {
        type Value = Paired<T>;

        fn expecting(&self, formatter: &mut Formatter) -> FmtResult {
            write!(formatter, "a tuple of size 2")
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let forward = seq
                .next_element()?
                .ok_or_else(|| Error::invalid_length(0, &"a tuple of size 2"))?;
            let backward = seq
                .next_element()?
                .ok_or_else(|| Error::invalid_length(1, &"a tuple of size 2"))?;
            Ok(Paired::new(forward, backward))
        }
    }

    impl<'de, T> Deserialize<'de> for Paired<T>
    where
        T: Deserialize<'de>,
    {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            Ok(deserializer.deserialize_seq(PairedVisitor::default())?)
        }
    }
}
