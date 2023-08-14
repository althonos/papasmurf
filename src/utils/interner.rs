use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;
use std::ops::Index;
use std::sync::Arc;
use std::sync::RwLock;

use super::Rc;

/// A generic interner for borrowed data.
#[derive(Debug)]
pub struct Interner<T>
where
    T: ?Sized,
{
    data: Rc<RwLock<HashMap<Arc<T>, Arc<T>>>>,
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

    /// Intern a value or get a reference to the previously interned value.
    pub fn intern(&self, obj: &T) -> Arc<T>
    where
        T: Hash + Eq,
        for<'a> &'a T: Into<Arc<T>>,
    {
        {
            let r = self.data.read().expect("failed to acquired lock");
            if let Some(arc) = r.get(obj) {
                return arc.clone();
            }
        }

        let mut w = self.data.write().expect("failed to acquired lock");
        let arc: Arc<T> = obj.into();
        w.insert(arc.clone(), arc.clone());
        arc
    }
}
