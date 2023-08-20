mod interner;
mod ordered_set;
mod paired;

pub use self::interner::Interner;
pub use self::ordered_set::OrderedSet;
pub use self::paired::Paired;

/// The reference counter type used in this crate.
pub type Rc<T> = std::sync::Arc<T>;
