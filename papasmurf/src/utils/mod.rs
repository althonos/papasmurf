mod interner;
mod ordered_set;
mod paired;

pub use self::interner::Interner;
pub use self::ordered_set::OrderedSet;
pub use self::paired::Paired;

pub type Rc<T> = std::sync::Arc<T>;
