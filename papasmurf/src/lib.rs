#![doc = include_str!("../README.md")]

mod db;
mod errors;
mod primer;
mod seq;
mod utils;

pub mod io;
pub mod matrix;

pub use self::db::Builder;
pub use self::db::Database;
pub use self::db::Mapper;
pub use self::db::MapperResult;
pub use self::db::Region;
pub use self::errors::Error;
pub use self::primer::Primer;
pub use self::utils::Paired;
