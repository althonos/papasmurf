use std::fmt::Display;
use std::fmt::Error as FmtError;
use std::fmt::Formatter;

#[derive(Debug)]
pub enum Error {
    InvalidDna,
    InvalidDimensions,
    InvalidKmerLength,
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), FmtError> {
        match *self {
            Error::InvalidDna => f.write_str("invalid DNA"),
            Error::InvalidDimensions => f.write_str("invalid dimensions"),
            Error::InvalidKmerLength => f.write_str("invalid k-mer length"),
        }
    }
}

impl std::error::Error for Error {}
