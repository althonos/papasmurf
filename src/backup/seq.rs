use std::str::FromStr;

use crate::alphabet::Dna;


/// An owned DNA sequence.
#[derive(Debug, Default, Clone, PartialEq, Hash, Eq)]
pub struct Sequence {
    data: Vec<Dna>,
    ambiguous: bool,
}

impl Sequence {
    /// Create a new sequence.
    pub fn new(data: Vec<Dna>) -> Self {
        let ambiguous = data.iter().any(|s| s.is_ambiguous());
        Self { data, ambiguous }
    }

    /// Get the reverse-complement of this sequence.
    pub fn to_reverse_complement(&self) -> Sequence {
        let rev = self.data.iter().rev().map(|s| s.complement()).collect();
        Sequence {
            ambiguous: self.ambiguous,
            data: rev,
        }
    }

    /// The length of the sequence
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check whether the sequence is ambiguous.
    pub fn is_ambiguous(&self) -> bool {
        self.ambiguous
    }

    pub fn as_slice(&self) -> &[Dna] {
        self.data.as_slice()
    }

    pub fn as_seq(&self) -> Seq<'_> {
        Seq {
            ambiguous: self.ambiguous,
            data: self.data.as_slice()
        }
    }
}

impl From<&[Dna]> for Sequence {
    fn from(data: &[Dna]) -> Self {
        let data = unsafe {
            std::slice::from_raw_parts(data.as_ptr() as *const Dna, data.len())
        };
        Self {
            data: data.to_owned(),
            ambiguous: false,
        }
    }
}

impl From<Vec<Dna>> for Sequence {
    fn from(data: Vec<Dna>) -> Self {
        Self::new(data)
    }
}

impl FromStr for Sequence {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s
            .chars()
            .map(|c| Dna::from_char(c))
            .collect::<Option<Vec<Dna>>>()
        {
            Some(data) => Ok(Sequence::new(data)),
            None => Err(()),
        }
    }
}


#[derive(Debug, Default, Clone, PartialEq, Hash, Eq)]
pub struct Seq<'a> {
    data: &'a [Dna],
    ambiguous: bool
}

impl<'a> Seq<'a> {
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn subseq(&self, start: usize, stop: usize) -> Seq<'a> {
        Seq {
            data: &self.data[start..stop],
            ambiguous: self.ambiguous,
        }
    }
}

impl<'a> AsRef<[Dna]> for Seq<'a> {
    fn as_ref(&self) -> &[Dna] {
        self.data
    }
}

// // --- 

pub trait Mismatches<T: ?Sized> {
    /// Count the number of mismatches between two sequences.
    fn mismatches(&self, other: &T) -> usize;
}

// impl Mismatches<[Dna]> for [Dna] {
//     fn mismatches(&self, other: &[Dna]) -> usize {
//         assert_eq!(self.len(), other.len());
//         self.iter()
//             .zip(other.iter())
//             .filter(|(x, y)| x != y)
//             .count()
//     }
// }

// impl Mismatches<[Dna]> for [Dna] {
//     fn mismatches(&self, other: &[Dna]) -> usize {
//         assert_eq!(self.len(), other.len());
//         self.iter()
//             .zip(other.iter())
//             .filter(|(x, y)| !y.matches(x))
//             .count()
//     }
// }

// impl Mismatches<[Dna]> for [Dna] {
//     fn mismatches(&self, other: &[Dna]) -> usize {
//         other.mismatches(self)
//     }
// }

impl Mismatches<Sequence> for Sequence {
    fn mismatches(&self, other: &Sequence) -> usize {
        self.as_seq().mismatches(&other.as_seq())
    }
}

impl<'a> Mismatches<Seq<'a>> for Sequence {
    fn mismatches(&self, other: &Seq<'a>) -> usize {
        self.as_seq().mismatches(other)
    }
}

impl<'a> Mismatches<Sequence> for Seq<'a> {
    fn mismatches(&self, other: &Sequence) -> usize {
        self.mismatches(&other.as_seq())
    }
}

impl<'a> Mismatches<Seq<'a>> for Seq<'a> {
    fn mismatches(&self, other: &Seq<'a>) -> usize {
        assert_eq!(self.len(), other.len());

        let mut mm = 0;
        if !self.ambiguous && !other.ambiguous {
            for i in 0..self.len() {
                if self.data[i] != other.data[i] {
                    mm += 1;
                }
            }
        } else {
            println!("{:?} {:?}", self.data, other.data);
            for i in 0..self.len() {
                if !self.data[i].matches(&other.data[i]) {
                    mm += 1;
                }
            }
        }
        mm
    }
}