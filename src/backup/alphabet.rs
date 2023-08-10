// --- Dna ------------------------------------------------------------

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Dna {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    #[default]
    N,
}

impl Dna {
    pub const fn from_ascii(c: u8) -> Option<Self> {
        match c.to_ascii_uppercase() {
            b'A' => Some(Dna::A),
            b'C' => Some(Dna::C),
            b'G' => Some(Dna::G),
            b'T' | b'U' => Some(Dna::T),
            b'R' => Some(Dna::R),
            b'Y' => Some(Dna::Y),
            b'S' => Some(Dna::S),
            b'W' => Some(Dna::W),
            b'K' => Some(Dna::K),
            b'M' => Some(Dna::M),
            b'B' => Some(Dna::B),
            b'D' => Some(Dna::D),
            b'H' => Some(Dna::H),
            b'V' => Some(Dna::V),
            b'N' => Some(Dna::N),
            _ => None,
        }
    }

    pub const fn from_char(c: char) -> Option<Self> {
        if c.is_ascii() {
            Self::from_ascii(c as u8)
        } else {
            None
        }
    }

    pub const fn to_dna(&self) -> Option<Dna> {
        match self {
            Dna::A => Some(Dna::A),
            Dna::T => Some(Dna::T),
            Dna::G => Some(Dna::G),
            Dna::C => Some(Dna::C),
            _ => None,
        }
    }

    pub const fn is_ambiguous(&self) -> bool {
        match self {
            Dna::A | Dna::C | Dna::G | Dna::T => false,
            _ => true,
        }
    }

    #[inline]
    pub const fn matches(&self, other: &Dna) -> bool {
        match self {
            Dna::A => *other as u8 == Dna::A as u8,
            Dna::T => *other as u8 == Dna::T as u8,
            Dna::G => *other as u8 == Dna::G as u8,
            Dna::C => *other as u8 == Dna::C as u8,

            Dna::R => *other as u8 == Dna::A as u8 || *other as u8 == Dna::G as u8,
            Dna::Y => *other as u8 == Dna::C as u8 || *other as u8 == Dna::T as u8,

            Dna::S => *other as u8 == Dna::G as u8 || *other as u8 == Dna::C as u8,
            Dna::W => *other as u8 == Dna::A as u8 || *other as u8 == Dna::T as u8,

            Dna::K => *other as u8 == Dna::G as u8 || *other as u8 == Dna::T as u8,
            Dna::M => *other as u8 == Dna::A as u8 || *other as u8 == Dna::C as u8,

            Dna::B => {
                *other as u8 == Dna::C as u8
                    || *other as u8 == Dna::G as u8
                    || *other as u8 == Dna::T as u8
            }
            Dna::D => {
                *other as u8 == Dna::A as u8
                    || *other as u8 == Dna::G as u8
                    || *other as u8 == Dna::T as u8
            }
            Dna::H => {
                *other as u8 == Dna::A as u8
                    || *other as u8 == Dna::C as u8
                    || *other as u8 == Dna::T as u8
            }
            Dna::V => {
                *other as u8 == Dna::A as u8
                    || *other as u8 == Dna::C as u8
                    || *other as u8 == Dna::G as u8
            }

            Dna::N => true,
        }
    }

    pub const fn complement(&self) -> Self {
        match self {
            Dna::A => Dna::T,
            Dna::C => Dna::G,
            Dna::G => Dna::C,
            Dna::T => Dna::A,

            Dna::Y => Dna::R,
            Dna::R => Dna::Y,

            Dna::S => Dna::S,
            Dna::W => Dna::W,

            Dna::M => Dna::K,
            Dna::K => Dna::M,

            Dna::B => Dna::V,
            Dna::D => Dna::H,
            Dna::H => Dna::D,
            Dna::V => Dna::B,

            Dna::N => Dna::N,
        }
    }
}

// // --- Dna ---------------------------------------------------------------------

// #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
// #[repr(u8)]
// pub enum Dna {
//     A = 0,
//     C = 1,
//     G = 2,
//     T = 3,
// }

// impl Dna {
//     pub const fn from_ascii(c: u8) -> Option<Self> {
//         match c.to_ascii_uppercase() {
//             b'A' => Some(Dna::A),
//             b'C' => Some(Dna::C),
//             b'G' => Some(Dna::G),
//             b'T' | b'U' => Some(Dna::T),
//             _ => None,
//         }
//     }

//     pub const fn from_char(c: char) -> Option<Self> {
//         if c.is_ascii() {
//             Self::from_ascii(c as u8)
//         } else {
//             None
//         }
//     }

//     pub const fn to_ambiguous(&self) -> Dna {
//         match self {
//             Dna::A => Dna::A,
//             Dna::T => Dna::T,
//             Dna::G => Dna::G,
//             Dna::C => Dna::C,
//         }
//     }

//     pub const fn complement(&self) -> Self {
//         match self {
//             Dna::A => Dna::T,
//             Dna::C => Dna::G,
//             Dna::G => Dna::C,
//             Dna::T => Dna::A,
//         }
//     }
// }
