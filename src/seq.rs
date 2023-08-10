use crate::utils::Rc;

pub fn is_ambiguous(c: char) -> bool {
    match c {
        'A' | 'T' | 'C' | 'G' | 'U' => false,
        'R' | 'Y' | 'S' | 'W' | 'M' | 'K' | 'B' | 'D' | 'H' | 'V' | 'N' => true,
        _ => unreachable!(),
    }
}

pub fn count_ambiguous(s: &str) -> usize {
    s.chars().filter(|&c| is_ambiguous(c)).count()
}

pub fn reverse_complement(s: &str) -> String {
    let mut rev = String::with_capacity(s.len());
    for x in s.chars().rev() {
        rev.push(match x {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',

            'Y' => 'R',
            'R' => 'Y',

            'S' => 'S',
            'W' => 'W',

            'M' => 'K',
            'K' => 'M',

            'B' => 'V',
            'D' => 'H',
            'H' => 'D',
            'V' => 'B',

            'N' => 'N',

            _ => panic!(),
        })
    }
    rev
}

pub fn dna_match(c1: char, c2: char) -> bool {
    match c1 {
        'A' => c2 == 'A',
        'T' => c2 == 'T',
        'C' => c2 == 'C',
        'G' => c2 == 'G',

        'R' => c2 == 'A' || c2 == 'G',
        'Y' => c2 == 'C' || c2 == 'T',

        'S' => c2 == 'G' || c2 == 'G',
        'W' => c2 == 'A' || c2 == 'T',

        'K' => c2 == 'G' || c2 == 'T',
        'M' => c2 == 'A' || c2 == 'C',

        'B' => c2 == 'C' || c2 == 'G' || c2 == 'T',
        'D' => c2 == 'A' || c2 == 'G' || c2 == 'T',
        'H' => c2 == 'A' || c2 == 'C' || c2 == 'T',
        'V' => c2 == 'A' || c2 == 'C' || c2 == 'G',

        'N' => c2 == 'A' || c2 == 'G' || c2 == 'T' || c2 == 'C',
        _ => unreachable!(),
    }
}

pub fn mismatches(s1: &str, s2: &str) -> usize {
    s1.chars()
        .zip(s2.chars())
        .filter(|(x, y)| !dna_match(*x, *y))
        .count()
}

pub struct DesambiguationIterator<'a> {
    sequence: &'a str,
    buffer: Rc<String>,
    pos: Vec<usize>,
    state: Vec<usize>,
    variants: Vec<&'static str>,
    remaining: usize,
}

impl<'a> DesambiguationIterator<'a> {
    pub fn new(sequence: &'a str) -> Self {
        let mut buffer = Rc::new(String::new());
        let mut pos = Vec::new();
        let mut state = Vec::new();
        let mut variants = Vec::new();
        let mut remaining = 1;

        for (i, c) in sequence
            .chars()
            .enumerate()
            .filter(|(_, c)| is_ambiguous(*c))
        {
            let var = match c {
                'R' => "AG",
                'Y' => "CT",
                'S' => "GC",
                'W' => "AT",
                'K' => "GT",
                'M' => "AC",
                'B' => "CGT",
                'D' => "AGT",
                'H' => "ACT",
                'V' => "ACG",
                'N' => "ACGT",
                _ => unreachable!(),
            };

            remaining *= var.len();
            variants.push(var);
            pos.push(i);
            state.push(var.len() - 1);
        }

        DesambiguationIterator {
            buffer,
            pos,
            state,
            variants,
            remaining,
            sequence,
        }
    }
}

impl<'a> Iterator for DesambiguationIterator<'a> {
    type Item = Rc<String>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        let buffer = match Rc::get_mut(&mut self.buffer) {
            Some(b) => b,
            None => {
                self.buffer = Rc::new(self.sequence.to_string());
                Rc::get_mut(&mut self.buffer).unwrap()
            }
        };

        if buffer.is_empty() {
            buffer.push_str(self.sequence);
        }

        for i in 0..self.pos.len() {
            self.state[i] += 1;
            if self.state[i] == self.variants[i].len() {
                self.state[i] = 0;
            } else {
                break;
            }
        }

        for i in 0..self.pos.len() {
            buffer.replace_range(
                self.pos[i]..self.pos[i] + 1,
                &self.variants[i][self.state[i]..self.state[i] + 1],
            );
        }

        self.remaining -= 1;
        Some(self.buffer.clone())
    }
}

pub fn disambiguate<'a>(s: &'a str) -> Vec<String> {
    DesambiguationIterator::new(s)
        .map(|s| s.clone().as_ref().to_owned())
        .collect()
}
