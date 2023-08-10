use std::io::BufRead;
use std::io::Read;
use std::io::BufReader;
use std::io::Write;

// --- FASTQ -------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct FastqRecord {
    pub id: String,
    pub sequence: String,
    pub strand: String,
    pub quality: String,
}

#[derive(Debug, Clone)]
pub struct FastqReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }
}

impl<R: Read> From<R> for FastqReader<BufReader<R>> {
    fn from(reader: R) -> Self {
        Self::new(BufReader::new(reader))
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<FastqRecord, std::io::Error>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut record = FastqRecord::default();
        match self.reader.read_line(&mut record.id) {
            Ok(0) => return None,
            Err(e) => return Some(Err(e)),
            Ok(_) => (),
        }
        record.id.pop();

        match self.reader.read_line(&mut record.sequence) {
            Ok(0) => return None,
            Err(e) => return Some(Err(e)),
            Ok(_) => (),
        }
        record.sequence.pop();

        match self.reader.read_line(&mut record.strand) {
            Ok(0) => return None,
            Err(e) => return Some(Err(e)),
            Ok(_) => (),
        }
        record.strand.pop();

        match self.reader.read_line(&mut record.quality) {
            Ok(0) => return None,
            Err(e) => return Some(Err(e)),
            Ok(_) => (),
        }
        record.quality.pop();

        Some(Ok(record))
    }
}

// --- FASTA -------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

#[derive(Debug)]
pub struct FastaReader<R: BufRead> {
    reader: R,
    cache: Option<Result<String, std::io::Error>>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(mut reader: R) -> Self {
        let mut buffer = String::new();
        if let Err(e) = reader.read_line(&mut buffer) {
            Self {
                reader,
                cache: Some(Err(e)),
            }
        } else {
            Self {
                reader,
                cache: Some(Ok(buffer.trim_start_matches('>').trim_end().to_string())),
            }
        }
    }
}

impl<R: Read> From<R> for FastaReader<BufReader<R>> {
    fn from(reader: R) -> Self {
        Self::new(BufReader::new(reader))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FastaRecord, std::io::Error>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut id = match self.cache.take()? {
            Ok(id) => id,
            Err(e) => return Some(Err(e)),
        };

        let mut sequence = String::new();
        let mut end = 0;
        let mut tmp = String::new();

        loop {
            tmp.clear();
            match self.reader.read_line(&mut sequence) {
                Err(e) => return Some(Err(e)),
                Ok(0) => break,
                Ok(n) => {
                    if sequence[end..].starts_with('>') {
                        self.cache = Some(Ok(sequence[end+1..].trim_end().to_string()));
                        sequence.truncate(end);
                        break;
                    } else {
                        if sequence.ends_with('\n') {
                            sequence.pop();
                        }
                        end = sequence.len();
                    }

                }
            }
        }

        Some(Ok(FastaRecord { id, sequence }))
    }
}