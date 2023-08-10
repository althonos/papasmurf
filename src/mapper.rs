use super::db::Database;
use super::matrix::DokMatrix;

#[derive(Debug, Clone)]
pub struct Mapper<'db> {
    pub db: &'db Database,
    pub expected: Vec<DokMatrix<f32>>,
}

impl<'db> Mapper<'db> {
    pub fn new(db: &'db Database) -> Self {
        let expected = db
            .regions
            .iter()
            .map(|region| DokMatrix::new(0, region.unique_pairs.len()))
            .collect();
        Self { expected, db }
    }
}
