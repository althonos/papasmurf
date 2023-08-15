use serde::Deserialize;
use serde::Serialize;

use super::MatrixDimensions;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CscMatrix<T> {
    pub(super) rows: usize,
    pub(super) data: Vec<T>,
    pub(super) col_index: Vec<usize>,
    pub(super) row_index: Vec<usize>,
}

impl<T> CscMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        CscMatrix {
            rows,
            data: Vec::new(),
            row_index: Vec::new(),
            col_index: vec![0; cols + 1],
        }
    }

    pub fn nnz(&self) -> usize {
        self.data.len()
    }
}

impl<T> Default for CscMatrix<T> {
    fn default() -> Self {
        Self {
            rows: 0,
            data: Vec::new(),
            col_index: Vec::new(),
            row_index: Vec::new(),
        }
    }
}

impl<T> MatrixDimensions for CscMatrix<T> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn columns(&self) -> usize {
        self.col_index.len() - 1
    }
}
