//! A small dense linear-algebra toolkit used only by `hodge.rs`. Kept
//! dependency-free (no external numerical crate) per the roadmap's MVP
//! priority of correctness over performance; revisit if graphs get large.
use crate::graph::SpliceGraph;

/// Dense matrix stored as rows of `f64`.
pub type Matrix = Vec<Vec<f64>>;

/// Builds the signed vertex-edge incidence matrix B1 (nodes x edges): for
/// each edge `from -> to`, column `edge_index` has -1 at row `from` and +1
/// at row `to`. This encodes the graph's orientation so that `B1^T phi`
/// (see `transpose_vector`) gives the edge-wise difference `phi[to] -
/// phi[from]` for any per-node potential `phi`.
pub fn incidence_matrix(graph: &SpliceGraph) -> Matrix {
    let mut matrix = vec![vec![0.0; graph.edges.len()]; graph.nodes.len()];
    for (edge_index, edge) in graph.edges.iter().enumerate() {
        matrix[edge.from][edge_index] = -1.0;
        matrix[edge.to][edge_index] = 1.0;
    }
    matrix
}

/// Computes `matrix * vector` (one output per row).
pub fn matrix_vector(matrix: &Matrix, vector: &[f64]) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
        .collect()
}

/// Computes `matrix^T * vector` (one output per column) without
/// materializing the transpose.
pub fn transpose_vector(matrix: &Matrix, vector: &[f64]) -> Vec<f64> {
    let columns = matrix.first().map_or(0, Vec::len);
    (0..columns)
        .map(|column| {
            matrix
                .iter()
                .zip(vector)
                .map(|(row, value)| row[column] * value)
                .sum()
        })
        .collect()
}

pub fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

pub fn norm(vector: &[f64]) -> f64 {
    dot(vector, vector).sqrt()
}

/// Solves the square linear system `matrix * x = values` for `x` using
/// Gaussian elimination with partial pivoting (choosing the largest
/// available pivot in each column improves numerical stability). Used by
/// `hodge::fit_potential` to solve the least-squares normal equations for
/// the per-node potential. Returns an error if the system is singular
/// (e.g. the graph has no edges reaching some free node).
pub fn solve(mut matrix: Matrix, mut values: Vec<f64>) -> Result<Vec<f64>, String> {
    let n = values.len();
    if matrix.len() != n || matrix.iter().any(|row| row.len() != n) {
        return Err("linear system has inconsistent dimensions".to_owned());
    }

    for pivot in 0..n {
        let best = (pivot..n)
            .max_by(|&left, &right| {
                matrix[left][pivot]
                    .abs()
                    .total_cmp(&matrix[right][pivot].abs())
            })
            .expect("pivot range is nonempty");
        if matrix[best][pivot].abs() < 1e-12 {
            return Err("projection system is singular at the configured tolerance".to_owned());
        }
        matrix.swap(pivot, best);
        values.swap(pivot, best);

        // Eliminate this column from all rows below the pivot.
        let pivot_tail = matrix[pivot][pivot..].to_vec();
        for row in (pivot + 1)..n {
            let factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (entry, pivot_entry) in matrix[row][pivot..].iter_mut().zip(&pivot_tail) {
                *entry -= factor * pivot_entry;
            }
            values[row] -= factor * values[pivot];
        }
    }

    // Matrix is now upper triangular; solve from the last row upward.
    let mut solution = vec![0.0; n];
    for row in (0..n).rev() {
        let known: f64 = ((row + 1)..n)
            .map(|column| matrix[row][column] * solution[column])
            .sum();
        solution[row] = (values[row] - known) / matrix[row][row];
    }
    Ok(solution)
}
