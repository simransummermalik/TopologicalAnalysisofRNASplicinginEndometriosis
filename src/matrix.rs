use crate::graph::SpliceGraph;

pub type Matrix = Vec<Vec<f64>>;

pub fn incidence_matrix(graph: &SpliceGraph) -> Matrix {
    let mut matrix = vec![vec![0.0; graph.edges.len()]; graph.nodes.len()];
    for (edge_index, edge) in graph.edges.iter().enumerate() {
        matrix[edge.from][edge_index] = -1.0;
        matrix[edge.to][edge_index] = 1.0;
    }
    matrix
}

pub fn matrix_vector(matrix: &Matrix, vector: &[f64]) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
        .collect()
}

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

        let pivot_tail = matrix[pivot][pivot..].to_vec();
        for row in (pivot + 1)..n {
            let factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (entry, pivot_entry) in matrix[row][pivot..].iter_mut().zip(&pivot_tail) {
                *entry -= factor * pivot_entry;
            }
            values[row] -= factor * values[pivot];
        }
    }

    let mut solution = vec![0.0; n];
    for row in (0..n).rev() {
        let known: f64 = ((row + 1)..n)
            .map(|column| matrix[row][column] * solution[column])
            .sum();
        solution[row] = (values[row] - known) / matrix[row][row];
    }
    Ok(solution)
}
