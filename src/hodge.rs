use crate::graph::SpliceGraph;
use crate::matrix::{Matrix, dot, incidence_matrix, matrix_vector, norm, solve, transpose_vector};

pub const TOLERANCE: f64 = 1e-9;

#[derive(Debug)]
pub struct Checks {
    pub reconstruction_error: f64,
    pub kernel_error: f64,
    pub orthogonality_error: f64,
    pub passed: bool,
}

#[derive(Debug)]
pub struct Decomposition {
    pub gradient: Vec<f64>,
    pub cycle: Vec<f64>,
    pub gradient_energy: f64,
    pub cycle_energy: f64,
    pub signal_energy: f64,
    pub cycle_fraction: f64,
    pub checks: Checks,
}

pub fn decompose(graph: &SpliceGraph) -> Result<Decomposition, String> {
    let incidence = incidence_matrix(graph);
    let potential = fit_potential(graph, &incidence)?;
    let gradient = transpose_vector(&incidence, &potential);
    let cycle: Vec<f64> = graph
        .signal
        .iter()
        .zip(&gradient)
        .map(|(signal, fitted)| signal - fitted)
        .collect();

    let signal_energy = dot(&graph.signal, &graph.signal);
    let gradient_energy = dot(&gradient, &gradient);
    let cycle_energy = dot(&cycle, &cycle);
    if signal_energy <= 0.0 {
        return Err("normalized signal has zero energy".to_owned());
    }

    let reconstructed: Vec<f64> = gradient.iter().zip(&cycle).map(|(a, b)| a + b).collect();
    let reconstruction_delta: Vec<f64> = graph
        .signal
        .iter()
        .zip(reconstructed)
        .map(|(expected, actual)| expected - actual)
        .collect();
    let reconstruction_error = norm(&reconstruction_delta);
    let kernel_error = norm(&matrix_vector(&incidence, &cycle));
    let orthogonality_error = dot(&gradient, &cycle).abs();
    let scale = norm(&graph.signal).max(1.0);
    let passed = reconstruction_error <= TOLERANCE * scale
        && kernel_error <= TOLERANCE * scale
        && orthogonality_error <= TOLERANCE * scale * scale;

    Ok(Decomposition {
        gradient,
        cycle,
        gradient_energy,
        cycle_energy,
        signal_energy,
        cycle_fraction: cycle_energy / signal_energy,
        checks: Checks {
            reconstruction_error,
            kernel_error,
            orthogonality_error,
            passed,
        },
    })
}

fn fit_potential(graph: &SpliceGraph, incidence: &Matrix) -> Result<Vec<f64>, String> {
    let components = connected_components(graph);
    let mut is_root = vec![false; graph.nodes.len()];
    for component in components {
        is_root[component[0]] = true;
    }

    let free_nodes: Vec<usize> = (0..graph.nodes.len())
        .filter(|&node| !is_root[node])
        .collect();
    if free_nodes.is_empty() {
        return Ok(vec![0.0; graph.nodes.len()]);
    }

    let right_full = matrix_vector(incidence, &graph.signal);
    let mut laplacian = vec![vec![0.0; free_nodes.len()]; free_nodes.len()];
    let mut right = vec![0.0; free_nodes.len()];

    for (i, &row_node) in free_nodes.iter().enumerate() {
        right[i] = right_full[row_node];
        for (j, &column_node) in free_nodes.iter().enumerate() {
            laplacian[i][j] = incidence[row_node]
                .iter()
                .zip(&incidence[column_node])
                .map(|(a, b)| a * b)
                .sum();
        }
    }

    let free_potential = solve(laplacian, right)?;
    let mut potential = vec![0.0; graph.nodes.len()];
    for (&node, value) in free_nodes.iter().zip(free_potential) {
        potential[node] = value;
    }
    Ok(potential)
}

fn connected_components(graph: &SpliceGraph) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); graph.nodes.len()];
    for edge in &graph.edges {
        adjacency[edge.from].push(edge.to);
        adjacency[edge.to].push(edge.from);
    }

    let mut visited = vec![false; graph.nodes.len()];
    let mut components = Vec::new();
    for start in 0..graph.nodes.len() {
        if visited[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        visited[start] = true;
        while let Some(node) = stack.pop() {
            component.push(node);
            for &neighbor in &adjacency[node] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}
