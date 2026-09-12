//! The core Hodge decomposition, matching `docs/math_spec.md`. Splits a
//! graph's normalized junction-usage signal `F` into:
//!   - a gradient part `F_grad = B1^T phi`, explained by a single value
//!     ("potential") assigned to each node, and
//!   - a cycle-space part `F_cycle = F - F_grad`, the leftover signal that
//!     can't be explained by any per-node potential (e.g. exon-skipping
//!     loops, back-splicing).
//! `cycle_fraction = ||F_cycle||^2 / ||F||^2` is the number this project
//! ultimately compares between endometriosis and control samples.
use crate::graph::SpliceGraph;
use crate::matrix::{Matrix, dot, incidence_matrix, matrix_vector, norm, solve, transpose_vector};

pub const TOLERANCE: f64 = 1e-9;

/// Results of the three numerical sanity checks every decomposition must
/// pass before its output is trusted (see docs/math_spec.md "Acceptance
/// checks").
#[derive(Debug)]
pub struct Checks {
    /// How far `gradient + cycle` is from the original signal.
    pub reconstruction_error: f64,
    /// How far `B1 * cycle` is from zero (cycle must lie in ker(B1)).
    pub kernel_error: f64,
    /// How far `gradient . cycle` is from zero (the two parts must be
    /// orthogonal, i.e. genuinely separate signal).
    pub orthogonality_error: f64,
    pub passed: bool,
}

/// The full per-sample, per-gene decomposition result.
#[derive(Debug)]
pub struct Decomposition {
    pub gradient: Vec<f64>,
    pub cycle: Vec<f64>,
    pub gradient_energy: f64,
    pub cycle_energy: f64,
    pub signal_energy: f64,
    /// Share of signal energy that lies in the cycle space; the project's
    /// headline per-sample measurement.
    pub cycle_fraction: f64,
    pub checks: Checks,
}

/// Runs the decomposition for one graph's normalized signal and verifies
/// the result numerically before returning it.
pub fn decompose(graph: &SpliceGraph) -> Result<Decomposition, String> {
    let incidence = incidence_matrix(graph);
    let potential = fit_potential(graph, &incidence)?;
    // F_grad = B1^T * phi: the edge-wise difference implied by the fitted
    // per-node potential.
    let gradient = transpose_vector(&incidence, &potential);
    // F_cycle = F - F_grad: whatever the potential couldn't explain.
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

    // Verify the three identities from docs/math_spec.md rather than
    // trusting the solver blindly; failing loudly here is safer than
    // reporting a cycle_fraction computed from a broken projection.
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

/// Finds the per-node potential `phi` minimizing `||F - B1^T phi||^2` by
/// solving its normal equations `(B1 B1^T) phi = B1 F`. `B1 B1^T` is the
/// graph Laplacian, which is singular (rank-deficient by one per connected
/// component: adding a constant to every potential in a component doesn't
/// change `B1^T phi`). To get a solvable system we pin one node per
/// component ("root") to potential zero and solve only for the rest
/// ("free" nodes) — the resulting F_grad is the same regardless of which
/// root is chosen, since it only depends on differences between potentials.
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
    // Restrict the Laplacian (B1 B1^T) and right-hand side (B1 F) to the
    // free nodes' rows/columns; the root nodes' potentials are fixed at 0.
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

/// Groups nodes into connected components using the graph's *undirected*
/// adjacency (edge direction doesn't matter for reachability here) via
/// iterative depth-first search. Each component needs its own root node in
/// `fit_potential` because potentials from different components can't be
/// compared or linked by any edge.
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
