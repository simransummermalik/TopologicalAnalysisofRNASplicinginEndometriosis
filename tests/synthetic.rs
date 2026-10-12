use std::collections::BTreeMap;
use std::path::PathBuf;

use splice_girl::analyze_junctions;
use splice_girl::matrix::incidence_matrix;
use splice_girl::parser::{parse_tsv, parse_tsv_text};

fn analyze_fixture(name: &str) -> splice_girl::GroupResult {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(name);
    let rows = parse_tsv(&path).expect("fixture should parse");
    let mut results = analyze_junctions(rows).expect("fixture should analyze");
    assert_eq!(results.len(), 1);
    results.remove(0)
}

fn connected_component_count(result: &splice_girl::GroupResult) -> usize {
    let mut adjacency = vec![Vec::new(); result.graph.nodes.len()];
    for edge in &result.graph.edges {
        adjacency[edge.from].push(edge.to);
        adjacency[edge.to].push(edge.from);
    }

    let mut visited = vec![false; adjacency.len()];
    let mut components = 0;
    for start in 0..adjacency.len() {
        if visited[start] {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(node) = stack.pop() {
            for &neighbor in &adjacency[node] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
    }
    components
}

fn numerical_rank(input: &[Vec<f64>]) -> usize {
    if input.is_empty() {
        return 0;
    }
    let mut matrix = input.to_vec();
    let rows = matrix.len();
    let columns = matrix[0].len();
    let mut rank = 0;

    for column in 0..columns {
        let Some(pivot) = (rank..rows).max_by(|&left, &right| {
            matrix[left][column]
                .abs()
                .total_cmp(&matrix[right][column].abs())
        }) else {
            break;
        };
        if matrix[pivot][column].abs() < 1e-12 {
            continue;
        }
        matrix.swap(rank, pivot);
        let pivot_value = matrix[rank][column];
        let pivot_row = matrix[rank].clone();
        for (row_index, row) in matrix.iter_mut().enumerate() {
            if row_index == rank {
                continue;
            }
            let factor = row[column] / pivot_value;
            for entry in column..columns {
                row[entry] -= factor * pivot_row[entry];
            }
        }
        rank += 1;
        if rank == rows {
            break;
        }
    }
    rank
}

#[test]
fn linear_graph_has_no_cycle_component() {
    let result = analyze_fixture("linear.tsv");
    assert!(result.decomposition.checks.passed);
    assert!(result.decomposition.cycle_fraction.abs() < 1e-12);
}

#[test]
fn shuffled_rows_preserve_graph_and_metrics() {
    let ordered = parse_tsv_text(
        "sample_id\tgene_id\tfrom\tto\tcount\n\
         sample_linear\tGENE_LINEAR\tE1\tE2\t100\n\
         sample_linear\tGENE_LINEAR\tE2\tE3\t80\n",
    )
    .expect("ordered input should parse");
    let shuffled = parse_tsv_text(
        "sample_id\tgene_id\tfrom\tto\tcount\n\
         sample_linear\tGENE_LINEAR\tE2\tE3\t80\n\
         sample_linear\tGENE_LINEAR\tE1\tE2\t100\n",
    )
    .expect("shuffled input should parse");

    let ordered = analyze_junctions(ordered)
        .expect("ordered input should analyze")
        .remove(0);
    let shuffled = analyze_junctions(shuffled)
        .expect("shuffled input should analyze")
        .remove(0);

    assert_eq!(ordered.graph.nodes, shuffled.graph.nodes);
    assert_eq!(ordered.graph.edges.len(), shuffled.graph.edges.len());
    for (left, right) in ordered.graph.edges.iter().zip(&shuffled.graph.edges) {
        assert_eq!((left.from, left.to), (right.from, right.to));
        assert!((left.count - right.count).abs() < 1e-12);
    }
    for (left, right) in ordered.graph.signal.iter().zip(&shuffled.graph.signal) {
        assert!((left - right).abs() < 1e-12);
    }
    assert!(
        (ordered.decomposition.cycle_fraction - shuffled.decomposition.cycle_fraction).abs()
            < 1e-12
    );
}

#[test]
fn incidence_dimensions_match_cycle_rank() {
    for fixture in [
        "linear.tsv",
        "cycle.tsv",
        "exon_skip.tsv",
        "disconnected.tsv",
    ] {
        let result = analyze_fixture(fixture);
        let incidence = incidence_matrix(&result.graph);
        let node_count = result.graph.nodes.len();
        let edge_count = result.graph.edges.len();
        assert_eq!(incidence.len(), node_count);
        assert!(incidence.iter().all(|row| row.len() == edge_count));

        let rank = numerical_rank(&incidence);
        let components = connected_component_count(&result);
        let expected_cycle_rank = edge_count + components - node_count;
        assert_eq!(edge_count - rank, expected_cycle_rank);
    }
}

#[test]
fn equal_weight_directed_cycle_is_entirely_cyclic() {
    let result = analyze_fixture("cycle.tsv");
    assert!(result.decomposition.checks.passed);
    assert!((result.decomposition.cycle_fraction - 1.0).abs() < 1e-12);
}

#[test]
fn branching_tree_has_no_cycle_component() {
    let result = analyze_fixture("branching_tree.tsv");
    assert!(result.decomposition.checks.passed);
    assert!(result.decomposition.cycle_fraction.abs() < 1e-12);
}

#[test]
fn exon_skipping_graph_has_nonzero_cycle_component_despite_no_directed_cycle() {
    // E1 -> E2 -> E3 plus the E1 -> E3 skip edge is acyclic as a directed
    // graph (no directed path leads back to E1), but as an undirected
    // triangle it has one independent cycle. With equal edge weights, hand
    // solving the normal equations (pinning E1's potential to zero) gives
    // cycle = (w/3, -w/3, w/3) for any w, so cycle_fraction is exactly 1/9
    // regardless of the total read support. M11 exists to confirm the
    // engine agrees with that independently derived answer.
    let result = analyze_fixture("exon_skip.tsv");
    assert!(result.decomposition.checks.passed);
    assert!((result.decomposition.cycle_fraction - 1.0 / 9.0).abs() < 1e-12);
}

#[test]
fn fixed_noise_fixture_matches_independent_expectation() {
    let result = analyze_fixture("noise.tsv");
    assert!(result.decomposition.checks.passed);
    let expected = 5618.0 / 45027.0;
    assert!((result.decomposition.cycle_fraction - expected).abs() < 1e-12);
}

#[test]
fn back_splice_weight_sweep_matches_reviewed_expectations() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("weight_sweep.tsv");
    let rows = parse_tsv(&path).expect("weight sweep should parse");
    let results = analyze_junctions(rows).expect("weight sweep should analyze");
    let by_sample: BTreeMap<_, _> = results
        .into_iter()
        .map(|result| (result.sample_id.clone(), result))
        .collect();
    let expected = [
        ("sample_weight_0_01", 4489.0 / 6667.0),
        ("sample_weight_0_05", 1681.0 / 2403.0),
        ("sample_weight_0_10", 49.0 / 67.0),
        ("sample_weight_0_20", 121.0 / 153.0),
        ("sample_weight_0_40", 8.0 / 9.0),
    ];
    for (sample, expected_fraction) in expected {
        let result = by_sample.get(sample).expect("expected sweep sample");
        assert!(result.decomposition.checks.passed);
        assert!((result.decomposition.cycle_fraction - expected_fraction).abs() < 1e-12);
    }
}

#[test]
fn incidence_columns_use_negative_source_and_positive_destination() {
    let result = analyze_fixture("linear.tsv");
    assert_eq!(result.graph.nodes, ["E1", "E2", "E3"]);
    assert_eq!(
        incidence_matrix(&result.graph),
        vec![vec![-1.0, 0.0], vec![1.0, -1.0], vec![0.0, 1.0]]
    );
}

#[test]
fn duplicate_directed_edges_are_rejected() {
    let input = "sample_id\tgene_id\tfrom\tto\tcount\n\
                 s1\tg1\tE1\tE2\t10\n\
                 s1\tg1\tE1\tE2\t5\n";
    let rows = parse_tsv_text(input).expect("duplicate rows have valid TSV syntax");
    let error = analyze_junctions(rows).expect_err("duplicate edges must be rejected");
    assert!(error.contains("duplicate edge E1 -> E2"));
}

#[test]
fn zero_total_support_is_rejected() {
    let input = "sample_id\tgene_id\tfrom\tto\tcount\n\
                 s1\tg1\tE1\tE2\t0\n";
    let rows = parse_tsv_text(input).expect("zero counts have valid TSV syntax");
    let error = analyze_junctions(rows).expect_err("zero-total groups must be rejected");
    assert!(error.contains("positive total support"));
}

#[test]
fn edge_case_fixtures_have_explicit_statuses() {
    let zero = parse_tsv(
        &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/zero_support.tsv"),
    )
    .expect("zero-support fixture should parse before analysis");
    assert!(analyze_junctions(zero).is_err());

    let duplicate = parse_tsv(
        &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/duplicate_edge.tsv"),
    )
    .expect("duplicate fixture should parse before graph validation");
    assert!(
        analyze_junctions(duplicate)
            .expect_err("duplicate edge should be rejected")
            .contains("duplicate edge")
    );

    let self_loop = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/self_loop.tsv");
    assert!(parse_tsv(&self_loop).is_err());

    let disconnected = analyze_fixture("disconnected.tsv");
    assert!(disconnected.decomposition.checks.passed);
    assert!(disconnected.decomposition.cycle_fraction.abs() < 1e-12);
}
