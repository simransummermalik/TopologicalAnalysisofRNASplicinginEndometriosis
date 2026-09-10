use std::path::PathBuf;

use splice_girl::analyze_junctions;
use splice_girl::matrix::incidence_matrix;
use splice_girl::parser::parse_tsv;

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

#[test]
fn linear_graph_has_no_cycle_component() {
    let result = analyze_fixture("linear.tsv");
    assert!(result.decomposition.checks.passed);
    assert!(result.decomposition.cycle_fraction.abs() < 1e-12);
}

#[test]
fn equal_weight_directed_cycle_is_entirely_cyclic() {
    let result = analyze_fixture("cycle.tsv");
    assert!(result.decomposition.checks.passed);
    assert!((result.decomposition.cycle_fraction - 1.0).abs() < 1e-12);
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
