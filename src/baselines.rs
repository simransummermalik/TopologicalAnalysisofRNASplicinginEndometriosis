//! Descriptive baselines on exactly the same eligible biological samples.
use crate::{
    GroupResult, hodge::connected_components, manifest::Condition, statistics::StatisticsRun,
};
use std::fmt::Write as _;

pub const METRICS: [&str; 4] = [
    "total_junction_support",
    "node_count",
    "junction_count",
    "cycle_rank",
];

pub fn measurements(result: &GroupResult) -> [f64; 4] {
    let graph = &result.graph;
    [
        graph.edges.iter().map(|edge| edge.count).sum(),
        graph.nodes.len() as f64,
        graph.edges.len() as f64,
        (graph.edges.len() + connected_components(graph).len() - graph.nodes.len()) as f64,
    ]
}

fn same_topology(left: &GroupResult, right: &GroupResult) -> bool {
    left.graph.nodes == right.graph.nodes
        && left.graph.edges.iter().map(|e| (e.from, e.to)).eq(right
            .graph
            .edges
            .iter()
            .map(|e| (e.from, e.to)))
}

fn field(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}
fn number(value: Option<f64>) -> String {
    value.map(|v| format!("{v:.10}")).unwrap_or_default()
}

/// Missing groups remain empty, never zero. Different observed topologies are
/// flagged, not silently matched, imputed, or declared biologically comparable.
pub fn tables(run_id: &str, results: &[GroupResult], run: &StatisticsRun) -> (String, String) {
    let mut samples = String::from(
        "run_id\tgene_id\tsample_id\tcondition\ttotal_junction_support\tnode_count\tjunction_count\tcycle_rank\n",
    );
    let mut comparisons = String::from(
        "run_id\tgene_id\tmetric\tcase_n\tcontrol_n\tcase_mean\tcontrol_mean\tmean_difference_case_minus_control\tobserved_topology\n",
    );
    for summary in &run.comparisons {
        let gene = &summary.comparison.gene_id;
        let mut groups: [Vec<[f64; 4]>; 2] = [Vec::new(), Vec::new()];
        let mut first = None;
        let mut identical = true;
        for sample in &summary.comparison.included {
            let result = results
                .iter()
                .find(|r| r.gene_id == *gene && r.sample_id == sample.sample_id)
                .expect("comparison includes only analyzed results");
            if let Some(previous) = first {
                identical &= same_topology(previous, result);
            } else {
                first = Some(result);
            }
            let values = measurements(result);
            let group = if sample.condition == Condition::Case {
                0
            } else {
                1
            };
            groups[group].push(values);
            writeln!(
                samples,
                "{}\t{}\t{}\t{}\t{:.10}\t{}\t{}\t{}",
                field(run_id),
                field(gene),
                field(&sample.sample_id),
                sample.condition,
                values[0],
                values[1],
                values[2],
                values[3]
            )
            .unwrap();
        }
        let topology = if first.is_none() {
            "no_eligible_samples"
        } else if groups.iter().map(Vec::len).sum::<usize>() < 2 {
            "single_sample"
        } else if identical {
            "same"
        } else {
            "different_review_required"
        };
        for (index, name) in METRICS.iter().enumerate() {
            let mean = |group: usize| {
                (!groups[group].is_empty()).then(|| {
                    groups[group].iter().map(|v| v[index]).sum::<f64>() / groups[group].len() as f64
                })
            };
            let case = mean(0);
            let control = mean(1);
            writeln!(
                comparisons,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                field(run_id),
                field(gene),
                name,
                groups[0].len(),
                groups[1].len(),
                number(case),
                number(control),
                number(case.zip(control).map(|(a, b)| a - b)),
                topology
            )
            .unwrap();
        }
    }
    (samples, comparisons)
}
