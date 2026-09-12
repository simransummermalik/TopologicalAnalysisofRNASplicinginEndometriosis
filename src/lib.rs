//! Library entry point tying the pipeline together: parsed junctions ->
//! per-(sample, gene) graphs -> Hodge decomposition. `main.rs` is a thin
//! CLI wrapper around `analyze_junctions`.
pub mod graph;
pub mod hodge;
pub mod matrix;
pub mod parser;

use std::collections::BTreeMap;

use graph::SpliceGraph;
use hodge::{Decomposition, decompose};
use parser::Junction;

/// The graph and decomposition computed for one sample's one gene.
#[derive(Debug)]
pub struct GroupResult {
    pub sample_id: String,
    pub gene_id: String,
    pub graph: SpliceGraph,
    pub decomposition: Decomposition,
}

/// Groups junction rows by (sample_id, gene_id) — the decomposition is
/// always computed within a single sample/gene, never pooled across
/// samples — then builds a graph and runs the decomposition for each group.
pub fn analyze_junctions(junctions: Vec<Junction>) -> Result<Vec<GroupResult>, String> {
    let mut groups: BTreeMap<(String, String), Vec<Junction>> = BTreeMap::new();

    for junction in junctions {
        groups
            .entry((junction.sample_id.clone(), junction.gene_id.clone()))
            .or_default()
            .push(junction);
    }

    let mut results = Vec::with_capacity(groups.len());
    for ((sample_id, gene_id), rows) in groups {
        let graph = SpliceGraph::from_junctions(&rows)?;
        let decomposition = decompose(&graph)?;
        results.push(GroupResult {
            sample_id,
            gene_id,
            graph,
            decomposition,
        });
    }

    Ok(results)
}
