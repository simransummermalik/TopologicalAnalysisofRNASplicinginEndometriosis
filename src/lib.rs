pub mod graph;
pub mod hodge;
pub mod matrix;
pub mod parser;

use std::collections::BTreeMap;

use graph::SpliceGraph;
use hodge::{Decomposition, decompose};
use parser::Junction;

#[derive(Debug)]
pub struct GroupResult {
    pub sample_id: String,
    pub gene_id: String,
    pub graph: SpliceGraph,
    pub decomposition: Decomposition,
}

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
