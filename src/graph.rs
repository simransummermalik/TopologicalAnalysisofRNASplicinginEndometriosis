//! Turns the junction rows for one (sample, gene) group into a directed
//! graph with a normalized edge signal, ready for the incidence-matrix and
//! Hodge-decomposition steps in `matrix.rs` / `hodge.rs`.
use std::collections::{BTreeMap, BTreeSet};

use crate::parser::Junction;

/// A directed edge between two node indices into `SpliceGraph::nodes`.
#[derive(Clone, Debug)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub count: f64,
}

/// A splice graph for one sample and gene: nodes are exons/splice sites,
/// edges are observed junctions, and `signal[i]` is the normalized usage
/// (`count / total`) of `edges[i]`.
#[derive(Debug)]
pub struct SpliceGraph {
    pub nodes: Vec<String>,
    pub edges: Vec<Edge>,
    pub signal: Vec<f64>,
}

impl SpliceGraph {
    /// Builds a graph from one sample/gene's junction rows. Node and edge
    /// order is derived from sorted identifiers (not input row order) so
    /// the resulting incidence matrix and results are deterministic
    /// regardless of how the TSV was written.
    pub fn from_junctions(rows: &[Junction]) -> Result<Self, String> {
        let total: f64 = rows.iter().map(|row| row.count).sum();
        if !total.is_finite() || total <= 0.0 {
            return Err("each sample/gene group must have positive total support".to_owned());
        }

        let mut node_names = BTreeSet::new();
        let mut edge_keys = BTreeSet::new();
        for row in rows {
            node_names.insert(row.from.clone());
            node_names.insert(row.to.clone());
            // A duplicate edge would need an aggregation rule (sum? keep
            // latest?) that hasn't been decided yet (see ROADMAP.md's
            // "duplicate edge" open decision), so reject rather than guess.
            if !edge_keys.insert((row.from.clone(), row.to.clone())) {
                return Err(format!(
                    "duplicate edge {} -> {} in one sample/gene group",
                    row.from, row.to
                ));
            }
        }

        // BTreeSet iterates in sorted order, giving each node a stable index.
        let nodes: Vec<String> = node_names.into_iter().collect();
        let node_index: BTreeMap<&str, usize> = nodes
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect();

        // Sort edges too, so edge order (and therefore incidence-matrix
        // columns) never depends on the order rows appeared in the file.
        let mut sorted_rows: Vec<&Junction> = rows.iter().collect();
        sorted_rows.sort_by(|left, right| (&left.from, &left.to).cmp(&(&right.from, &right.to)));

        let edges: Vec<Edge> = sorted_rows
            .iter()
            .map(|row| Edge {
                from: node_index[row.from.as_str()],
                to: node_index[row.to.as_str()],
                count: row.count,
            })
            .collect();
        // Normalize within the sample/gene group so signal values are
        // comparable across samples with different sequencing depth.
        let signal = edges.iter().map(|edge| edge.count / total).collect();

        Ok(Self {
            nodes,
            edges,
            signal,
        })
    }
}
