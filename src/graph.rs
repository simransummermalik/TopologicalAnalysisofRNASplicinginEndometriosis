use std::collections::{BTreeMap, BTreeSet};

use crate::parser::Junction;

#[derive(Clone, Debug)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub count: f64,
}

#[derive(Debug)]
pub struct SpliceGraph {
    pub nodes: Vec<String>,
    pub edges: Vec<Edge>,
    pub signal: Vec<f64>,
}

impl SpliceGraph {
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
            if !edge_keys.insert((row.from.clone(), row.to.clone())) {
                return Err(format!(
                    "duplicate edge {} -> {} in one sample/gene group",
                    row.from, row.to
                ));
            }
        }

        let nodes: Vec<String> = node_names.into_iter().collect();
        let node_index: BTreeMap<&str, usize> = nodes
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect();

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
        let signal = edges.iter().map(|edge| edge.count / total).collect();

        Ok(Self {
            nodes,
            edges,
            signal,
        })
    }
}
