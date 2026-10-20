//! Joins a parsed sample manifest (`manifest.rs`) against the per-sample
//! decomposition results (`analyze_junctions` in `lib.rs`) for one gene,
//! producing a sample-preserving case/control comparison. This is the R10
//! layer: it adds no new math and does not touch `graph.rs`/`hodge.rs` —
//! it only decides which already-computed samples are eligible for one
//! gene's comparison and reports why any sample was left out.
use std::collections::{BTreeMap, BTreeSet};

use crate::GroupResult;
use crate::manifest::{Condition, SampleRecord};

/// One sample eligible for this gene's comparison.
#[derive(Debug)]
pub struct IncludedSample {
    pub sample_id: String,
    pub condition: Condition,
    pub cycle_fraction: f64,
}

/// One sample left out of this gene's comparison, with a stated reason.
#[derive(Debug)]
pub struct ExcludedSample {
    pub sample_id: String,
    pub reason: String,
}

/// The full result of comparing one gene across every sample the manifest
/// and the junction data agree on. Every manifest sample_id and every
/// junction-data sample_id for this gene ends up in exactly one of
/// `included` or `excluded` — nothing is silently dropped.
#[derive(Debug)]
pub struct GeneComparison {
    pub gene_id: String,
    pub included: Vec<IncludedSample>,
    pub excluded: Vec<ExcludedSample>,
}

/// Compares one gene across samples, preserving sample identity throughout.
///
/// A sample is included only if it appears in both the manifest and the
/// junction-derived results for `gene_id`, and its decomposition's checks
/// passed. Every other case is an explicit exclusion:
/// - present in the results but missing from the manifest,
/// - present in the manifest but with no junction data for this gene
///   (never conflated with a measured zero), or
/// - present in both but with a failed decomposition check.
pub fn compare_gene(
    manifest: &[SampleRecord],
    results: &[GroupResult],
    gene_id: &str,
) -> GeneComparison {
    let manifest_by_sample: BTreeMap<&str, &SampleRecord> = manifest
        .iter()
        .map(|record| (record.sample_id.as_str(), record))
        .collect();

    let gene_results_by_sample: BTreeMap<&str, &GroupResult> = results
        .iter()
        .filter(|result| result.gene_id == gene_id)
        .map(|result| (result.sample_id.as_str(), result))
        .collect();

    let all_sample_ids: BTreeSet<&str> = manifest_by_sample
        .keys()
        .chain(gene_results_by_sample.keys())
        .copied()
        .collect();

    let mut included = Vec::new();
    let mut excluded = Vec::new();

    for sample_id in all_sample_ids {
        match (
            manifest_by_sample.get(sample_id),
            gene_results_by_sample.get(sample_id),
        ) {
            (Some(record), Some(result)) if result.decomposition.checks.passed => {
                included.push(IncludedSample {
                    sample_id: sample_id.to_owned(),
                    condition: record.condition,
                    cycle_fraction: result.decomposition.cycle_fraction,
                });
            }
            (Some(_), Some(_)) => excluded.push(ExcludedSample {
                sample_id: sample_id.to_owned(),
                reason: format!("decomposition checks failed for gene '{gene_id}'"),
            }),
            (None, Some(_)) => excluded.push(ExcludedSample {
                sample_id: sample_id.to_owned(),
                reason: "no manifest entry for this sample_id".to_owned(),
            }),
            (Some(_), None) => excluded.push(ExcludedSample {
                sample_id: sample_id.to_owned(),
                reason: format!("no junction data for gene '{gene_id}'"),
            }),
            (None, None) => unreachable!("sample_id came from one of the two source maps"),
        }
    }

    GeneComparison {
        gene_id: gene_id.to_owned(),
        included,
        excluded,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Edge, SpliceGraph};
    use crate::hodge::{Checks, Decomposition};
    use crate::manifest::parse_manifest_text;
    use crate::parser::parse_tsv_text;

    fn sample(sample_id: &str, condition: Condition) -> SampleRecord {
        SampleRecord {
            sample_id: sample_id.to_owned(),
            condition,
        }
    }

    /// A minimal, hand-built result standing in for whatever `analyze_junctions`
    /// would produce, with `checks.passed` set directly — this lets the
    /// failed-checks exclusion be tested without needing a real fixture that
    /// happens to fail its numerical identities.
    fn fake_result(
        sample_id: &str,
        gene_id: &str,
        cycle_fraction: f64,
        passed: bool,
    ) -> GroupResult {
        GroupResult {
            sample_id: sample_id.to_owned(),
            gene_id: gene_id.to_owned(),
            graph: SpliceGraph {
                nodes: vec!["E1".to_owned(), "E2".to_owned()],
                edges: vec![Edge {
                    from: 0,
                    to: 1,
                    count: 1.0,
                }],
                signal: vec![1.0],
            },
            decomposition: Decomposition {
                gradient: vec![1.0],
                cycle: vec![0.0],
                gradient_energy: 1.0,
                cycle_energy: 0.0,
                signal_energy: 1.0,
                cycle_fraction,
                checks: Checks {
                    reconstruction_error: 0.0,
                    kernel_error: 0.0,
                    orthogonality_error: 0.0,
                    passed,
                },
            },
        }
    }

    #[test]
    fn includes_samples_present_in_both_sources() {
        let manifest = vec![
            sample("case1", Condition::Case),
            sample("control1", Condition::Control),
        ];
        let results = vec![
            fake_result("case1", "GENE_A", 0.8, true),
            fake_result("control1", "GENE_A", 0.1, true),
        ];

        let comparison = compare_gene(&manifest, &results, "GENE_A");

        assert!(comparison.excluded.is_empty());
        assert_eq!(comparison.included.len(), 2);
        assert_eq!(comparison.included[0].sample_id, "case1");
        assert_eq!(comparison.included[0].condition, Condition::Case);
        assert!((comparison.included[0].cycle_fraction - 0.8).abs() < 1e-12);
    }

    #[test]
    fn excludes_a_manifest_sample_with_no_data_for_the_gene() {
        let manifest = vec![sample("case1", Condition::Case)];
        let results = vec![fake_result("case1", "GENE_B", 0.5, true)];

        let comparison = compare_gene(&manifest, &results, "GENE_A");

        assert!(comparison.included.is_empty());
        assert_eq!(comparison.excluded.len(), 1);
        assert_eq!(comparison.excluded[0].sample_id, "case1");
        assert!(comparison.excluded[0].reason.contains("no junction data"));
    }

    #[test]
    fn excludes_a_result_sample_with_no_manifest_entry() {
        let manifest: Vec<SampleRecord> = Vec::new();
        let results = vec![fake_result("mystery_sample", "GENE_A", 0.5, true)];

        let comparison = compare_gene(&manifest, &results, "GENE_A");

        assert!(comparison.included.is_empty());
        assert_eq!(comparison.excluded.len(), 1);
        assert!(comparison.excluded[0].reason.contains("no manifest entry"));
    }

    #[test]
    fn excludes_a_sample_whose_decomposition_checks_failed() {
        let manifest = vec![sample("case1", Condition::Case)];
        let results = vec![fake_result("case1", "GENE_A", 0.5, false)];

        let comparison = compare_gene(&manifest, &results, "GENE_A");

        assert!(comparison.included.is_empty());
        assert_eq!(comparison.excluded.len(), 1);
        assert!(comparison.excluded[0].reason.contains("checks failed"));
    }

    #[test]
    fn end_to_end_with_real_parsing_for_one_gene_among_several() {
        let manifest = parse_manifest_text(
            "sample_id\tcondition\n\
             sample_linear\tcase\n\
             sample_cycle\tcontrol\n",
        )
        .expect("manifest should parse");

        let junctions = parse_tsv_text(
            "sample_id\tgene_id\tfrom\tto\tcount\n\
             sample_linear\tGENE_LINEAR\tE1\tE2\t100\n\
             sample_linear\tGENE_LINEAR\tE2\tE3\t80\n\
             sample_cycle\tGENE_LINEAR\tE1\tE2\t100\n\
             sample_cycle\tGENE_LINEAR\tE2\tE3\t100\n\
             sample_cycle\tGENE_OTHER\tE1\tE2\t50\n",
        )
        .expect("junctions should parse");
        let results = crate::analyze_junctions(junctions).expect("junctions should analyze");

        let comparison = compare_gene(&manifest, &results, "GENE_LINEAR");

        // Both samples have GENE_LINEAR data and a manifest entry; GENE_OTHER
        // is irrelevant to this comparison and must not appear anywhere.
        assert_eq!(comparison.included.len(), 2);
        assert!(comparison.excluded.is_empty());
    }
}
