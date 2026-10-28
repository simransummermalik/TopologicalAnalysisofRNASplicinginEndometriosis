//! Joins a parsed sample manifest (`manifest.rs`) against the per-sample
//! decomposition results (`analyze_junctions` in `lib.rs`) for one gene,
//! producing sample-preserving case/control comparisons. This layer adds no
//! inferential statistics and does not touch `graph.rs`/`hodge.rs`: it decides
//! which already-computed samples are eligible, reports every exclusion, and
//! can calculate transparent descriptive summaries for every observed gene.
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
    /// Known case/control label when the sample came from the manifest. This
    /// remains `None` for a junction-data sample that has no manifest row.
    pub condition: Option<Condition>,
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

/// Descriptive cycle-fraction values for one condition. These values keep the
/// biological sample as the unit: `n` is the number of eligible samples, not
/// the number of junctions.
#[derive(Clone, Debug, PartialEq)]
pub struct ConditionSummary {
    pub n: usize,
    pub mean: f64,
    pub median: f64,
    pub minimum: f64,
    pub maximum: f64,
}

/// Whether a gene has enough eligible data for a descriptive case/control
/// effect. Missing groups remain visible rather than being treated as zeros.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ComparisonStatus {
    Complete,
    NoEligibleSamples,
    MissingCase,
    MissingControl,
}

impl std::fmt::Display for ComparisonStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match self {
            Self::Complete => "complete",
            Self::NoEligibleSamples => "no_eligible_samples",
            Self::MissingCase => "missing_case",
            Self::MissingControl => "missing_control",
        };
        write!(f, "{text}")
    }
}

/// One gene's sample-preserving comparison plus descriptive group summaries.
/// `mean_difference` is defined only when both groups contain an eligible
/// sample and is always case mean minus control mean.
#[derive(Debug)]
pub struct GeneComparisonSummary {
    pub comparison: GeneComparison,
    pub case: Option<ConditionSummary>,
    pub control: Option<ConditionSummary>,
    pub mean_difference: Option<f64>,
    pub status: ComparisonStatus,
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
            (Some(record), Some(_)) => excluded.push(ExcludedSample {
                sample_id: sample_id.to_owned(),
                condition: Some(record.condition),
                reason: format!("decomposition checks failed for gene '{gene_id}'"),
            }),
            (None, Some(_)) => excluded.push(ExcludedSample {
                sample_id: sample_id.to_owned(),
                condition: None,
                reason: "no manifest entry for this sample_id".to_owned(),
            }),
            (Some(record), None) => excluded.push(ExcludedSample {
                sample_id: sample_id.to_owned(),
                condition: Some(record.condition),
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

/// Compares every gene represented in `results`, returning genes in stable
/// lexical order. This is descriptive by design: it does not assign p-values,
/// q-values, or biological importance before those rules are approved.
pub fn compare_all_genes(
    manifest: &[SampleRecord],
    results: &[GroupResult],
) -> Vec<GeneComparisonSummary> {
    let gene_ids: BTreeSet<&str> = results
        .iter()
        .map(|result| result.gene_id.as_str())
        .collect();

    gene_ids
        .into_iter()
        .map(|gene_id| summarize_comparison(compare_gene(manifest, results, gene_id)))
        .collect()
}

/// Adds descriptive case/control summaries to an existing per-gene
/// comparison without changing its inclusion or exclusion decisions.
pub fn summarize_comparison(comparison: GeneComparison) -> GeneComparisonSummary {
    let case_values: Vec<f64> = comparison
        .included
        .iter()
        .filter(|sample| sample.condition == Condition::Case)
        .map(|sample| sample.cycle_fraction)
        .collect();
    let control_values: Vec<f64> = comparison
        .included
        .iter()
        .filter(|sample| sample.condition == Condition::Control)
        .map(|sample| sample.cycle_fraction)
        .collect();

    let case = descriptive_summary(&case_values);
    let control = descriptive_summary(&control_values);
    let (status, mean_difference) = match (&case, &control) {
        (Some(case), Some(control)) => (ComparisonStatus::Complete, Some(case.mean - control.mean)),
        (None, Some(_)) => (ComparisonStatus::MissingCase, None),
        (Some(_), None) => (ComparisonStatus::MissingControl, None),
        (None, None) => (ComparisonStatus::NoEligibleSamples, None),
    };

    GeneComparisonSummary {
        comparison,
        case,
        control,
        mean_difference,
        status,
    }
}

fn descriptive_summary(values: &[f64]) -> Option<ConditionSummary> {
    if values.is_empty() {
        return None;
    }

    let mut ordered = values.to_vec();
    ordered.sort_by(f64::total_cmp);
    let n = ordered.len();
    let median = if n % 2 == 0 {
        (ordered[n / 2 - 1] + ordered[n / 2]) / 2.0
    } else {
        ordered[n / 2]
    };

    Some(ConditionSummary {
        n,
        mean: ordered.iter().sum::<f64>() / n as f64,
        median,
        minimum: ordered[0],
        maximum: ordered[n - 1],
    })
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

    #[test]
    fn summarizes_case_control_values_without_pooling_samples() {
        let comparison = GeneComparison {
            gene_id: "GENE_A".to_owned(),
            included: vec![
                IncludedSample {
                    sample_id: "case_1".to_owned(),
                    condition: Condition::Case,
                    cycle_fraction: 0.2,
                },
                IncludedSample {
                    sample_id: "case_2".to_owned(),
                    condition: Condition::Case,
                    cycle_fraction: 0.8,
                },
                IncludedSample {
                    sample_id: "control_1".to_owned(),
                    condition: Condition::Control,
                    cycle_fraction: 0.1,
                },
            ],
            excluded: Vec::new(),
        };

        let summary = summarize_comparison(comparison);

        assert_eq!(summary.status, ComparisonStatus::Complete);
        assert_eq!(summary.case.as_ref().unwrap().n, 2);
        assert!((summary.case.as_ref().unwrap().mean - 0.5).abs() < 1e-12);
        assert!((summary.case.as_ref().unwrap().median - 0.5).abs() < 1e-12);
        assert!((summary.mean_difference.unwrap() - 0.4).abs() < 1e-12);
    }

    #[test]
    fn missing_group_is_explicit_and_has_no_effect_value() {
        let comparison = GeneComparison {
            gene_id: "GENE_A".to_owned(),
            included: vec![IncludedSample {
                sample_id: "case_1".to_owned(),
                condition: Condition::Case,
                cycle_fraction: 0.2,
            }],
            excluded: Vec::new(),
        };

        let summary = summarize_comparison(comparison);

        assert_eq!(summary.status, ComparisonStatus::MissingControl);
        assert!(summary.control.is_none());
        assert!(summary.mean_difference.is_none());
    }

    #[test]
    fn all_gene_comparisons_are_sorted_and_keep_exclusions() {
        let manifest = vec![
            sample("case1", Condition::Case),
            sample("control1", Condition::Control),
        ];
        let results = vec![
            fake_result("case1", "GENE_Z", 0.7, true),
            fake_result("control1", "GENE_Z", 0.2, true),
            fake_result("case1", "GENE_A", 0.3, true),
        ];

        let summaries = compare_all_genes(&manifest, &results);

        assert_eq!(summaries.len(), 2);
        assert_eq!(summaries[0].comparison.gene_id, "GENE_A");
        assert_eq!(summaries[0].status, ComparisonStatus::MissingControl);
        assert_eq!(summaries[0].comparison.excluded.len(), 1);
        assert_eq!(summaries[1].comparison.gene_id, "GENE_Z");
        assert_eq!(summaries[1].status, ComparisonStatus::Complete);
    }
}
