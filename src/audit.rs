//! Compact input and analysis audit used by `splice-girl validate`.
//!
//! Validation deliberately runs the same parser, graph builder, Hodge solver,
//! and mathematical checks as a real analysis, then reports counts and sample
//! mismatches without printing every edge component.

use std::collections::{BTreeMap, BTreeSet};

use crate::GroupResult;
use crate::manifest::{Condition, SampleRecord};
use crate::parser::Junction;

#[derive(Debug, PartialEq, Eq)]
pub struct ManifestAudit {
    pub manifest_samples: usize,
    pub case_samples: usize,
    pub control_samples: usize,
    pub junction_samples_without_metadata: Vec<String>,
    pub manifest_samples_without_junctions: Vec<String>,
}

#[derive(Debug, PartialEq)]
pub struct ValidationSummary {
    pub junction_rows: usize,
    pub samples: usize,
    pub genes: usize,
    pub sample_gene_groups: usize,
    pub total_read_support: f64,
    pub minimum_group_support: f64,
    pub maximum_group_support: f64,
    pub checks_passed: usize,
    pub checks_failed: usize,
    pub manifest: Option<ManifestAudit>,
}

pub fn summarize_validation(
    junctions: &[Junction],
    results: &[GroupResult],
    manifest: Option<&[SampleRecord]>,
) -> ValidationSummary {
    let junction_samples: BTreeSet<&str> = junctions
        .iter()
        .map(|junction| junction.sample_id.as_str())
        .collect();
    let genes: BTreeSet<&str> = junctions
        .iter()
        .map(|junction| junction.gene_id.as_str())
        .collect();
    let mut support_by_group: BTreeMap<(&str, &str), f64> = BTreeMap::new();
    for junction in junctions {
        *support_by_group
            .entry((junction.sample_id.as_str(), junction.gene_id.as_str()))
            .or_default() += junction.count;
    }
    let total_read_support = support_by_group.values().sum();
    let minimum_group_support = support_by_group
        .values()
        .copied()
        .min_by(f64::total_cmp)
        .unwrap_or(0.0);
    let maximum_group_support = support_by_group
        .values()
        .copied()
        .max_by(f64::total_cmp)
        .unwrap_or(0.0);
    let checks_passed = results
        .iter()
        .filter(|result| result.decomposition.checks.passed)
        .count();
    let checks_failed = results.len() - checks_passed;

    let manifest = manifest.map(|records| {
        let manifest_samples: BTreeSet<&str> = records
            .iter()
            .map(|record| record.sample_id.as_str())
            .collect();
        ManifestAudit {
            manifest_samples: records.len(),
            case_samples: records
                .iter()
                .filter(|record| record.condition == Condition::Case)
                .count(),
            control_samples: records
                .iter()
                .filter(|record| record.condition == Condition::Control)
                .count(),
            junction_samples_without_metadata: junction_samples
                .difference(&manifest_samples)
                .map(|sample| (*sample).to_owned())
                .collect(),
            manifest_samples_without_junctions: manifest_samples
                .difference(&junction_samples)
                .map(|sample| (*sample).to_owned())
                .collect(),
        }
    });

    ValidationSummary {
        junction_rows: junctions.len(),
        samples: junction_samples.len(),
        genes: genes.len(),
        sample_gene_groups: support_by_group.len(),
        total_read_support,
        minimum_group_support,
        maximum_group_support,
        checks_passed,
        checks_failed,
        manifest,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analyze_junctions;
    use crate::manifest::parse_manifest_text;
    use crate::parser::parse_tsv_text;

    fn sample_data() -> Vec<Junction> {
        parse_tsv_text(
            "sample_id\tgene_id\tfrom\tto\tcount\n\
             s1\tGENE_A\tE1\tE2\t10\n\
             s1\tGENE_A\tE2\tE3\t20\n\
             s2\tGENE_B\tE1\tE2\t5\n",
        )
        .unwrap()
    }

    #[test]
    fn summarizes_rows_samples_genes_support_and_checks() {
        let junctions = sample_data();
        let results = analyze_junctions(junctions.clone()).unwrap();

        let summary = summarize_validation(&junctions, &results, None);

        assert_eq!(summary.junction_rows, 3);
        assert_eq!(summary.samples, 2);
        assert_eq!(summary.genes, 2);
        assert_eq!(summary.sample_gene_groups, 2);
        assert_eq!(summary.total_read_support, 35.0);
        assert_eq!(summary.minimum_group_support, 5.0);
        assert_eq!(summary.maximum_group_support, 30.0);
        assert_eq!(summary.checks_passed, 2);
        assert_eq!(summary.checks_failed, 0);
        assert!(summary.manifest.is_none());
    }

    #[test]
    fn reports_manifest_mismatches_in_both_directions() {
        let junctions = sample_data();
        let results = analyze_junctions(junctions.clone()).unwrap();
        let manifest = parse_manifest_text(
            "sample_id\tcondition\n\
             s1\tcase\n\
             manifest_only\tcontrol\n",
        )
        .unwrap();

        let summary = summarize_validation(&junctions, &results, Some(&manifest));
        let audit = summary.manifest.unwrap();

        assert_eq!(audit.manifest_samples, 2);
        assert_eq!(audit.case_samples, 1);
        assert_eq!(audit.control_samples, 1);
        assert_eq!(audit.junction_samples_without_metadata, vec!["s2"]);
        assert_eq!(
            audit.manifest_samples_without_junctions,
            vec!["manifest_only"]
        );
    }
}
