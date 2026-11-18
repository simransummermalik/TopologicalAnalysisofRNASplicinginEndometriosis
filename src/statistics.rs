//! Explicit-family statistics orchestration. Scientific choices are supplied
//! by the caller; no gene list, sidedness, or replication threshold is inferred.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use crate::GroupResult;
use crate::compare::{GeneComparisonSummary, IncludedSample, compare_gene, summarize_comparison};
use crate::fdr::{FdrFamily, GeneTest, benjamini_hochberg};
use crate::manifest::{Condition, SampleRecord};
use crate::permutation::{PermutationResult, Sidedness, permutation_test, set_permutation_test};

#[derive(Debug)]
pub struct StatsOptions {
    pub family_id: String,
    pub alternative: Sidedness,
    pub min_per_group: usize,
    pub rank_by: Option<crate::ranking::RankBy>,
}

#[derive(Debug)]
pub struct GeneStatistics {
    pub gene_id: String,
    pub status: &'static str,
    pub reason: String,
    pub permutation: Option<PermutationResult>,
    pub adjusted_p_value: Option<f64>,
}

/// One exact permutation test for the whole requested family (gene-set level).
/// One relabeling is applied to every gene at once, so correlation between
/// genes is preserved. Two-sided: the statistic is the mean over genes of the
/// absolute case-minus-control effect, so opposite changes do not cancel.
/// One-sided: the mean signed effect in the stated direction.
#[derive(Debug)]
pub struct SetTest {
    /// Tested genes whose included samples equal the common sample set.
    pub genes_used: usize,
    /// Tested genes left out because they lack one or more common samples.
    pub genes_left_out: usize,
    /// Descriptive per-sample mean cycle fraction over the genes used, in
    /// sample-ID order. For the two-sided test this is not the test statistic.
    pub samples: Vec<IncludedSample>,
    /// Descriptive direction: mean over genes of the signed effect, which
    /// equals the case-minus-control difference of the per-sample scores.
    pub signed_effect: Option<f64>,
    pub permutation: Option<PermutationResult>,
    pub reason: String,
}

#[derive(Debug)]
pub struct StatisticsRun {
    pub options: StatsOptions,
    pub family: FdrFamily,
    pub comparisons: Vec<GeneComparisonSummary>,
    pub genes: Vec<GeneStatistics>,
    pub annotations: std::collections::BTreeMap<String, crate::annotations::GeneAnnotation>,
    pub set_test: SetTest,
}

pub fn alternative_name(alternative: Sidedness) -> &'static str {
    match alternative {
        Sidedness::TwoSided => "two-sided",
        Sidedness::OneSidedCaseGreater => "case-greater",
        Sidedness::OneSidedCaseLess => "case-less",
    }
}

/// Name of the gene-set statistic used for an alternative.
pub fn set_statistic_name(alternative: Sidedness) -> &'static str {
    match alternative {
        Sidedness::TwoSided => "mean_absolute_effect",
        _ => "mean_signed_effect",
    }
}

pub fn parse_alternative(value: &str) -> Result<Sidedness, String> {
    match value {
        "two-sided" => Ok(Sidedness::TwoSided),
        "case-greater" => Ok(Sidedness::OneSidedCaseGreater),
        "case-less" => Ok(Sidedness::OneSidedCaseLess),
        _ => Err("--alternative must be two-sided, case-greater, or case-less; choose it before examining results".to_owned()),
    }
}

/// Accept a TSV containing one unique gene_id column, including richer
/// candidate-gene tables. Extra columns are preserved in the checksummed input.
pub fn parse_gene_family(path: &Path) -> Result<Vec<String>, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("could not read gene family {}: {error}", path.display()))?;
    parse_gene_family_text(&text)
        .map_err(|error| format!("invalid gene family {}: {error}", path.display()))
}

pub fn parse_gene_family_text(text: &str) -> Result<Vec<String>, String> {
    let mut lines = text.lines();
    let columns: Vec<_> = lines
        .next()
        .ok_or("empty file; add a gene_id header and at least one gene")?
        .split('\t')
        .collect();
    let indices: Vec<_> = columns
        .iter()
        .enumerate()
        .filter(|(_, value)| **value == "gene_id")
        .map(|(index, _)| index)
        .collect();
    if indices.len() != 1 {
        return Err("header must contain exactly one gene_id column".to_owned());
    }
    let mut genes = BTreeSet::new();
    for (offset, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != columns.len() {
            return Err(format!(
                "line {}: expected {} tab-separated fields",
                offset + 2,
                columns.len()
            ));
        }
        let gene = fields[indices[0]];
        if gene.trim().is_empty() || gene != gene.trim() {
            return Err(format!(
                "line {}: gene_id must be nonempty with no surrounding whitespace",
                offset + 2
            ));
        }
        if !genes.insert(gene.to_owned()) {
            return Err(format!(
                "line {}: duplicate gene_id '{gene}'; supply one row per requested gene",
                offset + 2
            ));
        }
    }
    if genes.is_empty() {
        return Err("no gene rows; supply at least one requested gene".to_owned());
    }
    Ok(genes.into_iter().collect())
}

/// Reject known repeated subjects when rich metadata provides subject_id.
/// Absence of that column does not establish independence: --design independent
/// remains a caller declaration, and confounders require scientific review.
pub fn validate_independent_metadata(path: &Path) -> Result<(), String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("could not read manifest {}: {error}", path.display()))?;
    let mut lines = text.lines();
    let columns: Vec<_> = lines.next().unwrap_or("").split('\t').collect();
    let subject_columns: Vec<_> = columns
        .iter()
        .enumerate()
        .filter(|(_, name)| **name == "subject_id")
        .map(|(index, _)| index)
        .collect();
    if subject_columns.len() > 1 {
        return Err("manifest has duplicate subject_id columns; keep exactly one".to_owned());
    }
    if let Some(&index) = subject_columns.first() {
        let mut seen = BTreeSet::new();
        for (offset, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let value = line.split('\t').nth(index).unwrap_or("").trim();
            if value.is_empty() || matches!(value, "NA" | "unknown" | "pending") {
                return Err(format!(
                    "manifest {} line {}: subject_id is missing; resolve subject identities before independent testing",
                    path.display(),
                    offset + 2
                ));
            }
            if !seen.insert(value) {
                return Err(format!(
                    "subject_id '{value}' appears more than once; --design independent does not support repeated subjects or paired samples"
                ));
            }
        }
    }
    Ok(())
}

pub fn test_gene_family(
    manifest: &[SampleRecord],
    results: &[GroupResult],
    gene_ids: &[String],
    options: StatsOptions,
) -> Result<StatisticsRun, String> {
    if options.min_per_group < 2 {
        return Err(
            "--min-per-group must be at least 2; one sample per group is descriptive only"
                .to_owned(),
        );
    }
    if gene_ids.is_empty() {
        return Err("the requested gene family is empty".to_owned());
    }
    let mut comparisons = Vec::new();
    let mut genes = Vec::new();
    for gene_id in gene_ids {
        let summary = summarize_comparison(compare_gene(manifest, results, gene_id));
        let cases = summary.case.as_ref().map_or(0, |group| group.n);
        let controls = summary.control.as_ref().map_or(0, |group| group.n);
        let mut row = GeneStatistics {
            gene_id: gene_id.clone(),
            status: "untested",
            reason: String::new(),
            permutation: None,
            adjusted_p_value: None,
        };
        if !results.iter().any(|result| result.gene_id == *gene_id) {
            row.reason = "no junction data for this requested gene".to_owned();
        } else if cases < options.min_per_group || controls < options.min_per_group {
            row.reason = format!(
                "requires at least {} eligible samples per group; case={cases}, control={controls}",
                options.min_per_group
            );
        } else {
            match permutation_test(&summary.comparison.included, options.alternative) {
                Ok(test) => {
                    row.status = "tested";
                    row.permutation = Some(test);
                }
                Err(reason) => row.reason = reason,
            }
        }
        comparisons.push(summary);
        genes.push(row);
    }
    let tests: Vec<_> = genes
        .iter()
        .map(|row| GeneTest {
            gene_id: row.gene_id.clone(),
            p_value: row.permutation.as_ref().map(|test| test.p_value),
        })
        .collect();
    let family = benjamini_hochberg(&options.family_id, &tests)?;
    for (row, adjusted) in genes.iter_mut().zip(&family.genes) {
        row.adjusted_p_value = adjusted.adjusted_p_value;
    }
    let set_test = test_gene_set(&comparisons, &genes, options.alternative);
    Ok(StatisticsRun {
        options,
        family,
        comparisons,
        genes,
        annotations: std::collections::BTreeMap::new(),
        set_test,
    })
}

/// Gene-set level test over the genes that passed the per-gene testability
/// rule. A gene enters only when its included samples equal the union of
/// included samples across tested genes, so every relabeling is applied to the
/// same samples for every gene and no value is imputed.
pub fn test_gene_set(
    comparisons: &[GeneComparisonSummary],
    genes: &[GeneStatistics],
    alternative: Sidedness,
) -> SetTest {
    let tested: Vec<&GeneComparisonSummary> = comparisons
        .iter()
        .zip(genes)
        .filter(|(_, row)| row.permutation.is_some())
        .map(|(summary, _)| summary)
        .collect();
    let mut conditions = BTreeMap::new();
    for summary in &tested {
        for sample in &summary.comparison.included {
            conditions.insert(sample.sample_id.clone(), sample.condition);
        }
    }
    let order: Vec<&String> = conditions.keys().collect();
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    for summary in &tested {
        let values: BTreeMap<&str, f64> = summary
            .comparison
            .included
            .iter()
            .map(|sample| (sample.sample_id.as_str(), sample.cycle_fraction))
            .collect();
        if values.len() == order.len() {
            matrix.push(order.iter().map(|id| values[id.as_str()]).collect());
        }
    }
    let mut result = SetTest {
        genes_used: matrix.len(),
        genes_left_out: tested.len() - matrix.len(),
        samples: Vec::new(),
        signed_effect: None,
        permutation: None,
        reason: String::new(),
    };
    if matrix.is_empty() {
        result.reason = if tested.is_empty() {
            "no requested gene was testable".to_owned()
        } else {
            "no tested gene has every common sample".to_owned()
        };
        return result;
    }
    result.samples = order
        .iter()
        .enumerate()
        .map(|(column, id)| IncludedSample {
            sample_id: (*id).clone(),
            condition: conditions[*id],
            cycle_fraction: matrix.iter().map(|gene| gene[column]).sum::<f64>()
                / matrix.len() as f64,
        })
        .collect();
    let mean = |condition: Condition| {
        let scores: Vec<f64> = result
            .samples
            .iter()
            .filter(|sample| sample.condition == condition)
            .map(|sample| sample.cycle_fraction)
            .collect();
        scores.iter().sum::<f64>() / scores.len() as f64
    };
    let is_case: Vec<bool> = result
        .samples
        .iter()
        .map(|sample| sample.condition == Condition::Case)
        .collect();
    match set_permutation_test(&matrix, &is_case, alternative) {
        Ok(test) => {
            result.signed_effect = Some(mean(Condition::Case) - mean(Condition::Control));
            result.permutation = Some(test);
        }
        Err(reason) => result.reason = reason,
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{analyze_junctions, manifest::parse_manifest_text, parser::parse_tsv_text};

    #[test]
    fn family_parser_accepts_extra_columns_and_rejects_bad_rows() {
        assert_eq!(
            parse_gene_family_text("gene_symbol\tgene_id\nB\tG2\nA\tG1\n").unwrap(),
            ["G1", "G2"]
        );
        for invalid in [
            "",
            "gene\nG1\n",
            "gene_id\tgene_id\nA\tB\n",
            "gene_id\n",
            "gene_id\nG1\nG1\n",
            "gene_id\n G1\n",
            "gene_id\tx\nG1\n",
        ] {
            assert!(parse_gene_family_text(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn explicit_thresholds_keep_requested_missing_genes_untested() {
        let manifest = parse_manifest_text("sample_id\tcondition\na\tcase\nb\tcontrol\n").unwrap();
        let results = analyze_junctions(
            parse_tsv_text(
                "sample_id\tgene_id\tfrom\tto\tcount\na\tG\tE1\tE2\t1\nb\tG\tE1\tE2\t1\n",
            )
            .unwrap(),
        )
        .unwrap();
        let genes = vec!["G".to_owned(), "absent".to_owned()];
        let options = StatsOptions {
            family_id: "toy".to_owned(),
            alternative: Sidedness::TwoSided,
            min_per_group: 2,
            rank_by: None,
        };
        let run = test_gene_family(&manifest, &results, &genes, options).unwrap();
        assert_eq!(run.family.tested_count, 0);
        assert!(run.genes[0].reason.contains("case=1, control=1"));
        assert!(run.genes[1].reason.contains("no junction data"));
        assert!(
            run.genes
                .iter()
                .all(|row| row.permutation.is_none() && row.adjusted_p_value.is_none())
        );
        assert!(run.comparisons[0].mean_difference.is_some());
        assert_eq!(run.set_test.genes_used, 0);
        assert!(run.set_test.permutation.is_none());
        assert_eq!(run.set_test.reason, "no requested gene was testable");
    }

    #[test]
    fn set_test_averages_complete_genes_and_matches_mean_gene_effect() {
        // Cases a,b; controls c,d. CYC has a cycle (fraction 1) only in a and b.
        // LINE is acyclic everywhere. PART lacks sample d and must be left out.
        let mut manifest = "sample_id\tcondition\n".to_owned();
        let mut junctions = "sample_id\tgene_id\tfrom\tto\tcount\n".to_owned();
        for (sample, condition) in [
            ("a", "case"),
            ("b", "case"),
            ("c", "control"),
            ("d", "control"),
        ] {
            manifest.push_str(&format!("{sample}\t{condition}\n"));
            junctions.push_str(&format!(
                "{sample}\tCYC\tE1\tE2\t5\n{sample}\tCYC\tE2\tE3\t5\n{sample}\tLINE\tE1\tE2\t5\n"
            ));
            if condition == "case" {
                junctions.push_str(&format!("{sample}\tCYC\tE3\tE1\t5\n"));
            }
            if sample != "d" {
                junctions.push_str(&format!("{sample}\tPART\tE1\tE2\t5\n"));
            }
        }
        let manifest = parse_manifest_text(&manifest).unwrap();
        let results = analyze_junctions(parse_tsv_text(&junctions).unwrap()).unwrap();
        let genes: Vec<String> = ["CYC", "LINE", "PART"].map(str::to_owned).to_vec();
        // Two samples per group is enough to show the arithmetic; rows are
        // built directly so the set calculation is tested on its own.
        let comparisons: Vec<_> = genes
            .iter()
            .map(|gene| summarize_comparison(compare_gene(&manifest, &results, gene)))
            .collect();
        let rows: Vec<_> = comparisons
            .iter()
            .map(|summary| GeneStatistics {
                gene_id: summary.comparison.gene_id.clone(),
                status: "tested",
                reason: String::new(),
                permutation: permutation_test(&summary.comparison.included, Sidedness::TwoSided)
                    .ok(),
                adjusted_p_value: None,
            })
            .collect();
        let set = test_gene_set(&comparisons, &rows, Sidedness::TwoSided);
        assert_eq!((set.genes_used, set.genes_left_out), (2, 1));
        let scores: Vec<_> = set
            .samples
            .iter()
            .map(|s| (s.sample_id.as_str(), s.cycle_fraction))
            .collect();
        assert_eq!(scores.len(), 4);
        for (id, score) in scores {
            let expected = if id == "a" || id == "b" { 0.5 } else { 0.0 };
            assert!((score - expected).abs() < 1e-9, "{id}: {score}");
        }
        let test = set.permutation.unwrap();
        // The signed effect equals the mean of the complete genes' effects
        // (1 and 0); the two-sided statistic is the mean of their magnitudes.
        let mean_gene_effect = (comparisons[0].mean_difference.unwrap()
            + comparisons[1].mean_difference.unwrap())
            / 2.0;
        assert!((set.signed_effect.unwrap() - mean_gene_effect).abs() < 1e-9);
        assert!((test.observed_effect - 0.5).abs() < 1e-9);
        // C(4,2) = 6 labelings; the observed split and its complement are extreme.
        assert_eq!((test.permutations_considered, test.extreme_count), (6, 2));
    }
}
