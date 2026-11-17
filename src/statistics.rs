//! Explicit-family statistics orchestration. Scientific choices are supplied
//! by the caller; no gene list, sidedness, or replication threshold is inferred.

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

use crate::GroupResult;
use crate::compare::{GeneComparisonSummary, compare_gene, summarize_comparison};
use crate::fdr::{FdrFamily, GeneTest, benjamini_hochberg};
use crate::manifest::SampleRecord;
use crate::permutation::{PermutationResult, Sidedness, permutation_test};

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

#[derive(Debug)]
pub struct StatisticsRun {
    pub options: StatsOptions,
    pub family: FdrFamily,
    pub comparisons: Vec<GeneComparisonSummary>,
    pub genes: Vec<GeneStatistics>,
}

pub fn alternative_name(alternative: Sidedness) -> &'static str {
    match alternative {
        Sidedness::TwoSided => "two-sided",
        Sidedness::OneSidedCaseGreater => "case-greater",
        Sidedness::OneSidedCaseLess => "case-less",
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
    Ok(StatisticsRun {
        options,
        family,
        comparisons,
        genes,
    })
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
    }
}
