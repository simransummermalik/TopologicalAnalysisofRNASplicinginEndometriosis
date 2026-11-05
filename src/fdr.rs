//! Benjamini–Hochberg adjusted p-values for an explicitly supplied family.
//! The caller chooses eligible tests before calling this library function.
//! It does not select genes, perform tests, or assign significance labels.

use std::collections::BTreeSet;

#[derive(Clone, Debug)]
pub struct GeneTest {
    pub gene_id: String,
    /// None means untested, not a measured p-value of zero.
    pub p_value: Option<f64>,
}

#[derive(Debug, PartialEq)]
pub struct AdjustedGene {
    pub gene_id: String,
    pub p_value: Option<f64>,
    /// BH-adjusted p-value; None for an untested gene. This is not a Storey
    /// q-value estimate or the probability that this individual gene is false.
    pub adjusted_p_value: Option<f64>,
}

#[derive(Debug, PartialEq)]
pub struct FdrFamily {
    pub family_id: String,
    pub method: &'static str,
    pub tested_count: usize,
    pub untested_count: usize,
    /// Original input order, with identifiers and raw values preserved.
    pub genes: Vec<AdjustedGene>,
}

/// Adjust the supplied valid tests, carrying untested rows for audit only.
/// The denominator is the number of Some(p) rows; no missing value is imputed.
/// For sorted p-values, a reverse cumulative minimum of m*p/rank enforces
/// monotonicity. Values are capped at one and returned in original order.
/// Scientific approval of the family and BH assumptions belongs upstream.
pub fn benjamini_hochberg(family_id: &str, tests: &[GeneTest]) -> Result<FdrFamily, String> {
    if family_id.trim().is_empty() {
        return Err("FDR family ID cannot be empty; name the set of eligible tests".to_owned());
    }
    let mut seen = BTreeSet::new();
    let mut tested = Vec::new();
    for (index, test) in tests.iter().enumerate() {
        if test.gene_id.trim().is_empty() {
            return Err(format!("FDR input row {} has an empty gene ID", index + 1));
        }
        if !seen.insert(test.gene_id.as_str()) {
            return Err(format!(
                "duplicate FDR gene ID '{}'; supply one test per gene in this family",
                test.gene_id
            ));
        }
        if let Some(p_value) = test.p_value {
            if !p_value.is_finite() || !(0.0..=1.0).contains(&p_value) {
                return Err(format!(
                    "gene '{}': p-value must be finite and between 0 and 1, found {p_value}",
                    test.gene_id
                ));
            }
            tested.push((index, p_value));
        }
    }
    tested.sort_by(|(left_index, left), (right_index, right)| {
        left.total_cmp(right)
            .then_with(|| left_index.cmp(right_index))
    });
    let tested_count = tested.len();
    let mut adjusted = vec![None; tests.len()];
    let mut ceiling = 1.0_f64;
    for (rank, &(index, p_value)) in tested.iter().enumerate().rev() {
        ceiling = ceiling.min(p_value * tested_count as f64 / (rank + 1) as f64);
        adjusted[index] = Some(ceiling);
    }
    Ok(FdrFamily {
        family_id: family_id.to_owned(),
        method: "Benjamini-Hochberg",
        tested_count,
        untested_count: tests.len() - tested_count,
        genes: tests
            .iter()
            .zip(adjusted)
            .map(|(test, adjusted_p_value)| AdjustedGene {
                gene_id: test.gene_id.clone(),
                p_value: test.p_value,
                adjusted_p_value,
            })
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(values: &[Option<f64>]) -> Vec<GeneTest> {
        values
            .iter()
            .enumerate()
            .map(|(index, &p_value)| GeneTest {
                gene_id: format!("G{index}"),
                p_value,
            })
            .collect()
    }

    fn close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
    }

    #[test]
    fn matches_hand_calculated_unsorted_family() {
        let input = rows(&[Some(0.01), Some(0.04), Some(0.03), Some(0.002)]);
        let result = benjamini_hochberg("toy", &input).unwrap();
        for (gene, expected) in result.genes.iter().zip([0.02, 0.04, 0.04, 0.008]) {
            close(gene.adjusted_p_value.unwrap(), expected);
        }
        assert_eq!(result.family_id, "toy");
        assert_eq!(result.method, "Benjamini-Hochberg");
        for (gene, original) in result.genes.iter().zip(&input) {
            assert_eq!(gene.gene_id, original.gene_id);
            assert_eq!(gene.p_value, original.p_value);
        }
    }

    #[test]
    fn reverse_minimum_prevents_nonmonotone_adjustments() {
        let result =
            benjamini_hochberg("toy", &rows(&[Some(0.01), Some(0.011), Some(0.012)])).unwrap();
        for gene in result.genes {
            close(gene.adjusted_p_value.unwrap(), 0.012);
        }
    }

    #[test]
    fn tied_values_are_equal_and_row_order_does_not_change_answers() {
        let mut input = rows(&[Some(0.02), Some(0.02), Some(0.8)]);
        let original = benjamini_hochberg("toy", &input).unwrap();
        close(original.genes[0].adjusted_p_value.unwrap(), 0.03);
        close(original.genes[1].adjusted_p_value.unwrap(), 0.03);
        input.reverse();
        let reversed = benjamini_hochberg("toy", &input).unwrap();
        for (left, right) in original.genes.iter().zip(reversed.genes.iter().rev()) {
            assert_eq!(left, right);
        }
    }

    #[test]
    fn untested_genes_remain_missing_and_do_not_expand_the_test_denominator() {
        let result = benjamini_hochberg("toy", &rows(&[None, Some(0.01), Some(0.04)])).unwrap();
        assert_eq!((result.tested_count, result.untested_count), (2, 1));
        assert_eq!(result.genes[0].p_value, None);
        assert_eq!(result.genes[0].adjusted_p_value, None);
        close(result.genes[1].adjusted_p_value.unwrap(), 0.02);
        close(result.genes[2].adjusted_p_value.unwrap(), 0.04);
    }

    #[test]
    fn rejects_invalid_p_values_with_the_gene_id() {
        for p_value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
            let error = benjamini_hochberg("toy", &rows(&[Some(p_value)])).unwrap_err();
            assert!(error.contains("G0"));
            assert!(error.contains("finite and between 0 and 1"));
        }
    }

    #[test]
    fn rejects_unnamed_families_and_empty_or_duplicate_genes() {
        assert!(benjamini_hochberg(" ", &[]).is_err());
        let mut input = rows(&[Some(0.1), Some(0.2)]);
        input[1].gene_id = input[0].gene_id.clone();
        assert!(
            benjamini_hochberg("toy", &input)
                .unwrap_err()
                .contains("duplicate")
        );
        input[1].gene_id = " ".to_owned();
        assert!(
            benjamini_hochberg("toy", &input)
                .unwrap_err()
                .contains("empty gene")
        );
    }

    #[test]
    fn handles_empty_untested_singleton_and_boundary_families() {
        let empty = benjamini_hochberg("toy", &[]).unwrap();
        assert_eq!(empty.tested_count, 0);
        let missing = benjamini_hochberg("toy", &rows(&[None, None])).unwrap();
        assert_eq!((missing.tested_count, missing.untested_count), (0, 2));
        for value in [0.0, 0.1, 1.0] {
            close(
                benjamini_hochberg("toy", &rows(&[Some(value)]))
                    .unwrap()
                    .genes[0]
                    .adjusted_p_value
                    .unwrap(),
                value,
            );
        }
        let result = benjamini_hochberg("toy", &rows(&[Some(0.0), Some(1.0), Some(0.9)])).unwrap();
        assert_eq!(
            result
                .genes
                .iter()
                .map(|gene| gene.adjusted_p_value.unwrap())
                .collect::<Vec<_>>(),
            [0.0, 1.0, 1.0]
        );
    }

    #[test]
    fn accepts_exact_permutation_results_without_changing_their_raw_values() {
        use crate::compare::IncludedSample;
        use crate::manifest::Condition;
        use crate::permutation::{Sidedness, permutation_test};
        let samples: Vec<_> = (0..8)
            .map(|index| IncludedSample {
                sample_id: format!("S{index}"),
                condition: if index < 4 {
                    Condition::Case
                } else {
                    Condition::Control
                },
                cycle_fraction: if index < 4 { 1.0 } else { 0.0 },
            })
            .collect();
        let p = permutation_test(&samples, Sidedness::TwoSided)
            .unwrap()
            .p_value;
        close(p, 2.0 / 70.0);
        let adjusted =
            benjamini_hochberg("synthetic-two-test-family", &rows(&[Some(p), Some(1.0)])).unwrap();
        assert_eq!(adjusted.genes[0].p_value, Some(p));
        close(adjusted.genes[0].adjusted_p_value.unwrap(), 4.0 / 70.0);
        close(adjusted.genes[1].adjusted_p_value.unwrap(), 1.0);
    }
}
