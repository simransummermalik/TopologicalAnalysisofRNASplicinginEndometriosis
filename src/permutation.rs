//! Exact permutation testing for one gene's case/control comparison, per
//! `docs/statistics_spec.md` S5. Takes the already-eligible samples from
//! `compare::GeneComparison` (inclusion/exclusion is decided there, not
//! here) and tests whether the observed case-minus-control effect is
//! unusual among every possible relabeling of the same group sizes.
//!
//! The permutation space is enumerated exactly, not sampled: S5's own
//! reasoning is that exhaustive enumeration is feasible for this project's
//! sample sizes and removes any question of whether enough permutations
//! were run. `MAX_PERMUTATIONS` makes that assumption's limit explicit
//! rather than letting a larger future cohort silently hang or exhaust
//! memory.

use std::collections::BTreeSet;

use crate::compare::IncludedSample;
use crate::hodge::TOLERANCE;
use crate::manifest::Condition;

/// Exhaustive enumeration stops being the right tool above this many
/// permutations. S5 states that a cohort this large would need the spec
/// revised to add a sampling procedure and a saved seed; this is where that
/// revision would need to happen, rather than letting enumeration run
/// unbounded.
const MAX_PERMUTATIONS: u128 = 1_000_000;

/// Which direction(s) count as "at least as extreme" as the observed effect.
/// `docs/statistics_spec.md` S5 leaves the choice between these open as a
/// scientific decision, not a mechanical one — this type exists so the code
/// never silently picks a side of that decision on its own.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sidedness {
    TwoSided,
    OneSidedCaseGreater,
    OneSidedCaseLess,
}

/// The result of an exact permutation test for one gene.
#[derive(Clone, Debug, PartialEq)]
pub struct PermutationResult {
    /// Case mean minus control mean under the real, observed labeling.
    pub observed_effect: f64,
    /// Total permutations in the exhaustive space, `C(n, k)` for `n`
    /// eligible samples and `k` of them actually labeled case. This is
    /// computed per call, not assumed to be the full-cohort `70`, since a
    /// specific gene may have fewer eligible samples than the full cohort.
    pub permutations_considered: usize,
    /// How many of those permutations were at least as extreme as the
    /// observed effect, counting the observed labeling itself.
    pub extreme_count: usize,
    pub p_value: f64,
}

/// Runs the exact permutation test described in S5 on one gene's already
/// eligible samples (typically `GeneComparison::included`, or the case/
/// control-only subset of it). Returns an error if either group is empty
/// (no comparison is possible) or if the exhaustive permutation space would
/// exceed `MAX_PERMUTATIONS`.
pub fn permutation_test(
    samples: &[IncludedSample],
    sidedness: Sidedness,
) -> Result<PermutationResult, String> {
    let n = samples.len();
    let values: Vec<f64> = samples.iter().map(|sample| sample.cycle_fraction).collect();
    let case_indices: Vec<usize> = samples
        .iter()
        .enumerate()
        .filter(|(_, sample)| sample.condition == Condition::Case)
        .map(|(index, _)| index)
        .collect();
    let k = case_indices.len();

    if k == 0 || k == n {
        return Err(
            "permutation test requires at least one eligible sample in each group".to_owned(),
        );
    }

    let observed_effect = group_difference(&values, &case_indices, n);
    let permutations_considered = binomial_coefficient(n, k)?;

    let tolerance = TOLERANCE * observed_effect.abs().max(1.0);
    let mut extreme_count = 0;
    for case_subset in combinations(n, k) {
        let effect = group_difference(&values, &case_subset, n);
        if is_at_least_as_extreme(effect, observed_effect, sidedness, tolerance) {
            extreme_count += 1;
        }
    }
    debug_assert!(extreme_count >= 1, "the observed labeling must count");

    Ok(PermutationResult {
        observed_effect,
        permutations_considered,
        extreme_count,
        p_value: extreme_count as f64 / permutations_considered as f64,
    })
}

/// Case mean minus control mean, where `case_indices` names which positions
/// in `values` (length `n`) are labeled case for this particular relabeling.
fn group_difference(values: &[f64], case_indices: &[usize], n: usize) -> f64 {
    let case_set: BTreeSet<usize> = case_indices.iter().copied().collect();
    let case_mean = mean_of(case_indices.iter().map(|&index| values[index]));
    let control_mean = mean_of(
        (0..n)
            .filter(|index| !case_set.contains(index))
            .map(|index| values[index]),
    );
    case_mean - control_mean
}

fn mean_of(values: impl Iterator<Item = f64>) -> f64 {
    let (sum, count) = values.fold((0.0, 0usize), |(sum, count), value| {
        (sum + value, count + 1)
    });
    sum / count as f64
}

fn is_at_least_as_extreme(
    effect: f64,
    observed: f64,
    sidedness: Sidedness,
    tolerance: f64,
) -> bool {
    match sidedness {
        Sidedness::TwoSided => effect.abs() >= observed.abs() - tolerance,
        Sidedness::OneSidedCaseGreater => effect >= observed - tolerance,
        Sidedness::OneSidedCaseLess => effect <= observed + tolerance,
    }
}

/// Computes `C(n, k)` using the standard incremental form, which stays
/// integral at every step rather than risking truncation from dividing too
/// early. Fails loudly instead of enumerating past `MAX_PERMUTATIONS`.
fn binomial_coefficient(n: usize, k: usize) -> Result<usize, String> {
    let k = k.min(n - k);
    let mut result: u128 = 1;
    for i in 1..=k {
        result = result * (n - k + i) as u128 / i as u128;
        if result > MAX_PERMUTATIONS {
            return Err(format!(
                "permutation space C({n}, {k}) exceeds the {MAX_PERMUTATIONS}-permutation \
                 enumeration limit; exhaustive enumeration is no longer feasible for this \
                 cohort size, and docs/statistics_spec.md S5 must be revised to add a sampling \
                 procedure and a saved seed before this gene can be tested"
            ));
        }
    }
    Ok(result as usize)
}

/// Every distinct way to choose `k` indices out of `0..n`, as sorted index
/// lists. Hand-rolled rather than pulled from a crate, matching this
/// project's dependency-free convention (see `src/matrix.rs`).
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(k);
    combinations_from(0, n, k, &mut current, &mut result);
    result
}

fn combinations_from(
    start: usize,
    n: usize,
    k: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == k {
        result.push(current.clone());
        return;
    }
    for index in start..n {
        current.push(index);
        combinations_from(index + 1, n, k, current, result);
        current.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(sample_id: &str, condition: Condition, cycle_fraction: f64) -> IncludedSample {
        IncludedSample {
            sample_id: sample_id.to_owned(),
            condition,
            cycle_fraction,
        }
    }

    #[test]
    fn combinations_of_four_choose_two_match_hand_enumeration() {
        let mut found = combinations(4, 2);
        found.sort();
        let mut expected = vec![
            vec![0, 1],
            vec![0, 2],
            vec![0, 3],
            vec![1, 2],
            vec![1, 3],
            vec![2, 3],
        ];
        expected.sort();
        assert_eq!(found, expected);
    }

    #[test]
    fn binomial_coefficient_matches_known_small_values() {
        assert_eq!(binomial_coefficient(8, 4).unwrap(), 70);
        assert_eq!(binomial_coefficient(4, 2).unwrap(), 6);
        assert_eq!(binomial_coefficient(5, 0).unwrap(), 1);
    }

    #[test]
    fn binomial_coefficient_refuses_an_infeasible_enumeration() {
        // C(50, 25) is on the order of 1.26e13, far past MAX_PERMUTATIONS.
        let error = binomial_coefficient(50, 25).unwrap_err();
        assert!(error.contains("exceeds"));
        assert!(error.contains("statistics_spec.md"));
    }

    /// A fully hand-enumerable case: 2 case, 2 control, so all C(4,2) = 6
    /// relabelings can be checked by hand rather than trusted from the code.
    /// Values: case = [0.9, 0.8], control = [0.1, 0.2].
    /// Observed effect = mean(0.9, 0.8) - mean(0.1, 0.2) = 0.85 - 0.15 = 0.70.
    /// Hand-computed effect for every one of the 6 relabelings of {0,1,2,3}
    /// into two groups of 2 (effect = mean(subset) - mean(complement)):
    ///   {0,1}: 0.70   {0,2}: 0.00   {0,3}: 0.10
    ///   {1,2}: -0.10  {1,3}: 0.00   {2,3}: -0.70
    /// Two of six (±0.70) are at least as extreme two-sided; one of six
    /// (0.70 itself) is at least as extreme one-sided (case greater).
    fn small_known_difference_fixture() -> Vec<IncludedSample> {
        vec![
            sample("case_a", Condition::Case, 0.9),
            sample("case_b", Condition::Case, 0.8),
            sample("control_a", Condition::Control, 0.1),
            sample("control_b", Condition::Control, 0.2),
        ]
    }

    #[test]
    fn two_sided_p_value_matches_hand_enumeration() {
        let result =
            permutation_test(&small_known_difference_fixture(), Sidedness::TwoSided).unwrap();
        assert!((result.observed_effect - 0.70).abs() < 1e-9);
        assert_eq!(result.permutations_considered, 6);
        assert_eq!(result.extreme_count, 2);
        assert!((result.p_value - 1.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn one_sided_p_value_matches_hand_enumeration() {
        let result = permutation_test(
            &small_known_difference_fixture(),
            Sidedness::OneSidedCaseGreater,
        )
        .unwrap();
        assert_eq!(result.extreme_count, 1);
        assert!((result.p_value - 1.0 / 6.0).abs() < 1e-9);
    }

    /// The full eight-sample cohort shape from `docs/statistics_spec.md`,
    /// with case values strictly greater than every control value so the
    /// observed labeling is the unique most extreme split in both
    /// directions. This independently reproduces the spec's documented
    /// p-value floor (1/70 one-sided, 2/70 two-sided) rather than just
    /// asserting a number copied from the doc.
    fn full_cohort_known_difference_fixture() -> Vec<IncludedSample> {
        vec![
            sample("case_1", Condition::Case, 0.9),
            sample("case_2", Condition::Case, 0.8),
            sample("case_3", Condition::Case, 0.7),
            sample("case_4", Condition::Case, 0.6),
            sample("control_1", Condition::Control, 0.4),
            sample("control_2", Condition::Control, 0.3),
            sample("control_3", Condition::Control, 0.2),
            sample("control_4", Condition::Control, 0.1),
        ]
    }

    #[test]
    fn full_cohort_reaches_the_documented_p_value_floor() {
        let one_sided = permutation_test(
            &full_cohort_known_difference_fixture(),
            Sidedness::OneSidedCaseGreater,
        )
        .unwrap();
        assert_eq!(one_sided.permutations_considered, 70);
        assert_eq!(one_sided.extreme_count, 1);
        assert!((one_sided.p_value - 1.0 / 70.0).abs() < 1e-9);

        let two_sided =
            permutation_test(&full_cohort_known_difference_fixture(), Sidedness::TwoSided).unwrap();
        assert_eq!(two_sided.extreme_count, 2);
        assert!((two_sided.p_value - 2.0 / 70.0).abs() < 1e-9);
    }

    #[test]
    fn no_difference_fixture_gives_p_value_of_one() {
        let samples = vec![
            sample("case_1", Condition::Case, 0.5),
            sample("case_2", Condition::Case, 0.5),
            sample("control_1", Condition::Control, 0.5),
            sample("control_2", Condition::Control, 0.5),
        ];
        let result = permutation_test(&samples, Sidedness::TwoSided).unwrap();
        assert!((result.observed_effect).abs() < 1e-9);
        assert_eq!(result.extreme_count, result.permutations_considered);
        assert!((result.p_value - 1.0).abs() < 1e-9);
    }

    #[test]
    fn p_value_does_not_depend_on_input_order() {
        let mut shuffled = full_cohort_known_difference_fixture();
        shuffled.reverse();
        shuffled.swap(0, 3);

        let ordered =
            permutation_test(&full_cohort_known_difference_fixture(), Sidedness::TwoSided).unwrap();
        let reordered = permutation_test(&shuffled, Sidedness::TwoSided).unwrap();

        assert!((ordered.p_value - reordered.p_value).abs() < 1e-9);
        assert_eq!(
            ordered.permutations_considered,
            reordered.permutations_considered
        );
    }

    #[test]
    fn rejects_a_comparison_with_an_empty_group() {
        let samples = vec![
            sample("case_1", Condition::Case, 0.5),
            sample("case_2", Condition::Case, 0.6),
        ];
        let error = permutation_test(&samples, Sidedness::TwoSided).unwrap_err();
        assert!(error.contains("at least one eligible sample in each group"));
    }
}
