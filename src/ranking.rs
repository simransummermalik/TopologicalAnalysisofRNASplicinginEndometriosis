//! Optional ordering of tested results. Ordering is not a significance decision.
use crate::statistics::StatisticsRun;

#[derive(Clone, Copy, Debug)]
pub enum RankBy {
    AbsoluteEffect,
    AdjustedP,
}

impl RankBy {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "absolute-effect" => Ok(Self::AbsoluteEffect),
            "adjusted-p" => Ok(Self::AdjustedP),
            _ => Err("--rank-by must be absolute-effect or adjusted-p".to_owned()),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::AbsoluteEffect => "absolute-effect",
            Self::AdjustedP => "adjusted-p",
        }
    }
}

/// Tested genes first; effect magnitude breaks adjusted-p ties, then gene ID.
/// Untested genes are retained alphabetically and must not receive a rank.
pub fn ordered_genes(run: &StatisticsRun, rule: RankBy) -> Vec<usize> {
    let mut indices: Vec<_> = (0..run.genes.len()).collect();
    indices.sort_by(|&a, &b| {
        let left = &run.genes[a];
        let right = &run.genes[b];
        let tested = |index: usize| run.genes[index].permutation.is_some();
        let effect = |index: usize| run.comparisons[index].mean_difference.unwrap_or(0.0).abs();
        tested(b).cmp(&tested(a)).then_with(|| {
            if !tested(a) {
                return left.gene_id.cmp(&right.gene_id);
            }
            let primary = match rule {
                RankBy::AbsoluteEffect => effect(b).total_cmp(&effect(a)),
                RankBy::AdjustedP => left
                    .adjusted_p_value
                    .unwrap()
                    .total_cmp(&right.adjusted_p_value.unwrap()),
            };
            primary
                .then_with(|| effect(b).total_cmp(&effect(a)))
                .then_with(|| left.gene_id.cmp(&right.gene_id))
        })
    });
    indices
}
