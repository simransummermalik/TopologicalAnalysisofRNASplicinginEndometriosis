//! One-factor shared-graph sensitivity, using library calls rather than tools.
use super::*;
use crate::compare::{ComparisonStatus, compare_all_genes};
use crate::output::{write_compare_bundle, write_stats_bundle};
use crate::permutation::Sidedness;
use crate::statistics::{
    StatsOptions, alternative_name, test_gene_family, validate_independent_metadata,
};

/// Explicit statistical choices for optional per-scenario inference.
pub struct StatisticsSettings {
    pub family: String,
    pub alternative: Sidedness,
    pub min_per_group: usize,
}

pub struct SensitivityOptions {
    pub primary: Settings,
    pub edge_counts: Vec<u64>,
    pub sample_counts: Vec<usize>,
    pub gene_support: Vec<u64>,
    pub statistics: Option<StatisticsSettings>,
}

impl SensitivityOptions {
    /// Primary first, followed by unique one-factor changes in sorted order.
    pub fn scenarios(&self) -> Result<Vec<(String, Settings)>, String> {
        if self.edge_counts.is_empty()
            || self.sample_counts.is_empty()
            || self.gene_support.is_empty()
        {
            return Err("sensitivity requires nonempty edge, sample and coverage grids".to_owned());
        }
        let mut settings = vec![self.primary];
        for value in self.edge_counts.iter().copied().collect::<BTreeSet<_>>() {
            settings.push(Settings {
                minimum_edge_count: value,
                ..self.primary
            });
        }
        for value in self.sample_counts.iter().copied().collect::<BTreeSet<_>>() {
            settings.push(Settings {
                minimum_samples_per_edge: value,
                ..self.primary
            });
        }
        for value in self.gene_support.iter().copied().collect::<BTreeSet<_>>() {
            settings.push(Settings {
                minimum_gene_support: value,
                ..self.primary
            });
        }
        let mut seen = BTreeSet::new();
        Ok(settings
            .into_iter()
            .filter(|setting| seen.insert(*setting))
            .enumerate()
            .map(|(index, setting)| {
                (
                    if index == 0 {
                        "primary".to_owned()
                    } else {
                        format!("scenario_{index:03}")
                    },
                    setting,
                )
            })
            .collect())
    }
}

/// Export preparation, comparisons and optional statistics in each scenario.
/// Every requested gene remains in sensitivity.tsv, even if it never qualifies.
pub fn run_sensitivity(
    inputs: &Inputs,
    options: &SensitivityOptions,
    output: &Path,
    command: &[String],
) -> Result<usize, String> {
    let cohort = Cohort::load(inputs)?;
    let scenarios = options.scenarios()?;
    for (_, settings) in &scenarios {
        settings.validate(cohort.samples.len())?;
    }
    if let Some(stats) = &options.statistics {
        if inputs.genes.is_none() || stats.family.trim().is_empty() || stats.min_per_group < 2 {
            return Err(
                "statistics require --genes, a nonempty --family and --min-per-group >= 2"
                    .to_owned(),
            );
        }
        validate_independent_metadata(&inputs.manifest)?;
        if [Condition::Case, Condition::Control]
            .iter()
            .any(|condition| {
                cohort
                    .samples
                    .iter()
                    .filter(|sample| sample.condition == *condition)
                    .count()
                    < stats.min_per_group
            })
        {
            return Err("manifest does not meet --min-per-group".to_owned());
        }
    }
    let stats_json = options.statistics.as_ref().map(|stats| format!("{{\"family\":{},\"alternative\":{},\"design\":\"independent (caller-declared)\",\"min_per_group\":{},\"role\":\"sensitivity; no automatic primary selection or correction across scenarios\"}}", json_string(&stats.family), json_string(alternative_name(stats.alternative)), stats.min_per_group)).unwrap_or_else(|| "null".to_owned());
    let settings_json = format!(
        "{{\"primary\":{},\"scenarios\":[{}],\"statistics\":{}}}",
        options.primary.json(),
        scenarios
            .iter()
            .map(|(_, setting)| setting.json())
            .collect::<Vec<_>>()
            .join(","),
        stats_json
    );
    let mut run = NativeRun::new(
        output,
        "shared graph one-factor sensitivity",
        inputs,
        &settings_json,
        command,
    )?;
    let result = (|| {
        run.write("sensitivity.tsv", "scenario\tminimum_junction_count\tminimum_samples_per_edge\tminimum_gene_support\tgene_id\tstatus\tshared_edges\tcycle_rank\tcovered_samples\trequired_samples\tcase_n\tcontrol_n\tmean_difference_case_minus_control\n")?;
        for (name, settings) in &scenarios {
            let directory = output.join(name);
            run.scenarios.push(format!(
                "{{\"scenario\":{},\"status\":\"running\",\"settings\":{}}}",
                json_string(name),
                settings.json()
            ));
            run.save("running", None)?;
            let prepared = cohort.build(*settings)?;
            export_prepared(inputs, &cohort, &prepared, *settings, &directory, command)?;
            run.files
                .insert(PathBuf::from(name).join("run_manifest.json"));
            let results = crate::analyze_junctions(prepared.junctions.clone())?;
            if results
                .iter()
                .any(|result| !result.decomposition.checks.passed)
            {
                return Err(format!(
                    "{name}: mathematical checks failed; no complete analysis published"
                ));
            }
            let summaries = compare_all_genes(&cohort.samples, &results);
            if !results.is_empty() {
                write_compare_bundle(
                    &directory.join("results"),
                    &directory.join("junctions.tsv"),
                    &directory.join("samples.tsv"),
                    &results,
                    &summaries,
                    command,
                )?;
                run.files
                    .insert(PathBuf::from(name).join("results/run_manifest.json"));
            }
            let by_gene: BTreeMap<_, _> = summaries
                .iter()
                .map(|summary| (summary.comparison.gene_id.as_str(), summary))
                .collect();
            let mut rows = String::new();
            for gene in &prepared.genes {
                let summary = by_gene.get(gene.gene_id.as_str());
                if gene.status == "eligible"
                    && summary.is_none_or(|summary| summary.status != ComparisonStatus::Complete)
                {
                    return Err(format!(
                        "{name}: eligible gene lacks a complete Rust comparison"
                    ));
                }
                let cases = summary
                    .and_then(|summary| summary.case.as_ref())
                    .map_or(0, |group| group.n);
                let controls = summary
                    .and_then(|summary| summary.control.as_ref())
                    .map_or(0, |group| group.n);
                let effect = if gene.status == "eligible" {
                    summary
                        .and_then(|summary| summary.mean_difference)
                        .map(|effect| format!("{effect:.10}"))
                        .unwrap_or_default()
                } else {
                    String::new()
                };
                rows.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    name,
                    settings.minimum_edge_count,
                    settings.minimum_samples_per_edge,
                    settings.minimum_gene_support,
                    gene.gene_id,
                    gene.status,
                    gene.shared_edges,
                    gene.cycle_rank,
                    gene.covered_samples,
                    gene.required_samples,
                    cases,
                    controls,
                    effect
                ));
            }
            OpenOptions::new()
                .append(true)
                .open(output.join("sensitivity.tsv"))
                .and_then(|mut file| file.write_all(rows.as_bytes()))
                .map_err(|e| format!("could not append sensitivity table: {e}"))?;
            let stats_status = if let Some(stats) = &options.statistics {
                if prepared.eligible_genes.is_empty() {
                    "no_eligible_genes"
                } else {
                    let family = directory.join("eligible_genes.tsv");
                    let selected: BTreeSet<_> =
                        prepared.eligible_genes.iter().map(String::as_str).collect();
                    let results: Vec<_> = results
                        .into_iter()
                        .filter(|result| selected.contains(result.gene_id.as_str()))
                        .collect();
                    let mut statistics = test_gene_family(
                        &cohort.samples,
                        &results,
                        &prepared.eligible_genes,
                        StatsOptions {
                            family_id: stats.family.clone(),
                            alternative: stats.alternative,
                            min_per_group: stats.min_per_group,
                            rank_by: None,
                        },
                    )?;
                    statistics.annotations = crate::annotations::read_annotations(&family)?;
                    write_stats_bundle(
                        &directory.join("statistics"),
                        &directory.join("junctions.tsv"),
                        &directory.join("samples.tsv"),
                        &family,
                        &results,
                        &statistics,
                        command,
                    )?;
                    run.files
                        .insert(PathBuf::from(name).join("statistics/run_manifest.json"));
                    "complete"
                }
            } else {
                "not_requested"
            };
            *run.scenarios.last_mut().expect("active scenario") = format!(
                "{{\"scenario\":{},\"status\":\"complete\",\"settings\":{},\"eligible_genes\":{},\"retained_rows\":{},\"statistics_status\":{}}}",
                json_string(name),
                settings.json(),
                prepared.eligible_genes.len(),
                prepared.junctions.len(),
                json_string(stats_status)
            );
            run.save("running", None)?;
        }
        run.summary = format!(
            "{{\"scenarios\":{},\"requested_genes\":{}}}",
            scenarios.len(),
            cohort.family.len()
        );
        Ok(())
    })();
    if result.is_err() {
        if let Some(last) = run.scenarios.last_mut() {
            *last = last.replace("\"status\":\"running\"", "\"status\":\"failed\"");
        }
    }
    run.finish(result)?;
    Ok(scenarios.len())
}
