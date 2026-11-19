//! Strict saved configuration, read-only preflight and a frozen analysis run.
//! Pipeline provenance is checked from small manifests/counts, never raw reads.
use super::*;
use serde::{Deserialize, Serialize};
use serde_json::json;

#[cfg(test)]
mod tests;

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Configuration {
    schema_version: u32,
    input_kind: InputKind,
    junctions: PathBuf,
    manifest: PathBuf,
    genes: Option<PathBuf>,
    output: PathBuf,
    complete_counts: bool,
    shared_graph: Thresholds,
    sensitivity: Option<Grid>,
    statistics: Option<Statistics>,
    provenance: Option<Provenance>,
}

#[derive(Deserialize, Serialize, PartialEq)]
#[serde(rename_all = "snake_case")]
enum InputKind {
    Synthetic,
    Pipeline,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Thresholds {
    minimum_edge_count: u64,
    minimum_samples_per_edge: usize,
    minimum_gene_support: u64,
}

impl Thresholds {
    fn settings(&self) -> Settings {
        Settings {
            minimum_edge_count: self.minimum_edge_count,
            minimum_samples_per_edge: self.minimum_samples_per_edge,
            minimum_gene_support: self.minimum_gene_support,
        }
    }
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Grid {
    edge_counts: Vec<u64>,
    sample_counts: Vec<usize>,
    gene_support: Vec<u64>,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Statistics {
    family: String,
    alternative: String,
    design: String,
    min_per_group: usize,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Provenance {
    batch_manifest: PathBuf,
    reference_manifest: PathBuf,
    gene_map_manifest: PathBuf,
    node_rule: String,
    strandedness: String,
    minimum_anchor: u64,
    exclude_duplicates: bool,
}

#[derive(Deserialize, Serialize, Clone)]
struct Identity {
    path: PathBuf,
    sha256: String,
}

#[derive(Deserialize)]
struct ReferenceManifest {
    status: String,
    method: String,
    inputs: Vec<Identity>,
    outputs: Vec<Identity>,
}
#[derive(Deserialize)]
struct MapManifest {
    schema_version: u32,
    status: String,
    node_rule: String,
    input: Identity,
    output: Identity,
}
#[derive(Deserialize)]
struct BatchManifest {
    pipeline: String,
    status: String,
    samples: Vec<BatchSample>,
    outputs: Vec<Identity>,
}
#[derive(Deserialize)]
struct BatchSample {
    sample_id: String,
    condition: String,
    status: String,
    run_manifest: Identity,
}
#[derive(Deserialize)]
struct SampleManifest {
    pipeline: String,
    status: String,
    settings: Methods,
    inputs: Vec<Identity>,
    outputs: Vec<Identity>,
    commands: Vec<Vec<String>>,
    stages: Vec<Stage>,
}
#[derive(Deserialize)]
struct Stage {
    status: String,
}
#[derive(Deserialize)]
struct Methods {
    sample_id: String,
    strandedness: String,
    minimum_unique_reads: u64,
    minimum_anchor: u64,
    exclude_duplicates: bool,
    count_unit: String,
}

struct Plan {
    base: PathBuf,
    configuration: Configuration,
    original: String,
    inputs: Inputs,
    output: PathBuf,
    options: SensitivityOptions,
    checked: BTreeMap<PathBuf, String>,
    report: serde_json::Value,
}

fn fail(message: impl Into<String>) -> Result<(), String> {
    Err(message.into())
}

fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, String> {
    let text =
        fs::read_to_string(path).map_err(|e| format!("could not read {}: {e}", path.display()))?;
    serde_json::from_str(&text).map_err(|e| format!("invalid JSON in {}: {e}", path.display()))
}

fn require(condition: bool, message: impl Into<String>) -> Result<(), String> {
    if condition { Ok(()) } else { fail(message) }
}

fn digest_matches(actual: &str, expected: &str) -> bool {
    expected.len() == 64
        && expected.bytes().all(|byte| byte.is_ascii_hexdigit())
        && actual == expected
}

fn resolve(base: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_owned()
    } else {
        base.join(path)
    }
}

fn checked_file(path: &Path, checked: &mut BTreeMap<PathBuf, String>) -> Result<PathBuf, String> {
    let canonical = path
        .canonicalize()
        .map_err(|e| format!("cannot access {}: {e}", path.display()))?;
    require(
        canonical.is_file(),
        format!("not a regular file: {}", canonical.display()),
    )?;
    checked.insert(canonical.clone(), sha256_file(&canonical)?);
    Ok(canonical)
}

impl Configuration {
    fn load(path: &Path) -> Result<Plan, String> {
        let path = path
            .canonicalize()
            .map_err(|e| format!("cannot access config {}: {e}", path.display()))?;
        let base = path.parent().expect("canonical config parent");
        let original = fs::read_to_string(&path).map_err(|e| e.to_string())?;
        let configuration: Self = serde_json::from_str(&original)
            .map_err(|e| format!("invalid configuration {}: {e}", path.display()))?;
        require(
            configuration.schema_version == 1,
            "only configuration schema_version 1 is supported",
        )?;
        require(
            !configuration.output.as_os_str().is_empty(),
            "configuration output cannot be empty",
        )?;
        let output = resolve(base, &configuration.output);
        require(
            output.symlink_metadata().is_err(),
            format!(
                "output path {} already exists; choose a new folder",
                output.display()
            ),
        )?;
        // Check that the nearest existing parent is a directory, without creating it.
        let mut ancestor = output.parent().ok_or("output needs a parent directory")?;
        while !ancestor.exists() {
            ancestor = ancestor.parent().ok_or("output has no existing parent")?;
        }
        require(ancestor.is_dir(), "output parent is not a directory")?;
        let mut checked = BTreeMap::new();
        checked_file(&path, &mut checked)?;
        let inputs = Inputs {
            junctions: checked_file(&resolve(base, &configuration.junctions), &mut checked)?,
            manifest: checked_file(&resolve(base, &configuration.manifest), &mut checked)?,
            genes: configuration
                .genes
                .as_ref()
                .map(|file| checked_file(&resolve(base, file), &mut checked))
                .transpose()?,
            complete_counts: configuration.complete_counts,
        };
        let cohort = Cohort::load(&inputs)?;
        let primary = configuration.shared_graph.settings();
        let statistics = configuration.statistics.as_ref().map(|stats| {
            require(stats.design == "independent", "configuration supports only independent statistics")?;
            require(inputs.genes.is_some() && !stats.family.trim().is_empty() && stats.min_per_group >= 2,
                    "statistics require a gene family, nonempty family label and min_per_group >= 2")?;
            crate::statistics::validate_independent_metadata(&inputs.manifest)?;
            require([Condition::Case, Condition::Control].iter().all(|condition| cohort.samples.iter().filter(|sample| sample.condition == *condition).count() >= stats.min_per_group), "cohort does not meet min_per_group")?;
            Ok::<_, String>(StatisticsSettings { family: stats.family.clone(), alternative: crate::statistics::parse_alternative(&stats.alternative)?, min_per_group: stats.min_per_group })
        }).transpose()?;
        let options = SensitivityOptions {
            primary,
            edge_counts: configuration.sensitivity.as_ref().map_or_else(
                || vec![primary.minimum_edge_count],
                |grid| grid.edge_counts.clone(),
            ),
            sample_counts: configuration.sensitivity.as_ref().map_or_else(
                || vec![primary.minimum_samples_per_edge],
                |grid| grid.sample_counts.clone(),
            ),
            gene_support: configuration.sensitivity.as_ref().map_or_else(
                || vec![primary.minimum_gene_support],
                |grid| grid.gene_support.clone(),
            ),
            statistics,
        };
        let scenarios = options.scenarios()?;
        for (_, settings) in &scenarios {
            settings.validate(cohort.samples.len())?;
        }
        match configuration.input_kind {
            InputKind::Synthetic => require(
                configuration.provenance.is_none(),
                "synthetic configurations must have provenance null; use pipeline mode for biological exports",
            )?,
            InputKind::Pipeline => verify_pipeline(
                configuration
                    .provenance
                    .as_ref()
                    .ok_or("pipeline input_kind requires provenance")?,
                base,
                &inputs,
                &cohort,
                &mut checked,
            )?,
        }
        let prepared = cohort.build(primary)?;
        let report = json!({
            "status": "PASS", "schema_version": 1, "configuration": path,
            "input_kind": configuration.input_kind,
            "input_identity": if configuration.input_kind == InputKind::Synthetic { "synthetic caller declaration; no biological provenance verified" } else { "completed batch/sample counts and matching recorded reference/map/method identities verified" },
            "cohort_samples": cohort.samples.len(), "requested_genes": cohort.family.len(),
            "primary_eligible_genes": prepared.eligible_genes.len(), "scenarios": scenarios.len(),
            "output": output, "checked_files": checked.iter().map(|(path, sha)| json!({"path":path,"sha256":sha})).collect::<Vec<_>>(),
            "limitations": ["Checks do not approve biological thresholds, donor independence or study conclusions.", "Reference sequence/index hashes are compared from recorded manifests; large raw/reference files are not reread.", "Human usability review is separate from this software check."]
        });
        Ok(Plan {
            base: base.to_path_buf(),
            configuration,
            original,
            inputs,
            output,
            options,
            checked,
            report,
        })
    }
}

fn unique_output<'a>(files: &'a [Identity], name: &str) -> Result<&'a Identity, String> {
    let matches: Vec<_> = files
        .iter()
        .filter(|record| record.path.file_name().is_some_and(|file| file == name))
        .collect();
    require(
        matches.len() == 1,
        format!("provenance needs exactly one {name} output"),
    )?;
    Ok(matches[0])
}

fn command_argument<'a>(commands: &'a [Vec<String>], flag: &str) -> Result<&'a str, String> {
    let values: Vec<_> = commands
        .iter()
        .flat_map(|command| command.windows(2))
        .filter(|pair| pair[0] == flag)
        .map(|pair| pair[1].as_str())
        .collect();
    require(
        values.len() == 1,
        format!("provenance needs exactly one recorded {flag} argument"),
    )?;
    Ok(values[0])
}

fn command_value(commands: &[Vec<String>], flag: &str) -> Result<PathBuf, String> {
    command_argument(commands, flag).map(PathBuf::from)
}

fn input_hash<'a>(files: &'a [Identity], path: &Path) -> Result<&'a str, String> {
    let found: Vec<_> = files.iter().filter(|record| record.path == path).collect();
    require(
        found.len() == 1,
        format!(
            "missing or ambiguous recorded input identity: {}",
            path.display()
        ),
    )?;
    Ok(&found[0].sha256)
}

fn index_hashes(files: &[Identity]) -> Result<Vec<String>, String> {
    let mut hashes: Vec<_> = files
        .iter()
        .filter(|record| {
            matches!(
                record.path.extension().and_then(|x| x.to_str()),
                Some("ht2" | "ht2l")
            )
        })
        .map(|record| record.sha256.clone())
        .collect();
    require(
        hashes.len() == 8 && hashes.iter().all(|value| digest_matches(value, value)),
        "reference provenance requires eight valid index-part hashes",
    )?;
    hashes.sort();
    Ok(hashes)
}

fn verify_pipeline(
    provenance: &Provenance,
    base: &Path,
    inputs: &Inputs,
    cohort: &Cohort,
    checked: &mut BTreeMap<PathBuf, String>,
) -> Result<(), String> {
    require(
        matches!(
            provenance.node_rule.as_str(),
            "shared-splice-site" | "exon-interval"
        ),
        "unknown expected node_rule",
    )?;
    require(
        matches!(
            provenance.strandedness.as_str(),
            "RF" | "FR" | "F" | "R" | "unstranded"
        ) && provenance.minimum_anchor > 0,
        "invalid expected strandedness/minimum_anchor",
    )?;
    let reference_path = checked_file(&resolve(base, &provenance.reference_manifest), checked)?;
    let reference: ReferenceManifest = read_json(&reference_path)?;
    require(
        reference.status == "complete"
            && reference.method == "HISAT2 reference preparation"
            && reference.inputs.len() == 2,
        "reference manifest is incomplete or unsupported",
    )?;
    let expected_index = index_hashes(&reference.outputs)?;
    let map_path = checked_file(&resolve(base, &provenance.gene_map_manifest), checked)?;
    let map: MapManifest = read_json(&map_path)?;
    require(
        map.schema_version == 1
            && map.status == "complete"
            && map.node_rule == provenance.node_rule,
        "map provenance has wrong schema, completion status or node rule",
    )?;
    require(
        digest_matches(&map.input.sha256, &reference.inputs[1].sha256),
        "map GTF differs from the configured reference GTF",
    )?;
    let expected_hints = &unique_output(&reference.outputs, "splice_sites.tsv")?.sha256;
    let batch_path = checked_file(&resolve(base, &provenance.batch_manifest), checked)?;
    let batch: BatchManifest = read_json(&batch_path)?;
    require(
        batch.status == "complete"
            && batch.pipeline == "Splice Girl sequential Cleaver/HISAT2 batch",
        "batch is incomplete or has an unsupported pipeline",
    )?;
    let combined = unique_output(&batch.outputs, "junctions.tsv")?;
    require(
        digest_matches(&sha256_file(&inputs.junctions)?, &combined.sha256),
        "combined junction table does not match the completed batch output",
    )?;
    require(
        batch.samples.len() == cohort.samples.len(),
        "batch sample count does not match analysis manifest",
    )?;
    let mut seen = BTreeSet::new();
    for sample in &batch.samples {
        require(
            !sample.sample_id.is_empty()
                && sample.sample_id.as_bytes()[0].is_ascii_alphanumeric()
                && sample
                    .sample_id
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || b"._-".contains(&byte)),
            "invalid pipeline sample ID",
        )?;
        let expected = cohort
            .samples
            .iter()
            .find(|row| row.sample_id == sample.sample_id)
            .ok_or_else(|| {
                format!(
                    "batch sample {} is absent from analysis manifest",
                    sample.sample_id
                )
            })?;
        require(
            seen.insert(&sample.sample_id)
                && sample.status == "complete"
                && sample.condition == expected.condition.to_string(),
            format!(
                "sample {} is duplicate, incomplete or has a conflicting condition",
                sample.sample_id
            ),
        )?;
        let local = batch_path
            .parent()
            .unwrap()
            .join(&sample.sample_id)
            .join("run_manifest.json");
        let path = if local.exists() {
            local
        } else {
            resolve(batch_path.parent().unwrap(), &sample.run_manifest.path)
        };
        let path = checked_file(&path, checked)?;
        require(
            digest_matches(&sha256_file(&path)?, &sample.run_manifest.sha256),
            format!("sample {} run manifest checksum mismatch", sample.sample_id),
        )?;
        let run: SampleManifest = read_json(&path)?;
        require(
            run.status == "complete"
                && run.pipeline == "Splice Girl Cleaver/HISAT2"
                && !run.stages.is_empty()
                && run.stages.iter().all(|stage| stage.status == "complete"),
            format!("sample {} processing is incomplete", sample.sample_id),
        )?;
        let methods = run.settings;
        require(
            methods.sample_id == sample.sample_id
                && methods.strandedness == provenance.strandedness
                && methods.minimum_anchor == provenance.minimum_anchor
                && methods.minimum_unique_reads == 1
                && methods.exclude_duplicates == provenance.exclude_duplicates
                && methods.count_unit
                    == "primary uniquely mapped read alignments with NH:i:1 and XS strand",
            format!(
                "sample {} preprocessing settings differ from configuration (support must be 1)",
                sample.sample_id
            ),
        )?;
        let gene_map = command_value(&run.commands, "--gene-map")?;
        let hints = command_value(&run.commands, "--known-splicesite-infile")?;
        require(
            command_argument(&run.commands, "--min-anchor")? == methods.minimum_anchor.to_string()
                && command_argument(&run.commands, "--min-unique")? == "1",
            "recorded counter command contradicts anchor/support settings",
        )?;
        let duplicate_flags = run
            .commands
            .iter()
            .flatten()
            .filter(|argument| *argument == "--exclude-duplicates")
            .count();
        require(
            duplicate_flags == usize::from(methods.exclude_duplicates),
            "recorded counter command contradicts duplicate policy",
        )?;
        if methods.strandedness == "unstranded" {
            require(
                !run.commands
                    .iter()
                    .flatten()
                    .any(|argument| argument == "--rna-strandness"),
                "recorded alignment command contradicts unstranded settings",
            )?;
        } else {
            require(
                command_argument(&run.commands, "--rna-strandness")? == methods.strandedness,
                "recorded alignment command contradicts strandedness settings",
            )?;
        }
        let index_prefix = command_value(&run.commands, "-x")?;
        let ordinary = PathBuf::from(format!("{}.1.ht2", index_prefix.display()));
        let suffix = if run.inputs.iter().any(|record| record.path == ordinary) {
            "ht2"
        } else {
            "ht2l"
        };
        for part in 1..=8 {
            input_hash(
                &run.inputs,
                &PathBuf::from(format!("{}.{part}.{suffix}", index_prefix.display())),
            )?;
        }
        require(
            digest_matches(input_hash(&run.inputs, &gene_map)?, &map.output.sha256),
            format!("sample {} uses a different gene map", sample.sample_id),
        )?;
        require(
            digest_matches(input_hash(&run.inputs, &hints)?, expected_hints)
                && index_hashes(&run.inputs)? == expected_index,
            format!(
                "sample {} reference index or splice hints differ",
                sample.sample_id
            ),
        )?;
        let count_record = unique_output(&run.outputs, "junctions.tsv")?;
        let local = path.parent().unwrap().join("junctions.tsv");
        let counts = checked_file(&local, checked)?;
        require(
            digest_matches(&sha256_file(&counts)?, &count_record.sha256),
            format!(
                "sample {} junction output checksum mismatch",
                sample.sample_id
            ),
        )?;
        // Prove that the merged input contains exactly this completed sample's
        // rows, rather than trusting a completion flag or aggregate hash alone.
        let raw = Table::read(
            &counts,
            &["sample_id", "gene_id", "from", "to", "count"],
            true,
        )?;
        let mut actual = BTreeMap::new();
        for row in raw.rows {
            let count = row[4]
                .parse::<u64>()
                .map_err(|_| "sample export contains a noninteger count")?;
            require(
                row[0] == sample.sample_id
                    && actual
                        .insert((row[1].clone(), row[2].clone(), row[3].clone()), count)
                        .is_none(),
                "sample export has mismatched identity or duplicate edges",
            )?;
        }
        let expected: BTreeMap<_, _> = cohort
            .counts
            .iter()
            .flat_map(|(gene, edges)| {
                edges.iter().filter_map(|(edge, samples)| {
                    samples
                        .get(&sample.sample_id)
                        .map(|count| ((gene.clone(), edge.0.clone(), edge.1.clone()), *count))
                })
            })
            .collect();
        require(
            actual == expected,
            format!(
                "sample {} counts differ between the combined and per-sample exports",
                sample.sample_id
            ),
        )?;
    }
    Ok(())
}

impl Plan {
    fn verify_frozen(&self) -> Result<(), String> {
        for (path, digest) in &self.checked {
            require(
                sha256_file(path)? == *digest,
                format!("checked input changed during analysis: {}", path.display()),
            )?;
        }
        Ok(())
    }

    fn execute(self, command: &[String]) -> Result<(), String> {
        self.verify_frozen()?;
        let settings = serde_json::to_string(
            &json!({"configuration":self.configuration,"preflight":self.report}),
        )
        .map_err(|e| e.to_string())?;
        let mut run = NativeRun::new(
            &self.output,
            "configured shared graph analysis",
            &self.inputs,
            &settings,
            command,
        )?;
        let result = (|| {
            run.write("configuration.json", &self.original)?;
            run.write(
                "preflight.json",
                &serde_json::to_string_pretty(&self.report).map_err(|e| e.to_string())?,
            )?;
            let mut resolved =
                serde_json::to_value(&self.configuration).map_err(|e| e.to_string())?;
            resolved["junctions"] = json!(self.inputs.junctions);
            resolved["manifest"] = json!(self.inputs.manifest);
            resolved["genes"] = json!(self.inputs.genes);
            resolved["output"] = json!(self.output);
            if let Some(provenance) = &self.configuration.provenance {
                for (name, path) in [
                    ("batch_manifest", &provenance.batch_manifest),
                    ("reference_manifest", &provenance.reference_manifest),
                    ("gene_map_manifest", &provenance.gene_map_manifest),
                ] {
                    resolved["provenance"][name] = json!(resolve(&self.base, path));
                }
            }
            run.write(
                "resolved_configuration.json",
                &serde_json::to_string_pretty(&resolved).map_err(|e| e.to_string())?,
            )?;
            let analysis = self.output.join("analysis");
            run_sensitivity(&self.inputs, &self.options, &analysis, command)?;
            run.files
                .insert(PathBuf::from("analysis/run_manifest.json"));
            let primary = analysis.join("primary");
            let rows_path = primary.join("junctions.tsv");
            let rows = Table::read(
                &rows_path,
                &["sample_id", "gene_id", "from", "to", "count"],
                true,
            )?;
            if !rows.rows.is_empty() {
                let results = crate::analyze_junctions(crate::parser::parse_tsv(&rows_path)?)?;
                let samples = crate::manifest::parse_manifest(&self.inputs.manifest)?;
                let summaries = crate::compare::compare_all_genes(&samples, &results);
                crate::output::write_report_bundle(
                    &self.output.join("primary_report"),
                    &rows_path,
                    &self.inputs.manifest,
                    &results,
                    &summaries,
                    command,
                )?;
                run.files
                    .insert(PathBuf::from("primary_report/run_manifest.json"));
            }
            let synthetic = self.configuration.input_kind == InputKind::Synthetic;
            let label = if synthetic {
                "Synthetic demonstration — not biological findings"
            } else {
                "Pipeline identities verified; biological interpretation requires review"
            };
            let mut html = format!(
                "<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>Splice Girl configured analysis</title><body><h1>Splice Girl</h1><p>{label}</p><p>All requested genes and exclusions remain in the eligibility and sensitivity tables.</p><ul><li><a href=\"configuration.json\">Original configuration</a></li><li><a href=\"preflight.json\">Preflight and checked identities</a></li><li><a href=\"analysis/sensitivity.tsv\">Sensitivity results</a></li><li><a href=\"analysis/primary/gene_eligibility.tsv\">Primary gene eligibility</a></li><li><a href=\"analysis/primary/sample_eligibility.tsv\">Primary sample eligibility</a></li>"
            );
            if self.output.join("primary_report/report.html").is_file() {
                html.push_str("<li><a href=\"primary_report/report.html\">Primary descriptive report</a></li>");
            }
            if primary.join("statistics/report.html").is_file() {
                html.push_str("<li><a href=\"analysis/primary/statistics/report.html\">Primary statistics report</a></li>");
            }
            html.push_str(&format!(
                "</ul><p>A cycle component does not prove circular RNA or observed use of a zero-count edge. Software is licensed under {}. Please cite the version or commit used.</p></body></html>\n",
                env!("CARGO_PKG_LICENSE")
            ));
            run.write("report.html", &html)?;
            self.verify_frozen()?;
            run.summary = format!(
                "{{\"primary_eligible_genes\":{},\"scenarios\":{}}}",
                self.report["primary_eligible_genes"], self.report["scenarios"]
            );
            Ok(())
        })();
        run.finish(result)?;
        println!(
            "Configured analysis complete: {}",
            self.output.join("report.html").display()
        );
        Ok(())
    }
}

pub fn run_config_command(name: &str, args: &[String], command: &[String]) -> Result<(), String> {
    let help = "usage: splice-girl check --config <analysis.json>\n       splice-girl run --config <analysis.json>\n\ncheck validates configuration, inputs, eligibility and pipeline identities without writing files.\nrun repeats preflight, saves configuration/checksums, and generates report.html in a new output folder.\nPaths are relative to the configuration file, regardless of your working folder.\nUse splice-girl example --output demo-inputs for a labeled synthetic configuration.";
    if args == ["--help"] || args == ["-h"] {
        println!("{help}");
        return Ok(());
    }
    require(args.len() == 2 && args[0] == "--config", help)?;
    let plan = Configuration::load(Path::new(&args[1]))?;
    if name == "check" {
        println!(
            "Preflight: PASS\nInput kind: {}\nSamples: {}\nPrimary eligible genes: {}\nScenarios: {}\nNo output files written.",
            plan.report["input_kind"],
            plan.report["cohort_samples"],
            plan.report["primary_eligible_genes"],
            plan.report["scenarios"]
        );
        Ok(())
    } else {
        plan.execute(command)
    }
}
