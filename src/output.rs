//! Machine-readable result exports and reproducibility manifests (R13).
//!
//! The analysis code deliberately stays independent of file formats. This
//! module turns completed in-memory results into a small output bundle without
//! changing the graph construction or Hodge calculations.

use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::GroupResult;
use crate::compare::{GeneComparison, GeneComparisonSummary};
use crate::hodge::TOLERANCE;

const SAMPLE_METRICS: &str = "sample_metrics.tsv";
const JUNCTION_COMPONENTS: &str = "junction_components.tsv";
const COMPARISONS: &str = "comparisons.tsv";
const GENE_SUMMARIES: &str = "gene_summaries.tsv";
const RUN_MANIFEST: &str = "run_manifest.json";

/// Identifies the files written by one export operation.
#[derive(Debug)]
pub struct ExportSummary {
    pub run_id: String,
    pub output_dir: PathBuf,
    pub files: Vec<PathBuf>,
}

#[derive(Debug)]
struct InputRecord {
    role: &'static str,
    path: String,
    sha256: String,
}

/// Writes the outputs for `splice-girl analyze`.
pub fn write_analysis_bundle(
    output_dir: &Path,
    junctions_path: &Path,
    results: &[GroupResult],
    command: &[String],
) -> Result<ExportSummary, String> {
    let metadata = RunMetadata::new("analyze", command, &[("junctions", junctions_path)])?;
    let result_refs: Vec<&GroupResult> = results.iter().collect();
    let metrics = sample_metrics_tsv(&metadata.run_id, &result_refs);
    let components = junction_components_tsv(&metadata.run_id, &result_refs);
    let failures = failed_checks(results);
    let manifest = run_manifest_json(
        &metadata,
        &[SAMPLE_METRICS, JUNCTION_COMPONENTS, RUN_MANIFEST],
        results.len(),
        results.iter().map(|result| result.graph.edges.len()).sum(),
        0,
        0,
        &[],
        &failures,
    );

    write_bundle(
        output_dir,
        &metadata.run_id,
        &[
            (SAMPLE_METRICS, metrics),
            (JUNCTION_COMPONENTS, components),
            (RUN_MANIFEST, manifest),
        ],
    )
}

/// Writes the outputs for `splice-girl gene`. Metrics and junction components
/// are limited to the requested gene, while `comparisons.tsv` also records
/// every explicit exclusion.
pub fn write_gene_bundle(
    output_dir: &Path,
    junctions_path: &Path,
    sample_manifest_path: &Path,
    results: &[GroupResult],
    comparison: &GeneComparison,
    command: &[String],
) -> Result<ExportSummary, String> {
    let metadata = RunMetadata::new(
        "gene",
        command,
        &[
            ("junctions", junctions_path),
            ("sample_manifest", sample_manifest_path),
        ],
    )?;
    let selected: Vec<&GroupResult> = results
        .iter()
        .filter(|result| result.gene_id == comparison.gene_id)
        .collect();
    let metrics = sample_metrics_tsv(&metadata.run_id, &selected);
    let components = junction_components_tsv(&metadata.run_id, &selected);
    let comparisons = comparisons_tsv(&metadata.run_id, comparison);
    let failures: Vec<FailureRecord> = selected
        .iter()
        .filter(|result| !result.decomposition.checks.passed)
        .map(|result| FailureRecord {
            sample_id: result.sample_id.clone(),
            gene_id: result.gene_id.clone(),
            reason: "mathematical checks failed".to_owned(),
        })
        .collect();
    let exclusions: Vec<ExclusionRecord> = comparison
        .excluded
        .iter()
        .map(|sample| ExclusionRecord {
            gene_id: comparison.gene_id.clone(),
            sample_id: sample.sample_id.clone(),
            reason: sample.reason.clone(),
        })
        .collect();
    let manifest = run_manifest_json(
        &metadata,
        &[
            SAMPLE_METRICS,
            JUNCTION_COMPONENTS,
            COMPARISONS,
            RUN_MANIFEST,
        ],
        selected.len(),
        selected.iter().map(|result| result.graph.edges.len()).sum(),
        comparison.included.len(),
        comparison.excluded.len(),
        &exclusions,
        &failures,
    );

    write_bundle(
        output_dir,
        &metadata.run_id,
        &[
            (SAMPLE_METRICS, metrics),
            (JUNCTION_COMPONENTS, components),
            (COMPARISONS, comparisons),
            (RUN_MANIFEST, manifest),
        ],
    )
}

/// Writes a sample-preserving, descriptive comparison for every gene. The
/// bundle deliberately contains no p-values or ranking because those require
/// a separately approved statistical specification.
pub fn write_compare_bundle(
    output_dir: &Path,
    junctions_path: &Path,
    sample_manifest_path: &Path,
    results: &[GroupResult],
    summaries: &[GeneComparisonSummary],
    command: &[String],
) -> Result<ExportSummary, String> {
    let metadata = RunMetadata::new(
        "compare",
        command,
        &[
            ("junctions", junctions_path),
            ("sample_manifest", sample_manifest_path),
        ],
    )?;
    let result_refs: Vec<&GroupResult> = results.iter().collect();
    let metrics = sample_metrics_tsv(&metadata.run_id, &result_refs);
    let components = junction_components_tsv(&metadata.run_id, &result_refs);
    let comparisons = all_comparisons_tsv(&metadata.run_id, summaries);
    let gene_summaries = gene_summaries_tsv(&metadata.run_id, summaries);
    let failures = failed_checks(results);
    let exclusions: Vec<ExclusionRecord> = summaries
        .iter()
        .flat_map(|summary| {
            summary
                .comparison
                .excluded
                .iter()
                .map(|sample| ExclusionRecord {
                    gene_id: summary.comparison.gene_id.clone(),
                    sample_id: sample.sample_id.clone(),
                    reason: sample.reason.clone(),
                })
        })
        .collect();
    let included_samples: usize = summaries
        .iter()
        .map(|summary| summary.comparison.included.len())
        .sum();
    let excluded_samples: usize = summaries
        .iter()
        .map(|summary| summary.comparison.excluded.len())
        .sum();
    let manifest = run_manifest_json(
        &metadata,
        &[
            SAMPLE_METRICS,
            JUNCTION_COMPONENTS,
            COMPARISONS,
            GENE_SUMMARIES,
            RUN_MANIFEST,
        ],
        results.len(),
        results.iter().map(|result| result.graph.edges.len()).sum(),
        included_samples,
        excluded_samples,
        &exclusions,
        &failures,
    );

    write_bundle(
        output_dir,
        &metadata.run_id,
        &[
            (SAMPLE_METRICS, metrics),
            (JUNCTION_COMPONENTS, components),
            (COMPARISONS, comparisons),
            (GENE_SUMMARIES, gene_summaries),
            (RUN_MANIFEST, manifest),
        ],
    )
}

struct RunMetadata {
    run_id: String,
    unix_time_seconds: u64,
    mode: &'static str,
    command: Vec<String>,
    inputs: Vec<InputRecord>,
}

impl RunMetadata {
    fn new(
        mode: &'static str,
        command: &[String],
        inputs: &[(&'static str, &Path)],
    ) -> Result<Self, String> {
        let duration = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(|error| format!("system clock is before the Unix epoch: {error}"))?;
        let unix_time_seconds = duration.as_secs();
        let run_id = format!(
            "splice-girl-{unix_time_seconds}-{:09}",
            duration.subsec_nanos()
        );
        let mut input_records = Vec::with_capacity(inputs.len());
        for (role, path) in inputs {
            input_records.push(InputRecord {
                role,
                path: path.display().to_string(),
                sha256: sha256_file(path)?,
            });
        }
        Ok(Self {
            run_id,
            unix_time_seconds,
            mode,
            command: command.to_vec(),
            inputs: input_records,
        })
    }
}

fn write_bundle(
    output_dir: &Path,
    run_id: &str,
    files: &[(&str, String)],
) -> Result<ExportSummary, String> {
    if output_dir.exists() {
        return Err(format!(
            "output path {} already exists; choose a new folder so existing results are not overwritten",
            output_dir.display()
        ));
    }
    fs::create_dir_all(output_dir).map_err(|error| {
        format!(
            "could not create output folder {}: {error}",
            output_dir.display()
        )
    })?;

    let mut written = Vec::with_capacity(files.len());
    for (name, contents) in files {
        let path = output_dir.join(name);
        fs::write(&path, contents)
            .map_err(|error| format!("could not write {}: {error}", path.display()))?;
        written.push(path);
    }

    Ok(ExportSummary {
        run_id: run_id.to_owned(),
        output_dir: output_dir.to_path_buf(),
        files: written,
    })
}

fn sample_metrics_tsv(run_id: &str, results: &[&GroupResult]) -> String {
    let mut output = String::from(
        "run_id\tsample_id\tgene_id\tstatus\tnode_count\tedge_count\ttotal_support\tsignal_energy\tgradient_energy\tcycle_energy\tcycle_fraction\tchecks_status\treconstruction_error\tkernel_error\torthogonality_error\n",
    );
    for result in results {
        let decomposition = &result.decomposition;
        let total_support: f64 = result.graph.edges.iter().map(|edge| edge.count).sum();
        let checks_status = if decomposition.checks.passed {
            "PASS"
        } else {
            "FAIL"
        };
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.10}\t{:.10}\t{:.10}\t{:.10}\t{:.10}\t{}\t{:.10e}\t{:.10e}\t{:.10e}\n",
            tsv(run_id),
            tsv(&result.sample_id),
            tsv(&result.gene_id),
            if decomposition.checks.passed { "passed" } else { "failed" },
            result.graph.nodes.len(),
            result.graph.edges.len(),
            total_support,
            decomposition.signal_energy,
            decomposition.gradient_energy,
            decomposition.cycle_energy,
            decomposition.cycle_fraction,
            checks_status,
            decomposition.checks.reconstruction_error,
            decomposition.checks.kernel_error,
            decomposition.checks.orthogonality_error,
        ));
    }
    output
}

fn junction_components_tsv(run_id: &str, results: &[&GroupResult]) -> String {
    let mut output = String::from(
        "run_id\tsample_id\tgene_id\tfrom\tto\traw_count\tnormalized_signal\tgradient_component\tcycle_component\n",
    );
    for result in results {
        for (index, edge) in result.graph.edges.iter().enumerate() {
            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{:.10}\t{:.10}\t{:.10}\t{:.10}\n",
                tsv(run_id),
                tsv(&result.sample_id),
                tsv(&result.gene_id),
                tsv(&result.graph.nodes[edge.from]),
                tsv(&result.graph.nodes[edge.to]),
                edge.count,
                result.graph.signal[index],
                result.decomposition.gradient[index],
                result.decomposition.cycle[index],
            ));
        }
    }
    output
}

fn comparisons_tsv(run_id: &str, comparison: &GeneComparison) -> String {
    let mut output = String::from(
        "run_id\tgene_id\tsample_id\tcondition\tstatus\tcycle_fraction\texclusion_reason\n",
    );
    append_comparison_rows(&mut output, run_id, comparison);
    output
}

fn all_comparisons_tsv(run_id: &str, summaries: &[GeneComparisonSummary]) -> String {
    let mut output = String::from(
        "run_id\tgene_id\tsample_id\tcondition\tstatus\tcycle_fraction\texclusion_reason\n",
    );
    for summary in summaries {
        append_comparison_rows(&mut output, run_id, &summary.comparison);
    }
    output
}

fn append_comparison_rows(output: &mut String, run_id: &str, comparison: &GeneComparison) {
    for sample in &comparison.included {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\tincluded\t{:.10}\t\n",
            tsv(run_id),
            tsv(&comparison.gene_id),
            tsv(&sample.sample_id),
            sample.condition,
            sample.cycle_fraction,
        ));
    }
    for sample in &comparison.excluded {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\texcluded\t\t{}\n",
            tsv(run_id),
            tsv(&comparison.gene_id),
            tsv(&sample.sample_id),
            sample
                .condition
                .map(|condition| condition.to_string())
                .unwrap_or_default(),
            tsv(&sample.reason),
        ));
    }
}

fn gene_summaries_tsv(run_id: &str, summaries: &[GeneComparisonSummary]) -> String {
    let mut output = String::from(
        "run_id\tgene_id\tstatus\tcase_n\tcontrol_n\tcase_mean\tcontrol_mean\tmean_difference_case_minus_control\tcase_median\tcontrol_median\tcase_min\tcase_max\tcontrol_min\tcontrol_max\texcluded_count\tnote\n",
    );
    for summary in summaries {
        let note = match summary.status {
            crate::compare::ComparisonStatus::Complete => {
                "descriptive effect only; no p-value or q-value calculated"
            }
            crate::compare::ComparisonStatus::NoEligibleSamples => {
                "no sample was eligible for this gene"
            }
            crate::compare::ComparisonStatus::MissingCase => {
                "no eligible case sample; effect not calculated"
            }
            crate::compare::ComparisonStatus::MissingControl => {
                "no eligible control sample; effect not calculated"
            }
        };
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            tsv(run_id),
            tsv(&summary.comparison.gene_id),
            summary.status,
            summary.case.as_ref().map_or(0, |group| group.n),
            summary.control.as_ref().map_or(0, |group| group.n),
            optional_number(summary.case.as_ref().map(|group| group.mean)),
            optional_number(summary.control.as_ref().map(|group| group.mean)),
            optional_number(summary.mean_difference),
            optional_number(summary.case.as_ref().map(|group| group.median)),
            optional_number(summary.control.as_ref().map(|group| group.median)),
            optional_number(summary.case.as_ref().map(|group| group.minimum)),
            optional_number(summary.case.as_ref().map(|group| group.maximum)),
            optional_number(summary.control.as_ref().map(|group| group.minimum)),
            optional_number(summary.control.as_ref().map(|group| group.maximum)),
            summary.comparison.excluded.len(),
            note,
        ));
    }
    output
}

fn optional_number(value: Option<f64>) -> String {
    value
        .map(|number| format!("{number:.10}"))
        .unwrap_or_default()
}

#[derive(Debug)]
struct ExclusionRecord {
    gene_id: String,
    sample_id: String,
    reason: String,
}

#[derive(Debug)]
struct FailureRecord {
    sample_id: String,
    gene_id: String,
    reason: String,
}

fn failed_checks(results: &[GroupResult]) -> Vec<FailureRecord> {
    results
        .iter()
        .filter(|result| !result.decomposition.checks.passed)
        .map(|result| FailureRecord {
            sample_id: result.sample_id.clone(),
            gene_id: result.gene_id.clone(),
            reason: "mathematical checks failed".to_owned(),
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn run_manifest_json(
    metadata: &RunMetadata,
    output_files: &[&str],
    sample_gene_groups: usize,
    junction_rows: usize,
    included_samples: usize,
    excluded_samples: usize,
    exclusions: &[ExclusionRecord],
    failures: &[FailureRecord],
) -> String {
    let command = metadata
        .command
        .iter()
        .map(|argument| json_string(argument))
        .collect::<Vec<_>>()
        .join(", ");
    let inputs = metadata
        .inputs
        .iter()
        .map(|input| {
            format!(
                "    {{\"role\": {}, \"path\": {}, \"sha256\": {}}}",
                json_string(input.role),
                json_string(&input.path),
                json_string(&input.sha256)
            )
        })
        .collect::<Vec<_>>()
        .join(",\n");
    let outputs = output_files
        .iter()
        .map(|name| json_string(name))
        .collect::<Vec<_>>()
        .join(", ");
    let exclusions_json = exclusions
        .iter()
        .map(|record| {
            format!(
                "    {{\"gene_id\": {}, \"sample_id\": {}, \"reason\": {}}}",
                json_string(&record.gene_id),
                json_string(&record.sample_id),
                json_string(&record.reason)
            )
        })
        .collect::<Vec<_>>()
        .join(",\n");
    let failures_json = failures
        .iter()
        .map(|record| {
            format!(
                "    {{\"sample_id\": {}, \"gene_id\": {}, \"reason\": {}}}",
                json_string(&record.sample_id),
                json_string(&record.gene_id),
                json_string(&record.reason)
            )
        })
        .collect::<Vec<_>>()
        .join(",\n");

    format!(
        "{{\n  \"schema_version\": \"1.0\",\n  \"run_id\": {},\n  \"created_unix_seconds\": {},\n  \"tool\": {{\"name\": \"splice-girl\", \"version\": {}}},\n  \"mode\": {},\n  \"command\": [{}],\n  \"inputs\": [\n{}\n  ],\n  \"outputs\": [{}],\n  \"settings\": {{\n    \"normalization\": \"count / total_gene_support\",\n    \"solver\": \"dense normal equations with one pinned node per connected component and partial pivoting\",\n    \"precision\": \"f64\",\n    \"acceptance_tolerance\": {:.1e},\n    \"minimum_total_support\": \"greater than zero\",\n    \"random_seed\": null\n  }},\n  \"summary\": {{\n    \"sample_gene_groups\": {},\n    \"junction_rows\": {},\n    \"included_samples\": {},\n    \"excluded_samples\": {},\n    \"failure_count\": {}\n  }},\n  \"exclusions\": [\n{}\n  ],\n  \"failures\": [\n{}\n  ]\n}}\n",
        json_string(&metadata.run_id),
        metadata.unix_time_seconds,
        json_string(env!("CARGO_PKG_VERSION")),
        json_string(metadata.mode),
        command,
        inputs,
        outputs,
        TOLERANCE,
        sample_gene_groups,
        junction_rows,
        included_samples,
        excluded_samples,
        failures.len(),
        exclusions_json,
        failures_json,
    )
}

fn tsv(value: &str) -> String {
    value
        .chars()
        .map(|character| match character {
            '\t' | '\n' | '\r' => ' ',
            other => other,
        })
        .collect()
}

fn json_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len() + 2);
    escaped.push('"');
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            character if character <= '\u{1f}' => {
                escaped.push_str(&format!("\\u{:04x}", character as u32));
            }
            character => escaped.push(character),
        }
    }
    escaped.push('"');
    escaped
}

fn sha256_file(path: &Path) -> Result<String, String> {
    let mut file = File::open(path)
        .map_err(|error| format!("could not open {} for checksum: {error}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let read = file
            .read(&mut buffer)
            .map_err(|error| format!("could not checksum {}: {error}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(hasher.finalize_hex())
}

/// Small streaming SHA-256 implementation so result bundles remain
/// dependency-free and work on systems without an external checksum command.
struct Sha256 {
    state: [u32; 8],
    buffer: [u8; 64],
    buffer_len: usize,
    message_len_bits: u64,
}

impl Sha256 {
    fn new() -> Self {
        Self {
            state: [
                0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
                0x5be0cd19,
            ],
            buffer: [0; 64],
            buffer_len: 0,
            message_len_bits: 0,
        }
    }

    fn update(&mut self, mut input: &[u8]) {
        self.message_len_bits = self
            .message_len_bits
            .wrapping_add((input.len() as u64).wrapping_mul(8));

        if self.buffer_len > 0 {
            let needed = 64 - self.buffer_len;
            let take = needed.min(input.len());
            self.buffer[self.buffer_len..self.buffer_len + take].copy_from_slice(&input[..take]);
            self.buffer_len += take;
            input = &input[take..];
            if self.buffer_len == 64 {
                let block = self.buffer;
                self.compress(&block);
                self.buffer_len = 0;
            }
        }

        while input.len() >= 64 {
            let block: &[u8; 64] = input[..64]
                .try_into()
                .expect("a 64-byte slice always converts to a 64-byte array");
            self.compress(block);
            input = &input[64..];
        }

        if !input.is_empty() {
            self.buffer[..input.len()].copy_from_slice(input);
            self.buffer_len = input.len();
        }
    }

    fn finalize_hex(mut self) -> String {
        let message_len_bits = self.message_len_bits;
        self.buffer[self.buffer_len] = 0x80;
        self.buffer_len += 1;

        if self.buffer_len > 56 {
            self.buffer[self.buffer_len..].fill(0);
            let block = self.buffer;
            self.compress(&block);
            self.buffer = [0; 64];
            self.buffer_len = 0;
        }

        self.buffer[self.buffer_len..56].fill(0);
        self.buffer[56..64].copy_from_slice(&message_len_bits.to_be_bytes());
        let block = self.buffer;
        self.compress(&block);

        self.state
            .iter()
            .map(|word| format!("{word:08x}"))
            .collect()
    }

    fn compress(&mut self, block: &[u8; 64]) {
        const K: [u32; 64] = [
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
            0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
            0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
            0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
            0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
            0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
            0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
            0xc67178f2,
        ];
        let mut words = [0_u32; 64];
        for (index, chunk) in block.chunks_exact(4).enumerate() {
            words[index] = u32::from_be_bytes(chunk.try_into().expect("chunk length is four"));
        }
        for index in 16..64 {
            let s0 = words[index - 15].rotate_right(7)
                ^ words[index - 15].rotate_right(18)
                ^ (words[index - 15] >> 3);
            let s1 = words[index - 2].rotate_right(17)
                ^ words[index - 2].rotate_right(19)
                ^ (words[index - 2] >> 10);
            words[index] = words[index - 16]
                .wrapping_add(s0)
                .wrapping_add(words[index - 7])
                .wrapping_add(s1);
        }

        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = self.state;
        for index in 0..64 {
            let big_s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let choose = (e & f) ^ ((!e) & g);
            let temp1 = h
                .wrapping_add(big_s1)
                .wrapping_add(choose)
                .wrapping_add(K[index])
                .wrapping_add(words[index]);
            let big_s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let majority = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = big_s0.wrapping_add(majority);

            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }

        for (state, value) in self.state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
            *state = state.wrapping_add(value);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analyze_junctions;
    use crate::compare::{compare_all_genes, compare_gene};
    use crate::manifest::parse_manifest_text;
    use crate::parser::parse_tsv_text;

    fn temporary_output(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock should be valid")
            .as_nanos();
        std::env::temp_dir().join(format!("splice-girl-{name}-{unique}"))
    }

    #[test]
    fn sha256_matches_known_vectors() {
        let mut empty = Sha256::new();
        empty.update(b"");
        assert_eq!(
            empty.finalize_hex(),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );

        let mut abc = Sha256::new();
        abc.update(b"a");
        abc.update(b"bc");
        assert_eq!(
            abc.finalize_hex(),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );

        let mut two_blocks = Sha256::new();
        two_blocks.update(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq");
        assert_eq!(
            two_blocks.finalize_hex(),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
    }

    #[test]
    fn analysis_bundle_contains_metrics_components_and_manifest() {
        let input = Path::new("tests/fixtures/cycle.tsv");
        let junctions = crate::parser::parse_tsv(input).expect("fixture should parse");
        let results = analyze_junctions(junctions).expect("fixture should analyze");
        let output = temporary_output("analysis-export");
        let command = vec![
            "splice-girl".to_owned(),
            "analyze".to_owned(),
            input.display().to_string(),
        ];

        let summary = write_analysis_bundle(&output, input, &results, &command)
            .expect("analysis export should succeed");
        assert_eq!(summary.files.len(), 3);
        let metrics = fs::read_to_string(output.join(SAMPLE_METRICS)).unwrap();
        assert!(metrics.contains("cycle_fraction"));
        assert!(metrics.contains("sample_cycle\tGENE_CYCLE"));
        let components = fs::read_to_string(output.join(JUNCTION_COMPONENTS)).unwrap();
        assert!(components.contains("raw_count\tnormalized_signal"));
        let manifest = fs::read_to_string(output.join(RUN_MANIFEST)).unwrap();
        assert!(manifest.contains("\"mode\": \"analyze\""));
        assert!(manifest.contains("\"sha256\""));
        assert!(manifest.contains(&summary.run_id));

        fs::remove_dir_all(output).unwrap();
    }

    #[test]
    fn gene_bundle_records_inclusions_and_exclusions() {
        let junction_text = fs::read_to_string("examples/demo/junctions.tsv").unwrap();
        let manifest_text = fs::read_to_string("examples/demo/manifest.tsv").unwrap();
        let results = analyze_junctions(parse_tsv_text(&junction_text).unwrap()).unwrap();
        let manifest = parse_manifest_text(&manifest_text).unwrap();
        let comparison = compare_gene(&manifest, &results, "GENE_DEMO");
        let output = temporary_output("gene-export");
        let command = vec!["splice-girl".to_owned(), "gene".to_owned()];

        let summary = write_gene_bundle(
            &output,
            Path::new("examples/demo/junctions.tsv"),
            Path::new("examples/demo/manifest.tsv"),
            &results,
            &comparison,
            &command,
        )
        .expect("gene export should succeed");
        assert_eq!(summary.files.len(), 4);
        let comparisons = fs::read_to_string(output.join(COMPARISONS)).unwrap();
        assert!(comparisons.contains("case_01\tcase\tincluded"));
        assert!(comparisons.contains("control_03\tcontrol\texcluded"));
        assert!(comparisons.contains("orphan_01\t\texcluded"));
        let manifest = fs::read_to_string(output.join(RUN_MANIFEST)).unwrap();
        assert!(manifest.contains("\"included_samples\": 4"));
        assert!(manifest.contains("\"excluded_samples\": 2"));

        fs::remove_dir_all(output).unwrap();
    }

    #[test]
    fn all_gene_bundle_records_descriptive_effects_and_no_inference() {
        let junction_text = fs::read_to_string("examples/demo/all_genes.tsv").unwrap();
        let manifest_text = fs::read_to_string("examples/demo/manifest.tsv").unwrap();
        let results = analyze_junctions(parse_tsv_text(&junction_text).unwrap()).unwrap();
        let manifest = parse_manifest_text(&manifest_text).unwrap();
        let summaries = compare_all_genes(&manifest, &results);
        let output = temporary_output("compare-export");
        let command = vec!["splice-girl".to_owned(), "compare".to_owned()];

        let summary = write_compare_bundle(
            &output,
            Path::new("examples/demo/all_genes.tsv"),
            Path::new("examples/demo/manifest.tsv"),
            &results,
            &summaries,
            &command,
        )
        .expect("all-gene export should succeed");

        assert_eq!(summary.files.len(), 5);
        let genes = fs::read_to_string(output.join(GENE_SUMMARIES)).unwrap();
        assert!(genes.contains("mean_difference_case_minus_control"));
        assert!(genes.contains("GENE_HIGH_CYCLE\tcomplete\t2\t2"));
        assert!(genes.contains("GENE_LOW_CYCLE\tcomplete\t2\t2"));
        assert!(genes.contains("0.5555555556"));
        assert!(genes.contains("-0.5555555556"));
        assert!(genes.contains("no p-value or q-value calculated"));
        let manifest = fs::read_to_string(output.join(RUN_MANIFEST)).unwrap();
        assert!(manifest.contains("\"mode\": \"compare\""));
        assert!(manifest.contains("gene_summaries.tsv"));

        fs::remove_dir_all(output).unwrap();
    }

    #[test]
    fn refuses_to_overwrite_an_existing_output_folder() {
        let output = temporary_output("existing-export");
        fs::create_dir_all(&output).unwrap();
        let input = Path::new("tests/fixtures/linear.tsv");
        let junctions = crate::parser::parse_tsv(input).unwrap();
        let results = analyze_junctions(junctions).unwrap();
        let error = write_analysis_bundle(&output, input, &results, &[]).unwrap_err();
        assert!(error.contains("already exists"));

        fs::remove_dir_all(output).unwrap();
    }
}
