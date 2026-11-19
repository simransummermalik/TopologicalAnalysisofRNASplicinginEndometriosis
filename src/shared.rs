//! Native preparation of a common graph for every covered cohort sample.
//! These are explicit analysis choices, not implicit biological defaults.

mod cli;
mod sensitivity;

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::manifest::{Condition, SampleRecord, parse_manifest};
use crate::output::{json_string, sha256_file};
use crate::parser::Junction;

pub use cli::run_command;
pub use sensitivity::{SensitivityOptions, StatisticsSettings, run_sensitivity};

const HEADER: &str = "sample_id\tgene_id\tfrom\tto\tcount\n";
const MAX_EXACT_COUNT: u64 = 1 << 53;
type Edge = (String, String);
type EdgeCounts = BTreeMap<Edge, BTreeMap<String, u64>>;

/// Explicit support, recurrence and retained-coverage requirements.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Settings {
    pub minimum_edge_count: u64,
    pub minimum_samples_per_edge: usize,
    pub minimum_gene_support: u64,
}

impl Settings {
    fn validate(self, samples: usize) -> Result<(), String> {
        if self.minimum_edge_count == 0
            || self.minimum_samples_per_edge == 0
            || self.minimum_gene_support == 0
        {
            return Err("shared-graph thresholds must be positive integers".to_owned());
        }
        if self.minimum_samples_per_edge > samples {
            return Err("minimum samples per edge exceeds the manifest cohort size".to_owned());
        }
        Ok(())
    }

    fn json(self) -> String {
        format!(
            "{{\"minimum_edge_count\":{},\"minimum_samples_per_edge\":{},\"minimum_gene_support\":{},\"minimum_cycle_rank\":1,\"eligibility\":\"all manifest samples covered\"}}",
            self.minimum_edge_count, self.minimum_samples_per_edge, self.minimum_gene_support
        )
    }
}

/// Paths are resolved by the operating system relative to the caller's folder.
#[derive(Debug)]
pub struct Inputs {
    pub junctions: PathBuf,
    pub manifest: PathBuf,
    pub genes: Option<PathBuf>,
    /// A caller declaration after reviewing upstream complete count exports.
    pub complete_counts: bool,
}

struct Table {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn read(path: &Path, required: &[&str], exact: bool) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|e| format!("could not read {}: {e}", path.display()))?;
        let mut lines = text.lines();
        let header: Vec<String> = lines
            .next()
            .unwrap_or("")
            .split('\t')
            .map(str::to_owned)
            .collect();
        let unique: BTreeSet<_> = header.iter().collect();
        if unique.len() != header.len()
            || required
                .iter()
                .any(|name| !header.iter().any(|column| column == name))
            || (exact && header.iter().map(String::as_str).collect::<Vec<_>>() != required)
        {
            return Err(format!(
                "{}: expected unique columns {}",
                path.display(),
                required.join(", ")
            ));
        }
        let indices: Vec<_> = required
            .iter()
            .map(|name| header.iter().position(|value| value == name).unwrap())
            .collect();
        let mut rows = Vec::new();
        for (index, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let row: Vec<String> = line.split('\t').map(str::to_owned).collect();
            if row.len() != header.len() {
                return Err(format!(
                    "{} line {}: wrong number of fields",
                    path.display(),
                    index + 2
                ));
            }
            if indices
                .iter()
                .any(|&i| row[i].trim().is_empty() || row[i] != row[i].trim())
            {
                return Err(format!(
                    "{} line {}: empty or padded required value",
                    path.display(),
                    index + 2
                ));
            }
            rows.push(row);
        }
        Ok(Self { header, rows })
    }

    fn column(&self, name: &str) -> usize {
        self.header
            .iter()
            .position(|value| value == name)
            .expect("validated table column")
    }
}

/// Loaded counts and family metadata, shared across sensitivity scenarios.
pub struct Cohort {
    pub samples: Vec<SampleRecord>,
    counts: BTreeMap<String, EdgeCounts>,
    family_header: Vec<String>,
    family: BTreeMap<String, Vec<String>>,
    manifest_text: String,
}

impl Cohort {
    pub fn load(inputs: &Inputs) -> Result<Self, String> {
        if !inputs.complete_counts {
            return Err("--complete-counts is required: verify complete exports for every sample with minimum junction support 1 before inserting measured zeros".to_owned());
        }
        Table::read(&inputs.manifest, &["sample_id", "condition"], false)?;
        let mut samples = parse_manifest(&inputs.manifest)?;
        samples.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));
        if ![Condition::Case, Condition::Control]
            .iter()
            .all(|condition| samples.iter().any(|sample| sample.condition == *condition))
        {
            return Err("manifest must contain both case and control samples".to_owned());
        }
        let ids: BTreeSet<_> = samples
            .iter()
            .map(|sample| sample.sample_id.as_str())
            .collect();
        let raw = Table::read(
            &inputs.junctions,
            &["sample_id", "gene_id", "from", "to", "count"],
            true,
        )?;
        let mut observed = BTreeSet::new();
        let mut counts: BTreeMap<String, EdgeCounts> = BTreeMap::new();
        for (index, row) in raw.rows.into_iter().enumerate() {
            if !ids.contains(row[0].as_str()) {
                return Err(format!(
                    "junction line {}: sample {:?} is absent from manifest",
                    index + 2,
                    row[0]
                ));
            }
            if row[2] == row[3] {
                return Err(format!(
                    "junction line {}: self-loop is not supported",
                    index + 2
                ));
            }
            let count = if !row[4].is_empty() && row[4].bytes().all(|byte| byte.is_ascii_digit()) {
                row[4].parse::<u64>().ok().filter(|value| *value <= MAX_EXACT_COUNT)
            } else { None }.ok_or_else(|| format!("junction line {}: count must be a nonnegative integer within the exact f64 range (0..={MAX_EXACT_COUNT})", index + 2))?;
            if count > 0 {
                observed.insert(row[0].clone());
            }
            if counts
                .entry(row[1].clone())
                .or_default()
                .entry((row[2].clone(), row[3].clone()))
                .or_default()
                .insert(row[0].clone(), count)
                .is_some()
            {
                return Err(format!(
                    "junction line {}: duplicate sample/gene edge",
                    index + 2
                ));
            }
        }
        let absent: Vec<_> = ids
            .iter()
            .filter(|id| !observed.contains(**id))
            .copied()
            .collect();
        if !absent.is_empty() {
            return Err(format!(
                "manifest samples have no positive input counts: {}; check missing/failed exports",
                absent.join(", ")
            ));
        }
        let (family_header, family) = if let Some(path) = &inputs.genes {
            let table = Table::read(path, &["gene_id"], false)?;
            crate::statistics::parse_gene_family(path)?;
            let index = table.column("gene_id");
            let family = table
                .rows
                .into_iter()
                .map(|row| (row[index].clone(), row))
                .collect();
            (table.header, family)
        } else {
            (
                vec!["gene_id".to_owned()],
                counts
                    .keys()
                    .map(|gene| (gene.clone(), vec![gene.clone()]))
                    .collect(),
            )
        };
        let manifest_text = fs::read_to_string(&inputs.manifest).map_err(|e| e.to_string())?;
        Ok(Self {
            samples,
            counts,
            family_header,
            family,
            manifest_text,
        })
    }

    pub fn build(&self, settings: Settings) -> Result<Prepared, String> {
        settings.validate(self.samples.len())?;
        let mut prepared = Prepared {
            junctions: Vec::new(), genes: Vec::new(), eligible_genes: Vec::new(),
            edges: "gene_id\tfrom\tto\tsamples_meeting_support\tselected\n".to_owned(),
            sample_audit: "gene_id\tsample_id\tshared_graph_support\tcoverage_pass\tcycle_rank\tzero_edges\tstatus\n".to_owned(),
            family: format!("{}\n", self.family_header.join("\t")),
        };
        let empty = BTreeMap::new();
        for gene in self.family.keys() {
            let counts = self.counts.get(gene).unwrap_or(&empty);
            let mut edges = Vec::new();
            for (edge, values) in counts {
                let supporters = values
                    .values()
                    .filter(|&&value| value >= settings.minimum_edge_count)
                    .count();
                let selected = supporters >= settings.minimum_samples_per_edge;
                prepared.edges.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    gene, edge.0, edge.1, supporters, selected
                ));
                if selected {
                    edges.push(edge);
                }
            }
            let rank = cycle_rank(&edges);
            let mut covered_samples = 0;
            for sample in &self.samples {
                let values: Vec<_> = edges
                    .iter()
                    .map(|edge| counts[*edge].get(&sample.sample_id).copied().unwrap_or(0))
                    .collect();
                let total = values
                    .iter()
                    .try_fold(0u64, |sum, &value| {
                        sum.checked_add(value)
                            .filter(|total| *total <= MAX_EXACT_COUNT)
                    })
                    .ok_or_else(|| {
                        format!(
                            "{}/{gene}: retained support exceeds exact f64 integer range",
                            sample.sample_id
                        )
                    })?;
                let covered = total >= settings.minimum_gene_support;
                covered_samples += usize::from(covered);
                let status = if edges.is_empty() {
                    "no_shared_edges"
                } else if rank == 0 {
                    "no_cycle_in_shared_graph"
                } else if !covered {
                    "below_coverage"
                } else {
                    "included"
                };
                prepared.sample_audit.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    gene,
                    sample.sample_id,
                    total,
                    covered,
                    rank,
                    values.iter().filter(|&&value| value == 0).count(),
                    status
                ));
                if status == "included" {
                    for (edge, value) in edges.iter().zip(values) {
                        prepared.junctions.push(Junction {
                            sample_id: sample.sample_id.clone(),
                            gene_id: gene.clone(),
                            from: edge.0.clone(),
                            to: edge.1.clone(),
                            count: value as f64,
                        });
                    }
                }
            }
            let status = if edges.is_empty() {
                "no_shared_edges"
            } else if rank == 0 {
                "no_cycle_in_shared_graph"
            } else if covered_samples != self.samples.len() {
                "incomplete_coverage"
            } else {
                "eligible"
            };
            if status == "eligible" {
                prepared.eligible_genes.push(gene.clone());
                prepared
                    .family
                    .push_str(&format!("{}\n", self.family[gene].join("\t")));
            }
            prepared.genes.push(GeneEligibility {
                gene_id: gene.clone(),
                shared_edges: edges.len(),
                cycle_rank: rank,
                covered_samples,
                required_samples: self.samples.len(),
                status,
            });
        }
        Ok(prepared)
    }
}

fn cycle_rank(edges: &[&Edge]) -> usize {
    let names: BTreeSet<_> = edges.iter().flat_map(|edge| [&edge.0, &edge.1]).collect();
    let indices: BTreeMap<_, _> = names
        .into_iter()
        .enumerate()
        .map(|(i, name)| (name, i))
        .collect();
    let mut parent: Vec<_> = (0..indices.len()).collect();
    fn root(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }
    for edge in edges {
        let a = root(&mut parent, indices[&edge.0]);
        let b = root(&mut parent, indices[&edge.1]);
        parent[a] = b;
    }
    let components: BTreeSet<_> = (0..parent.len())
        .map(|node| root(&mut parent, node))
        .collect();
    edges.len() + components.len() - indices.len()
}

pub struct GeneEligibility {
    pub gene_id: String,
    pub shared_edges: usize,
    pub cycle_rank: usize,
    pub covered_samples: usize,
    pub required_samples: usize,
    pub status: &'static str,
}

/// Validated Rust analysis rows plus explicit edge/sample/gene decisions.
pub struct Prepared {
    pub junctions: Vec<Junction>,
    pub genes: Vec<GeneEligibility>,
    pub eligible_genes: Vec<String>,
    edges: String,
    sample_audit: String,
    family: String,
}

/// Prepare a folder using only the installed executable and supplied inputs.
pub fn prepare(
    inputs: &Inputs,
    settings: Settings,
    output: &Path,
    command: &[String],
) -> Result<usize, String> {
    let cohort = Cohort::load(inputs)?;
    let prepared = cohort.build(settings)?;
    export_prepared(inputs, &cohort, &prepared, settings, output, command)?;
    Ok(prepared.eligible_genes.len())
}

fn export_prepared(
    inputs: &Inputs,
    cohort: &Cohort,
    prepared: &Prepared,
    settings: Settings,
    output: &Path,
    command: &[String],
) -> Result<(), String> {
    let mut run = NativeRun::new(
        output,
        "shared graph preparation",
        inputs,
        &settings.json(),
        command,
    )?;
    let result = (|| {
        let mut junctions = HEADER.to_owned();
        for row in &prepared.junctions {
            junctions.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.0}\n",
                row.sample_id, row.gene_id, row.from, row.to, row.count
            ));
        }
        let mut genes =
            "gene_id\tshared_edges\tcycle_rank\tcovered_samples\trequired_samples\tstatus\n"
                .to_owned();
        for row in &prepared.genes {
            genes.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                row.gene_id,
                row.shared_edges,
                row.cycle_rank,
                row.covered_samples,
                row.required_samples,
                row.status
            ));
        }
        for (name, text) in [
            ("junctions.tsv", junctions.as_str()),
            ("edges.tsv", &prepared.edges),
            ("sample_eligibility.tsv", &prepared.sample_audit),
            ("gene_eligibility.tsv", &genes),
            ("eligible_genes.tsv", &prepared.family),
            ("samples.tsv", &cohort.manifest_text),
        ] {
            run.write(name, text)?;
        }
        run.summary = format!(
            "{{\"eligible_genes\":{},\"retained_rows\":{}}}",
            prepared.eligible_genes.len(),
            prepared.junctions.len()
        );
        Ok(())
    })();
    run.finish(result)
}

/// A claimed output directory with checksummed provenance, including failures.
struct NativeRun {
    root: PathBuf,
    base: String,
    files: BTreeSet<PathBuf>,
    scenarios: Vec<String>,
    summary: String,
}

impl NativeRun {
    fn new(
        output: &Path,
        mode: &str,
        inputs: &Inputs,
        settings: &str,
        command: &[String],
    ) -> Result<Self, String> {
        if output.exists() || output.symlink_metadata().is_ok() {
            return Err(format!(
                "output path {} already exists; choose a new folder",
                output.display()
            ));
        }
        let mut paths = vec![&inputs.junctions, &inputs.manifest];
        paths.extend(inputs.genes.iter());
        let records: Result<Vec<_>, _> = paths
            .into_iter()
            .map(|path| file_record(path, &path.display().to_string()))
            .collect();
        let executable =
            std::env::current_exe().map_err(|e| format!("could not identify executable: {e}"))?;
        let base = format!(
            "\"schema_version\":\"1.0\",\"tool\":{{\"name\":\"splice-girl\",\"version\":{},\"executable_sha256\":{}}},\"implementation\":\"native Rust\",\"mode\":{},\"command\":[{}],\"inputs\":[{}],\"settings\":{},\"complete_counts_declared\":true,\"approval\":\"explicit caller settings; not team approval\",\"anchor_sensitivity\":\"not evaluated; requires separate read-counter exports\"",
            json_string(env!("CARGO_PKG_VERSION")),
            json_string(&sha256_file(&executable)?),
            json_string(mode),
            command
                .iter()
                .map(|value| json_string(value))
                .collect::<Vec<_>>()
                .join(","),
            records?.join(","),
            settings
        );
        if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
            fs::create_dir_all(parent)
                .map_err(|e| format!("could not create {}: {e}", parent.display()))?;
        }
        fs::create_dir(output)
            .map_err(|e| format!("could not create {}: {e}", output.display()))?;
        let run = Self {
            root: output.to_path_buf(),
            base,
            files: BTreeSet::new(),
            scenarios: Vec::new(),
            summary: "null".to_owned(),
        };
        run.save("running", None)?;
        Ok(run)
    }

    fn write(&mut self, name: &str, text: &str) -> Result<(), String> {
        let path = self.root.join(name);
        self.files.insert(PathBuf::from(name));
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .and_then(|mut file| file.write_all(text.as_bytes()))
            .map_err(|e| format!("could not write {}: {e}", path.display()))
    }

    fn save(&self, status: &str, error: Option<&str>) -> Result<(), String> {
        let files: Result<Vec<_>, _> = self
            .files
            .iter()
            .filter(|name| self.root.join(name).is_file())
            .map(|name| file_record(&self.root.join(name), &name.display().to_string()))
            .collect();
        let json = format!(
            "{{{},\"status\":{},\"error\":{},\"summary\":{},\"scenarios\":[{}],\"outputs\":[{}]}}\n",
            self.base,
            json_string(status),
            error.map(json_string).unwrap_or_else(|| "null".to_owned()),
            self.summary,
            self.scenarios.join(","),
            files?.join(",")
        );
        fs::write(self.root.join("run_manifest.json"), json)
            .map_err(|e| format!("could not save run manifest: {e}"))
    }

    fn finish(&self, result: Result<(), String>) -> Result<(), String> {
        match result {
            Ok(()) => self.save("complete", None),
            Err(error) => {
                if let Err(record_error) = self.save("failed", Some(&error)) {
                    return Err(format!("{error}; additionally {record_error}"));
                }
                Err(error)
            }
        }
    }
}

fn file_record(path: &Path, label: &str) -> Result<String, String> {
    Ok(format!(
        "{{\"path\":{},\"sha256\":{}}}",
        json_string(label),
        json_string(&sha256_file(path)?)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn native_write_failure_preserves_files_and_records_failed_manifest() {
        let id = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "splice-girl-native-failure-{}-{id}",
            std::process::id()
        ));
        let input_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/shared_graph");
        let inputs = Inputs {
            junctions: input_root.join("junctions.tsv"),
            manifest: input_root.join("samples.tsv"),
            genes: Some(input_root.join("genes.tsv")),
            complete_counts: true,
        };
        let mut run = NativeRun::new(&root, "failure regression", &inputs, "{}", &[]).unwrap();
        run.write("preserved.tsv", "original\n").unwrap();
        let result = run.write("preserved.tsv", "replacement\n");
        assert!(result.is_err());
        assert!(run.finish(result).is_err());
        let manifest = fs::read_to_string(root.join("run_manifest.json")).unwrap();
        assert!(manifest.contains("\"status\":\"failed\""));
        assert!(!manifest.contains("\"status\":\"complete\""));
        assert_eq!(
            fs::read_to_string(root.join("preserved.tsv")).unwrap(),
            "original\n"
        );
        fs::remove_dir_all(root).unwrap();
    }
}
