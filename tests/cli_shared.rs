//! Installed-style acceptance: copy only the binary to a new working folder,
//! remove PATH, then use embedded example inputs and native analysis commands.
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

struct Workspace {
    root: PathBuf,
    binary: PathBuf,
}

static WORKSPACE_NUMBER: AtomicUsize = AtomicUsize::new(0);

impl Workspace {
    fn new() -> Self {
        let id = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let number = WORKSPACE_NUMBER.fetch_add(1, Ordering::Relaxed);
        let root = std::env::temp_dir().join(format!(
            "splice-girl-native-{}-{id}-{number}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        let binary = root.join("splice-girl");
        fs::copy(env!("CARGO_BIN_EXE_splice-girl"), &binary).unwrap();
        Self { root, binary }
    }

    fn run(&self, args: &[&str]) -> Output {
        Command::new(&self.binary)
            .args(args)
            .current_dir(&self.root)
            .env("PATH", "")
            .output()
            .unwrap()
    }

    fn ok(&self, args: &[&str]) {
        let result = self.run(args);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    fn example(&self) {
        self.ok(&["example", "--output", "inputs"]);
    }

    fn text(&self, file: &str) -> String {
        fs::read_to_string(self.root.join(file)).unwrap()
    }

    fn prepare_args(output: &str) -> Vec<&str> {
        vec![
            "prepare-shared",
            "--junctions",
            "inputs/junctions.tsv",
            "--manifest",
            "inputs/samples.tsv",
            "--genes",
            "inputs/genes.tsv",
            "--complete-counts",
            "--minimum-edge-count",
            "5",
            "--minimum-samples-per-edge",
            "4",
            "--minimum-gene-support",
            "50",
            "--output",
            output,
        ]
    }

    fn sensitivity_args(output: &str) -> Vec<&str> {
        vec![
            "sensitivity",
            "--junctions",
            "inputs/junctions.tsv",
            "--manifest",
            "inputs/samples.tsv",
            "--genes",
            "inputs/genes.tsv",
            "--complete-counts",
            "--primary-edge-count",
            "5",
            "--primary-sample-count",
            "4",
            "--primary-gene-support",
            "50",
            "--minimum-junction-counts",
            "3,5,10",
            "--minimum-samples-per-edge",
            "2,4,6",
            "--minimum-gene-support",
            "20,50,100",
            "--run-stats",
            "--family",
            "synthetic-native",
            "--alternative",
            "two-sided",
            "--design",
            "independent",
            "--min-per-group",
            "4",
            "--output",
            output,
        ]
    }
}

impl Drop for Workspace {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn table_value<'a>(table: &'a str, column: &str, gene: &str) -> &'a str {
    let mut lines = table.lines();
    let header: Vec<_> = lines.next().unwrap().split('\t').collect();
    let gene_index = header.iter().position(|&name| name == "gene_id").unwrap();
    let index = header.iter().position(|&name| name == column).unwrap();
    lines
        .find_map(|line| {
            let row: Vec<_> = line.split('\t').collect();
            (row[gene_index] == gene).then(|| row[index])
        })
        .unwrap()
}

#[test]
fn standalone_example_help_and_output_protection_need_no_repository_or_path() {
    let work = Workspace::new();
    for command in ["example", "prepare-shared", "sensitivity"] {
        work.ok(&[command, "--help"]);
    }
    work.example();
    assert!(
        work.text("inputs/README.md")
            .contains("splice-girl prepare-shared")
    );
    assert!(
        work.text("inputs/run_manifest.json")
            .contains("\"status\":\"complete\"")
    );
    let original = work.text("inputs/junctions.tsv");
    assert!(
        !work
            .run(&["example", "--output", "inputs"])
            .status
            .success()
    );
    assert_eq!(work.text("inputs/junctions.tsv"), original);
}

#[test]
fn native_preparation_keeps_zero_edges_and_excludes_missing_and_acyclic_genes() {
    let work = Workspace::new();
    work.example();
    work.ok(&Workspace::prepare_args("prepared"));
    let audit = work.text("prepared/gene_eligibility.tsv");
    for (gene, expected) in [
        ("TRIANGLE", "eligible"),
        ("TREE", "no_cycle_in_shared_graph"),
        ("LOW_COVERAGE", "incomplete_coverage"),
        ("ABSENT", "no_shared_edges"),
    ] {
        assert_eq!(table_value(&audit, "status", gene), expected);
    }
    let rows = work.text("prepared/junctions.tsv");
    assert!(rows.contains("S7\tTRIANGLE\tA\tC\t0\n"));
    assert!(!rows.contains("S7\tLOW_COVERAGE"));
    let family = work.text("prepared/eligible_genes.tsv");
    assert_eq!(
        family,
        "gene_id\tgene_symbol\nTRIANGLE\tSynthetic_triangle\n"
    );
    work.ok(&[
        "stats",
        "prepared/junctions.tsv",
        "prepared/samples.tsv",
        "--genes",
        "prepared/eligible_genes.tsv",
        "--family",
        "synthetic-native",
        "--alternative",
        "two-sided",
        "--design",
        "independent",
        "--min-per-group",
        "4",
        "--output",
        "statistics",
    ]);
    let statistics = work.text("statistics/statistics.tsv");
    assert_eq!(table_value(&statistics, "status", "TRIANGLE"), "tested");
    assert!(
        (table_value(&statistics, "p_value", "TRIANGLE")
            .parse::<f64>()
            .unwrap()
            - 2.0 / 70.0)
            .abs()
            < 1e-9
    );
    let manifest = work.text("prepared/run_manifest.json");
    assert!(manifest.contains("\"implementation\":\"native Rust\""));
    assert!(manifest.contains("\"status\":\"complete\""));
    let before = work.text("prepared/junctions.tsv");
    assert!(
        !work
            .run(&Workspace::prepare_args("prepared"))
            .status
            .success()
    );
    assert_eq!(work.text("prepared/junctions.tsv"), before);
}

#[test]
fn native_sensitivity_runs_seven_scenarios_and_skips_all_excluded_statistics() {
    let work = Workspace::new();
    work.example();
    work.ok(&Workspace::sensitivity_args("sensitivity"));
    let table = work.text("sensitivity/sensitivity.tsv");
    assert_eq!(table.lines().count(), 1 + 7 * 4);
    let primary = table
        .lines()
        .filter(|line| line.starts_with("primary\t"))
        .collect::<Vec<_>>()
        .join("\n");
    assert!(primary.contains("TRIANGLE\teligible\t3\t1\t8\t8\t4\t4\t-0.1141552511"));
    let record = work.text("sensitivity/run_manifest.json");
    assert_eq!(
        record
            .matches("\"statistics_status\":\"no_eligible_genes\"")
            .count(),
        3
    );
    assert!(record.contains("\"status\":\"complete\""));
    assert!(
        work.root
            .join("sensitivity/primary/statistics/report.html")
            .is_file()
    );
    assert!(
        !work
            .root
            .join("sensitivity/scenario_006/statistics")
            .exists()
    );
    for line in table.lines().skip(1) {
        let row: Vec<_> = line.split('\t').collect();
        if row[5] != "eligible" {
            assert_eq!(row[12], "");
        }
    }
    let metrics = work.text("sensitivity/primary/results/sample_metrics.tsv");
    let control = metrics
        .lines()
        .find(|line| line.contains("\tS7\tTRIANGLE\t"))
        .unwrap();
    let header: Vec<_> = metrics.lines().next().unwrap().split('\t').collect();
    let fraction = control
        .split('\t')
        .nth(
            header
                .iter()
                .position(|&name| name == "cycle_fraction")
                .unwrap(),
        )
        .unwrap()
        .parse::<f64>()
        .unwrap();
    assert!((fraction - 2.0 / 3.0).abs() < 1e-9);
}

#[test]
fn edge_selection_is_invariant_to_labels_and_input_row_order() {
    let work = Workspace::new();
    work.example();
    work.ok(&Workspace::prepare_args("first"));
    let raw = work.text("inputs/junctions.tsv");
    let mut rows: Vec<_> = raw.lines().skip(1).collect();
    rows.reverse();
    fs::write(
        work.root.join("inputs/junctions.tsv"),
        format!("{}\n{}\n", raw.lines().next().unwrap(), rows.join("\n")),
    )
    .unwrap();
    let manifest = work
        .text("inputs/samples.tsv")
        .replace("case", "temporary")
        .replace("control", "case")
        .replace("temporary", "control");
    fs::write(work.root.join("inputs/samples.tsv"), manifest).unwrap();
    work.ok(&Workspace::prepare_args("second"));
    for name in [
        "junctions.tsv",
        "edges.tsv",
        "sample_eligibility.tsv",
        "gene_eligibility.tsv",
        "eligible_genes.tsv",
    ] {
        assert_eq!(
            work.text(&format!("first/{name}")),
            work.text(&format!("second/{name}"))
        );
    }
}

#[test]
fn bad_options_missing_exports_and_invalid_counts_fail_before_creating_results() {
    let work = Workspace::new();
    work.example();
    let args = Workspace::prepare_args("bad");
    let mut missing = args.clone();
    missing.retain(|arg| *arg != "--complete-counts");
    let failure = work.run(&missing);
    assert!(!failure.status.success());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("--complete-counts"));
    assert!(!work.root.join("bad").exists());
    for count in ["-1", "NaN", "0.5", "9007199254740993"] {
        let raw = include_str!("../examples/shared_graph/junctions.tsv").replacen(
            "\t30\n",
            &format!("\t{count}\n"),
            1,
        );
        fs::write(work.root.join("inputs/junctions.tsv"), raw).unwrap();
        assert!(!work.run(&args).status.success());
        assert!(!work.root.join("bad").exists());
    }
    let raw = include_str!("../examples/shared_graph/junctions.tsv");
    fs::write(
        work.root.join("inputs/junctions.tsv"),
        raw.lines()
            .filter(|line| !line.starts_with("S7\t"))
            .collect::<Vec<_>>()
            .join("\n")
            + "\n",
    )
    .unwrap();
    let failure = work.run(&args);
    assert!(String::from_utf8_lossy(&failure.stderr).contains("no positive input counts"));
    fs::write(work.root.join("inputs/junctions.tsv"), raw).unwrap();
    let mut duplicate = args.clone();
    duplicate.extend(["--minimum-edge-count", "5"]);
    assert!(String::from_utf8_lossy(&work.run(&duplicate).stderr).contains("duplicate"));
    let mut invalid = Workspace::sensitivity_args("bad");
    let grid = invalid.iter().position(|&value| value == "3,5,10").unwrap();
    invalid[grid] = "3,nan";
    assert!(!work.run(&invalid).status.success());
    assert!(!work.root.join("bad").exists());
}

#[test]
fn independent_subject_validation_happens_before_any_sensitivity_export() {
    let work = Workspace::new();
    work.example();
    let manifest = work
        .text("inputs/samples.tsv")
        .replace("S7\tcontrol\tS7", "S7\tcontrol\tS6");
    fs::write(work.root.join("inputs/samples.tsv"), manifest).unwrap();
    let failure = work.run(&Workspace::sensitivity_args("bad"));
    assert!(!failure.status.success());
    assert!(String::from_utf8_lossy(&failure.stderr).contains("subject_id"));
    assert!(!work.root.join("bad").exists());
}

#[test]
fn duplicate_family_columns_are_rejected_before_sensitivity_export() {
    let work = Workspace::new();
    work.example();
    fs::write(
        work.root.join("inputs/genes.tsv"),
        "gene_id\tcategory\tcategory\nTRIANGLE\ta\tb\n",
    )
    .unwrap();
    let bad = work.run(&Workspace::sensitivity_args("bad"));
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("unique columns"));
    assert!(!work.root.join("bad").exists());
}

#[test]
fn optional_family_and_descriptive_sensitivity_work_without_statistics() {
    let work = Workspace::new();
    work.example();
    let mut args = Workspace::prepare_args("observed");
    let index = args.iter().position(|&arg| arg == "--genes").unwrap();
    args.drain(index..index + 2);
    work.ok(&args);
    assert!(
        !work
            .text("observed/gene_eligibility.tsv")
            .contains("ABSENT")
    );
    let mut args = Workspace::sensitivity_args("descriptive");
    let index = args.iter().position(|&arg| arg == "--run-stats").unwrap();
    let output = args.split_off(args.len() - 2);
    args.truncate(index);
    args.extend(output);
    work.ok(&args);
    let manifest = work.text("descriptive/run_manifest.json");
    assert_eq!(
        manifest
            .matches("\"statistics_status\":\"not_requested\"")
            .count(),
        7
    );
    assert!(!work.root.join("descriptive/primary/statistics").exists());
}
