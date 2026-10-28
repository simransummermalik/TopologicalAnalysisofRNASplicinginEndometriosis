use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

struct TestDirectory(PathBuf);

impl TestDirectory {
    fn new(label: &str) -> Self {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock should be valid")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "splice-girl-cli-workflows-{label}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&path).expect("temporary test directory should be created");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_splice-girl"))
        .args(arguments)
        .output()
        .expect("splice-girl command should start")
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn project_path(relative: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(relative)
        .display()
        .to_string()
}

#[test]
fn validate_reports_compact_input_and_math_summary() {
    let junctions = project_path("examples/demo/junctions.tsv");
    let output = run(&["validate", &junctions]);
    let stdout = text(&output.stdout);

    assert!(output.status.success(), "{}", text(&output.stderr));
    assert!(stdout.contains("validation summary"));
    assert!(stdout.contains("junction_rows: 12"));
    assert!(stdout.contains("samples: 5"));
    assert!(stdout.contains("mathematical_checks_failed: 0"));
    assert!(stdout.contains("validation_status: PASS"));
}

#[test]
fn validate_manifest_mismatches_are_visible_and_fail_the_audit() {
    let junctions = project_path("examples/demo/junctions.tsv");
    let manifest = project_path("examples/demo/manifest.tsv");
    let output = run(&["validate", &junctions, &manifest]);
    let stdout = text(&output.stdout);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stdout.contains("junction_samples_without_metadata: orphan_01"));
    assert!(stdout.contains("manifest_samples_without_junctions: control_03"));
    assert!(stdout.contains("validation_status: FAIL"));
    assert!(stderr.contains("sample/manifest mismatches"));
}

#[test]
fn compare_reports_every_gene_and_a_descriptive_effect() {
    let junctions = project_path("examples/demo/all_genes.tsv");
    let manifest = project_path("examples/demo/manifest.tsv");
    let output = run(&["compare", &junctions, &manifest]);
    let stdout = text(&output.stdout);

    assert!(output.status.success(), "{}", text(&output.stderr));
    assert!(stdout.contains("all-gene descriptive comparison"));
    assert!(stdout.contains("GENE_HIGH_CYCLE\tcomplete\t2\t2"));
    assert!(stdout.contains("GENE_LOW_CYCLE\tcomplete\t2\t2"));
    assert!(stdout.contains("0.5555555556"));
    assert!(stdout.contains("-0.5555555556"));
    assert!(stdout.contains("genes_with_both_groups: 2"));
    assert!(stdout.contains("no p-values, q-values, or ranking calculated"));
}

#[test]
fn compare_writes_a_reproducible_five_file_bundle() {
    let directory = TestDirectory::new("compare-output");
    let output_dir = directory.path().join("results");
    let junctions = project_path("examples/demo/all_genes.tsv");
    let manifest = project_path("examples/demo/manifest.tsv");
    let output = run(&[
        "compare",
        &junctions,
        &manifest,
        "--output",
        &output_dir.display().to_string(),
    ]);

    assert!(output.status.success(), "{}", text(&output.stderr));
    for name in [
        "sample_metrics.tsv",
        "junction_components.tsv",
        "comparisons.tsv",
        "gene_summaries.tsv",
        "run_manifest.json",
    ] {
        assert!(output_dir.join(name).is_file(), "missing {name}");
    }
}
