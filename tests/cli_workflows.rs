use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
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

#[test]
fn report_exports_the_same_values_as_compare_with_linked_graphs() {
    let directory = TestDirectory::new("report");
    let output_dir = directory.path().join("nested/report");
    let comparison_dir = directory.path().join("comparison");
    let junctions = project_path("examples/demo/all_genes.tsv");
    let manifest = project_path("examples/demo/manifest.tsv");
    let output = run(&[
        "report",
        &junctions,
        &manifest,
        "--output",
        output_dir.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{}", text(&output.stderr));
    let report = fs::read_to_string(output_dir.join("report.html")).unwrap();
    assert!(report.contains("GENE_HIGH_CYCLE"));
    assert!(report.contains("GENE_LOW_CYCLE"));
    assert!(report.contains("0.5555555556"));
    assert!(report.contains("-0.5555555556"));
    assert!(report.contains("control_03"));
    assert_eq!(report.matches("View decomposition</a>").count(), 8);
    assert_eq!(fs::read_dir(output_dir.join("graphs")).unwrap().count(), 8);
    for index in 1..=8 {
        let name = format!("graphs/graph-{index:06}.svg");
        assert!(report.contains(&format!("href=\"{name}\"")));
        let svg = fs::read_to_string(output_dir.join(&name)).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("3. Cycle component"));
    }
    let provenance = fs::read_to_string(output_dir.join("run_manifest.json")).unwrap();
    assert!(provenance.contains("\"mode\": \"report\""));
    assert!(provenance.contains("report.html"));
    assert!(provenance.contains("graphs/graph-000008.svg"));
    let comparison = run(&[
        "compare",
        &junctions,
        &manifest,
        "--output",
        comparison_dir.to_str().unwrap(),
    ]);
    assert!(comparison.status.success());
    for name in [
        "sample_metrics.tsv",
        "junction_components.tsv",
        "comparisons.tsv",
        "gene_summaries.tsv",
    ] {
        let without_run_id = |path| {
            fs::read_to_string(path)
                .unwrap()
                .lines()
                .map(|line| line.split_once('\t').unwrap().1.to_owned())
                .collect::<Vec<_>>()
        };
        assert_eq!(
            without_run_id(output_dir.join(name)),
            without_run_id(comparison_dir.join(name))
        );
    }
    let existing = run(&[
        "report",
        &junctions,
        &manifest,
        "--output",
        output_dir.to_str().unwrap(),
    ]);
    assert!(!existing.status.success());
    assert!(text(&existing.stderr).contains("already exists"));
    assert_eq!(
        report,
        fs::read_to_string(output_dir.join("report.html")).unwrap()
    );
}

#[test]
fn report_requires_an_explicit_destination() {
    let output = run(&[
        "report",
        &project_path("examples/demo/all_genes.tsv"),
        &project_path("examples/demo/manifest.tsv"),
    ]);
    assert!(!output.status.success());
    assert!(text(&output.stderr).contains("report needs --output"));
}

#[test]
fn report_keeps_missing_groups_visible_in_a_diagnostic_bundle() {
    let directory = TestDirectory::new("report-missing");
    let input = directory.path().join("input.tsv");
    let manifest = directory.path().join("manifest.tsv");
    fs::write(
        &input,
        "sample_id\tgene_id\tfrom\tto\tcount\ns1\tGENE\tE1\tE2\t4\n",
    )
    .unwrap();
    fs::write(&manifest, "sample_id\tcondition\ns1\tcase\ns2\tcontrol\n").unwrap();
    let output_dir = directory.path().join("results");
    let output = run(&[
        "report",
        input.to_str().unwrap(),
        manifest.to_str().unwrap(),
        "--output",
        output_dir.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(text(&output.stderr).contains("no gene has at least one eligible case"));
    let report = fs::read_to_string(output_dir.join("report.html")).unwrap();
    assert!(report.contains("missing_control"));
    assert!(report.contains("<strong>NA</strong>"));
    assert!(report.contains("no junction data for gene"));
}

#[test]
fn report_does_not_use_identifiers_as_paths_or_markup() {
    let directory = TestDirectory::new("report-escape");
    let input = directory.path().join("input.tsv");
    let manifest = directory.path().join("manifest.tsv");
    fs::write(&input, "sample_id\tgene_id\tfrom\tto\tcount\n<script>\t../bad&gene\tE1\tE2\t4\ncontrol\t../bad&gene\tE1\tE2\t4\n").unwrap();
    fs::write(
        &manifest,
        "sample_id\tcondition\n<script>\tcase\ncontrol\tcontrol\n",
    )
    .unwrap();
    let output_dir = directory.path().join("results");
    let output = run(&[
        "report",
        input.to_str().unwrap(),
        manifest.to_str().unwrap(),
        "--output",
        output_dir.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{}", text(&output.stderr));
    let report = fs::read_to_string(output_dir.join("report.html")).unwrap();
    assert!(report.contains("&lt;script&gt;"));
    assert!(report.contains("../bad&amp;gene"));
    assert!(!report.contains("<script>"));
    assert!(output_dir.join("graphs/graph-000001.svg").is_file());
    assert!(!directory.path().join("bad&gene").exists());
}

#[test]
fn interactive_menu_can_generate_a_report_then_return_to_menu() {
    let directory = TestDirectory::new("report-menu");
    let output_dir = directory.path().join("results");
    let mut child = Command::new(env!("CARGO_BIN_EXE_splice-girl"))
        .arg("--start")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let input = format!(
        "6\n{}\n{}\n{}\n\nq\n",
        project_path("examples/demo/all_genes.tsv"),
        project_path("examples/demo/manifest.tsv"),
        output_dir.display()
    );
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "{}", text(&output.stderr));
    assert!(output_dir.join("report.html").is_file());
    assert!(text(&output.stdout).contains("Open in your browser:"));
    assert!(!text(&output.stdout).contains("  error:"));
}
