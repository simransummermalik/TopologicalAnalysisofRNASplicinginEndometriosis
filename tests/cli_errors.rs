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
            "splice-girl-cli-errors-{label}-{}-{unique}",
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

fn write(path: &Path, contents: &str) {
    fs::write(path, contents).expect("test input should be written");
}

#[test]
fn missing_input_names_the_path_and_suggests_a_fix() {
    let directory = TestDirectory::new("missing");
    let missing = directory.path().join("missing.tsv");
    let output = run(&["analyze", &missing.display().to_string()]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&missing.display().to_string()));
    assert!(stderr.contains("check that the path exists"));
}

#[test]
fn empty_input_names_the_file_and_required_header() {
    let directory = TestDirectory::new("empty");
    let input = directory.path().join("empty.tsv");
    write(&input, "");
    let output = run(&["analyze", &input.display().to_string()]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&input.display().to_string()));
    assert!(stderr.contains("file is empty"));
    assert!(stderr.contains("sample_id\\tgene_id\\tfrom\\tto\\tcount"));
}

#[test]
fn malformed_header_names_the_file_and_expected_columns() {
    let directory = TestDirectory::new("header");
    let input = directory.path().join("bad_header.tsv");
    write(&input, "sample,gene,from,to,count\ns1,g1,E1,E2,1\n");
    let output = run(&["analyze", &input.display().to_string()]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&input.display().to_string()));
    assert!(stderr.contains("header must be exactly"));
    assert!(stderr.contains("check the column names, order, and tab separators"));
}

#[test]
fn invalid_count_names_the_file_line_column_value_and_fix() {
    let directory = TestDirectory::new("count");
    let input = directory.path().join("bad_count.tsv");
    write(
        &input,
        "sample_id\tgene_id\tfrom\tto\tcount\ns1\tg1\tE1\tE2\tmany\n",
    );
    let output = run(&["analyze", &input.display().to_string()]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&input.display().to_string()));
    assert!(stderr.contains("line 2, column 'count'"));
    assert!(stderr.contains("found 'many'"));
    assert!(stderr.contains("replace it with a numeric read count"));
}

#[test]
fn invalid_manifest_names_the_file_line_column_and_accepted_values() {
    let directory = TestDirectory::new("manifest");
    let manifest = directory.path().join("manifest.tsv");
    write(&manifest, "sample_id\tcondition\ns1\tpatient\n");
    let junctions = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/linear.tsv");
    let output = run(&[
        "gene",
        "GENE_LINEAR",
        &junctions.display().to_string(),
        &manifest.display().to_string(),
    ]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&manifest.display().to_string()));
    assert!(stderr.contains("line 2, column 'condition'"));
    assert!(stderr.contains("condition must be 'case' or 'control'"));
    assert!(stderr.contains("replace it with one of the accepted labels"));
}

#[test]
fn missing_sample_metadata_is_reported_without_silent_dropping() {
    let directory = TestDirectory::new("metadata");
    let junctions = directory.path().join("junctions.tsv");
    let manifest = directory.path().join("manifest.tsv");
    write(
        &junctions,
        "sample_id\tgene_id\tfrom\tto\tcount\ndata_only\tGENE_A\tE1\tE2\t10\n",
    );
    write(&manifest, "sample_id\tcondition\nmanifest_only\tcase\n");
    let output = run(&[
        "gene",
        "GENE_A",
        &junctions.display().to_string(),
        &manifest.display().to_string(),
    ]);
    let stdout = text(&output.stdout);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stdout.contains("no manifest entry for this sample_id"));
    assert!(stdout.contains("no junction data for gene 'GENE_A'"));
    assert!(stderr.contains("junction sample_id values match the manifest"));
}

#[test]
fn zero_support_names_the_sample_gene_and_correction() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/zero_support.tsv");
    let output = run(&["analyze", &fixture.display().to_string()]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&fixture.display().to_string()));
    assert!(stderr.contains("sample 'sample_zero', gene 'GENE_ZERO'"));
    assert!(stderr.contains("total read support is zero"));
    assert!(stderr.contains("junction count must be greater than zero"));
}

#[test]
fn comparison_requires_at_least_one_case_and_one_control() {
    let directory = TestDirectory::new("groups");
    let junctions = directory.path().join("junctions.tsv");
    let manifest = directory.path().join("manifest.tsv");
    write(
        &junctions,
        "sample_id\tgene_id\tfrom\tto\tcount\ncase_1\tGENE_A\tE1\tE2\t10\n",
    );
    write(&manifest, "sample_id\tcondition\ncase_1\tcase\n");
    let output = run(&[
        "gene",
        "GENE_A",
        &junctions.display().to_string(),
        &manifest.display().to_string(),
    ]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains("no eligible control samples remain"));
    assert!(stderr.contains("review the exclusions above"));
}

#[test]
fn existing_output_folder_is_not_overwritten() {
    let directory = TestDirectory::new("output");
    let output_directory = directory.path().join("results");
    fs::create_dir(&output_directory).expect("existing output directory should be created");
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/linear.tsv");
    let output = run(&[
        "analyze",
        &fixture.display().to_string(),
        "--output",
        &output_directory.display().to_string(),
    ]);
    let stderr = text(&output.stderr);

    assert!(!output.status.success());
    assert!(stderr.contains(&output_directory.display().to_string()));
    assert!(stderr.contains("already exists"));
    assert!(stderr.contains("choose a new folder"));
}
