use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

struct Fixture(PathBuf);

impl Fixture {
    fn new() -> Self {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root =
            std::env::temp_dir().join(format!("splice-girl-stats-{}-{stamp}", std::process::id()));
        fs::create_dir(&root).unwrap();
        let mut junctions = "sample_id\tgene_id\tfrom\tto\tcount\n".to_owned();
        let mut manifest = "sample_id\tcondition\tsubject_id\n".to_owned();
        for index in 0..8 {
            let sample = format!("S{index}");
            let condition = if index < 4 { "case" } else { "control" };
            manifest.push_str(&format!("{sample}\t{condition}\tperson_{index}\n"));
            junctions.push_str(&format!("{sample}\tHIGH\tE1\tE2\t10\n{sample}\tHIGH\tE2\tE3\t10\n{sample}\tFLAT\tE1\tE2\t10\n"));
            if index < 4 {
                junctions.push_str(&format!("{sample}\tHIGH\tE3\tE1\t10\n"));
            }
            if index == 0 || index == 4 {
                junctions.push_str(&format!("{sample}\tSMALL\tE1\tE2\t10\n"));
            }
        }
        // An unselected zero-support graph must not be analyzed or enter FDR.
        junctions.push_str("S0\tUNSELECTED\tE1\tE2\t0\n");
        fs::write(root.join("junctions.tsv"), junctions).unwrap();
        fs::write(root.join("manifest.tsv"), manifest).unwrap();
        fs::write(
            root.join("genes.tsv"),
            "gene_symbol\tgene_id\nx\tHIGH\ny\tFLAT\nz\tSMALL\nu\tMISSING\n",
        )
        .unwrap();
        Self(root)
    }

    fn arguments(&self, output: &str, alternative: &str, minimum: &str) -> Vec<String> {
        vec![
            "stats".to_owned(),
            self.0.join("junctions.tsv").display().to_string(),
            self.0.join("manifest.tsv").display().to_string(),
            "--genes".to_owned(),
            self.0.join("genes.tsv").display().to_string(),
            "--family".to_owned(),
            "toy<&>".to_owned(),
            "--alternative".to_owned(),
            alternative.to_owned(),
            "--design".to_owned(),
            "independent".to_owned(),
            "--min-per-group".to_owned(),
            minimum.to_owned(),
            "--output".to_owned(),
            self.0.join(output).display().to_string(),
        ]
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn run(arguments: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_splice-girl"))
        .args(arguments)
        .output()
        .unwrap()
}
fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

fn rows(path: &Path) -> Vec<BTreeMap<String, String>> {
    let contents = fs::read_to_string(path).unwrap();
    let mut lines = contents.lines();
    let header: Vec<_> = lines.next().unwrap().split('\t').collect();
    lines
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            assert_eq!(fields.len(), header.len());
            header
                .iter()
                .zip(fields)
                .map(|(key, value)| (key.to_string(), value.to_owned()))
                .collect()
        })
        .collect()
}

#[test]
fn stats_runs_exact_tests_bh_exports_reports_and_preserves_missing_genes() {
    let fixture = Fixture::new();
    let args = fixture.arguments("output", "two-sided", "2");
    let output = run(&args);
    assert!(output.status.success(), "{}", text(&output.stderr));
    let folder = fixture.0.join("output");
    let statistics = rows(&folder.join("statistics.tsv"));
    assert_eq!(
        statistics
            .iter()
            .map(|row| row["gene_id"].as_str())
            .collect::<Vec<_>>(),
        ["FLAT", "HIGH", "MISSING", "SMALL"]
    );
    let high = &statistics[1];
    assert_eq!(high["status"], "tested");
    assert_eq!(high["case_n"], "4");
    assert_eq!(high["control_n"], "4");
    assert_eq!(high["permutations_considered"], "70");
    assert_eq!(high["extreme_count"], "2");
    assert!(
        (high["mean_difference_case_minus_control"]
            .parse::<f64>()
            .unwrap()
            - 1.0)
            .abs()
            < 1e-9
    );
    assert!((high["p_value"].parse::<f64>().unwrap() - 2.0 / 70.0).abs() < 1e-9);
    assert!((high["adjusted_p_value"].parse::<f64>().unwrap() - 4.0 / 70.0).abs() < 1e-9);
    assert_eq!(high["tested_family_size"], "2");
    assert_eq!(statistics[0]["p_value"], "1.0000000000");
    for row in &statistics[2..] {
        assert_eq!(row["status"], "untested");
        assert!(row["p_value"].is_empty() && row["adjusted_p_value"].is_empty());
        assert!(!row["untested_reason"].is_empty());
    }
    let decisions = rows(&folder.join("comparisons.tsv"));
    assert_eq!(decisions.len(), 32);
    assert!(
        decisions
            .iter()
            .any(|row| row["gene_id"] == "MISSING" && row["status"] == "excluded")
    );
    let html = fs::read_to_string(folder.join("report.html")).unwrap();
    assert!(html.contains("Statistical tests"));
    assert!(html.contains("toy&lt;&amp;&gt;"));
    assert!(html.contains("0.0571428571"));
    assert!(html.contains("2 / 70"));
    assert!(html.contains("href=\"statistics.tsv\""));
    assert!(!html.contains("No p-values or adjusted p-values were calculated"));
    let manifest = fs::read_to_string(folder.join("run_manifest.json")).unwrap();
    assert!(manifest.contains("\"mode\": \"stats\""));
    assert!(manifest.contains("\"role\": \"gene_family\""));
    assert!(manifest.contains("\"tested_genes\": 2"));
    assert!(manifest.contains("\"alternative\": \"two-sided\""));
    assert!(
        !fs::read_to_string(folder.join("sample_metrics.tsv"))
            .unwrap()
            .contains("UNSELECTED")
    );
    assert_eq!(fs::read_dir(folder.join("graphs")).unwrap().count(), 18);
    let repeated = run(&args);
    assert!(!repeated.status.success());
    assert!(text(&repeated.stderr).contains("already exists"));
    assert_eq!(
        html,
        fs::read_to_string(folder.join("report.html")).unwrap()
    );
}

#[test]
fn sidedness_is_applied_instead_of_silently_using_two_sided() {
    let fixture = Fixture::new();
    for (direction, expected) in [("case-greater", 1.0 / 70.0), ("case-less", 1.0)] {
        let output = run(&fixture.arguments(direction, direction, "2"));
        assert!(output.status.success(), "{}", text(&output.stderr));
        let statistics = rows(&fixture.0.join(direction).join("statistics.tsv"));
        assert!((statistics[1]["p_value"].parse::<f64>().unwrap() - expected).abs() < 1e-9);
        assert_eq!(statistics[1]["alternative"], direction);
    }
}

#[test]
fn all_untestable_genes_save_diagnostics_but_exit_unsuccessfully() {
    let fixture = Fixture::new();
    let output = run(&fixture.arguments("output", "two-sided", "5"));
    assert!(!output.status.success());
    assert!(text(&output.stderr).contains("no requested gene was testable"));
    let statistics = rows(&fixture.0.join("output/statistics.tsv"));
    assert!(
        statistics
            .iter()
            .all(|row| row["p_value"].is_empty() && row["adjusted_p_value"].is_empty())
    );
    assert!(fixture.0.join("output/report.html").is_file());
}

#[test]
fn every_scientific_option_is_required_and_invalid_flags_fail_before_writing() {
    let fixture = Fixture::new();
    let base = fixture.arguments("output", "two-sided", "2");
    for option in [
        "--genes",
        "--family",
        "--alternative",
        "--design",
        "--min-per-group",
        "--output",
    ] {
        let mut args = base.clone();
        let index = args.iter().position(|value| value == option).unwrap();
        args.drain(index..index + 2);
        let output = run(&args);
        assert!(!output.status.success());
        assert!(text(&output.stderr).contains(&format!("stats requires {option}")));
    }
    for (option, value, message) in [
        ("--alternative", "either", "--alternative must"),
        ("--design", "paired", "supports only --design independent"),
        ("--min-per-group", "1", "at least 2"),
        ("--min-per-group", "abc", "integer"),
    ] {
        let mut args = base.clone();
        let index = args.iter().position(|entry| entry == option).unwrap();
        args[index + 1] = value.to_owned();
        let output = run(&args);
        assert!(!output.status.success());
        assert!(text(&output.stderr).contains(message));
    }
    let mut args = base;
    args.extend(["--family".to_owned(), "duplicate".to_owned()]);
    assert!(text(&run(&args).stderr).contains("duplicate stats option"));
    assert!(!fixture.0.join("output").exists());
}

#[test]
fn repeated_subjects_cannot_enter_an_independent_test() {
    let fixture = Fixture::new();
    let path = fixture.0.join("manifest.tsv");
    let manifest = fs::read_to_string(&path)
        .unwrap()
        .replace("person_1", "person_0");
    fs::write(path, manifest).unwrap();
    let output = run(&fixture.arguments("output", "two-sided", "2"));
    assert!(!output.status.success());
    assert!(text(&output.stderr).contains("appears more than once"));
    assert!(!fixture.0.join("output").exists());
}

#[test]
fn gene_family_input_errors_are_actionable_and_do_not_create_outputs() {
    let fixture = Fixture::new();
    for invalid in ["", "wrong_header\nHIGH\n", "gene_id\nHIGH\nHIGH\n"] {
        fs::write(fixture.0.join("genes.tsv"), invalid).unwrap();
        let output = run(&fixture.arguments("output", "two-sided", "2"));
        assert!(!output.status.success());
        assert!(text(&output.stderr).contains("invalid gene family"));
        assert!(!fixture.0.join("output").exists());
    }
}

#[test]
fn exact_enumeration_limit_is_reported_as_untested_without_a_fake_p_value() {
    let fixture = Fixture::new();
    let mut manifest = "sample_id\tcondition\n".to_owned();
    let mut junctions = "sample_id\tgene_id\tfrom\tto\tcount\n".to_owned();
    for index in 0..24 {
        let condition = if index < 12 { "case" } else { "control" };
        manifest.push_str(&format!("S{index}\t{condition}\n"));
        junctions.push_str(&format!("S{index}\tLARGE\tE1\tE2\t1\n"));
    }
    fs::write(fixture.0.join("manifest.tsv"), manifest).unwrap();
    fs::write(fixture.0.join("junctions.tsv"), junctions).unwrap();
    fs::write(fixture.0.join("genes.tsv"), "gene_id\nLARGE\n").unwrap();
    let output = run(&fixture.arguments("output", "two-sided", "2"));
    assert!(!output.status.success());
    let statistics = rows(&fixture.0.join("output/statistics.tsv"));
    assert_eq!(statistics[0]["status"], "untested");
    assert!(statistics[0]["p_value"].is_empty());
    assert!(statistics[0]["untested_reason"].contains("enumeration limit"));
}

#[test]
fn menu_option_seven_runs_the_same_statistics_workflow() {
    let fixture = Fixture::new();
    let mut child = Command::new(env!("CARGO_BIN_EXE_splice-girl"))
        .arg("--start")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let input = format!(
        "7\n{}\n{}\n{}\ntoy\ntwo-sided\nindependent\n2\n{}\n\nq\n",
        fixture.0.join("junctions.tsv").display(),
        fixture.0.join("manifest.tsv").display(),
        fixture.0.join("genes.tsv").display(),
        fixture.0.join("output").display()
    );
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    assert!(!text(&output.stdout).contains("  error:"));
    assert!(fixture.0.join("output/statistics.tsv").is_file());
}

#[test]
fn baselines_use_identical_samples_and_leave_missing_values_blank() {
    let fixture = Fixture::new();
    let output = run(&fixture.arguments("baselines", "two-sided", "2"));
    assert!(output.status.success(), "{}", text(&output.stderr));
    let folder = fixture.0.join("baselines");
    let baseline = rows(&folder.join("baseline_comparisons.tsv"));
    let support = baseline
        .iter()
        .find(|r| r["gene_id"] == "HIGH" && r["metric"] == "total_junction_support")
        .unwrap();
    assert_eq!(support["case_mean"], "30.0000000000");
    assert_eq!(support["control_mean"], "20.0000000000");
    assert_eq!(
        support["mean_difference_case_minus_control"],
        "10.0000000000"
    );
    assert_eq!(support["observed_topology"], "different_review_required");
    let cycle = baseline
        .iter()
        .find(|r| r["gene_id"] == "HIGH" && r["metric"] == "cycle_rank")
        .unwrap();
    assert_eq!(cycle["case_mean"], "1.0000000000");
    assert_eq!(cycle["control_mean"], "0.0000000000");
    assert!(
        baseline.iter().filter(|r| r["gene_id"] == "MISSING").all(
            |r| r["case_mean"].is_empty() && r["mean_difference_case_minus_control"].is_empty()
        )
    );
    assert_eq!(rows(&folder.join("baseline_metrics.tsv")).len(), 18);
    assert_eq!(
        fs::read(folder.join("gene_family.tsv")).unwrap(),
        fs::read(fixture.0.join("genes.tsv")).unwrap()
    );
    assert!(!folder.join("rankings.tsv").exists());
}

#[test]
fn requested_ranking_is_deterministic_and_does_not_rank_untested_genes() {
    let fixture = Fixture::new();
    for rule in ["absolute-effect", "adjusted-p"] {
        let mut args = fixture.arguments(rule, "two-sided", "2");
        args.extend(["--rank-by".to_owned(), rule.to_owned()]);
        let output = run(&args);
        assert!(output.status.success(), "{}", text(&output.stderr));
        let folder = fixture.0.join(rule);
        let rankings = rows(&folder.join("rankings.tsv"));
        assert_eq!(
            rankings
                .iter()
                .map(|r| r["gene_id"].as_str())
                .collect::<Vec<_>>(),
            ["HIGH", "FLAT", "MISSING", "SMALL"]
        );
        assert_eq!(rankings[0]["rank"], "1");
        assert_eq!(rankings[1]["rank"], "2");
        assert!(rankings[2..].iter().all(|r| r["rank"].is_empty()));
        let manifest = fs::read_to_string(folder.join("run_manifest.json")).unwrap();
        assert!(manifest.contains(&format!("\"ranking_rule\": \"{rule}\"")));
    }
    // Equal effect and p-value: lexical gene ID decides, regardless of input order.
    let junctions = fs::read_to_string(fixture.0.join("junctions.tsv")).unwrap();
    let mut tied = junctions.clone();
    for line in junctions
        .lines()
        .skip(1)
        .filter(|line| line.contains("\tHIGH\t"))
    {
        tied.push_str(&line.replace("\tHIGH\t", "\tALPHA\t"));
        tied.push('\n');
    }
    fs::write(fixture.0.join("junctions.tsv"), tied).unwrap();
    fs::write(fixture.0.join("genes.tsv"), "gene_id\nHIGH\nALPHA\n").unwrap();
    let mut args = fixture.arguments("ties", "two-sided", "2");
    args.extend(["--rank-by".to_owned(), "adjusted-p".to_owned()]);
    assert!(run(&args).status.success());
    assert_eq!(
        rows(&fixture.0.join("ties/rankings.tsv"))[0]["gene_id"],
        "ALPHA"
    );
    let mut args = fixture.arguments("invalid-rank", "two-sided", "2");
    args.extend(["--rank-by".to_owned(), "whatever".to_owned()]);
    assert!(text(&run(&args).stderr).contains("--rank-by must"));
    assert!(!fixture.0.join("invalid-rank").exists());
}
