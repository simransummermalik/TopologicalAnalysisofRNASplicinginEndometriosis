use super::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static NUMBER: AtomicUsize = AtomicUsize::new(0);
struct Fixture {
    root: PathBuf,
}

impl Fixture {
    fn new() -> Self {
        let time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let number = NUMBER.fetch_add(1, Ordering::Relaxed);
        let root = std::env::temp_dir().join(format!(
            "splice-girl-config-{}-{time}-{number}",
            std::process::id()
        ));
        fs::create_dir(&root).unwrap();
        for (name, text) in [
            (
                "junctions.tsv",
                include_str!("../../../examples/shared_graph/junctions.tsv"),
            ),
            (
                "samples.tsv",
                include_str!("../../../examples/shared_graph/samples.tsv"),
            ),
            (
                "genes.tsv",
                include_str!("../../../examples/shared_graph/genes.tsv"),
            ),
            (
                "analysis.json",
                include_str!("../../../examples/shared_graph/analysis.json"),
            ),
        ] {
            fs::write(root.join(name), text).unwrap();
        }
        Self { root }
    }
    fn config(&self) -> PathBuf {
        self.root.join("analysis.json")
    }
    fn change(&self, file: &str, update: impl FnOnce(&mut serde_json::Value)) {
        let path = self.root.join(file);
        let mut value: serde_json::Value = read_json(&path).unwrap();
        update(&mut value);
        fs::write(path, serde_json::to_string_pretty(&value).unwrap()).unwrap();
    }
    fn plan(&self) -> Result<Plan, String> {
        Configuration::load(&self.config())
    }
    fn fails(&self, fragment: &str) {
        let error = match self.plan() {
            Ok(_) => panic!("expected preflight failure"),
            Err(error) => error,
        };
        assert!(error.contains(fragment), "{error}");
        assert!(!self.root.join("configured-results").exists());
    }
    fn pipeline(&self) {
        let index: Vec<_> = (1..=8).map(|part| json!({"path":format!("/reference/genome.{part}.ht2"),"sha256":format!("{part:064x}")})).collect();
        let mut reference_outputs = index.clone();
        reference_outputs
            .push(json!({"path":"/reference/splice_sites.tsv","sha256":"d".repeat(64)}));
        fs::write(self.root.join("reference.json"), json!({"status":"complete","method":"HISAT2 reference preparation","inputs":[{"path":"genome.fa.gz","sha256":"a".repeat(64)},{"path":"genes.gtf.gz","sha256":"b".repeat(64)}],"outputs":reference_outputs}).to_string()).unwrap();
        fs::write(self.root.join("map.json"), json!({"schema_version":1,"status":"complete","node_rule":"shared-splice-site","input":{"path":"genes.gtf.gz","sha256":"b".repeat(64)},"output":{"path":"/reference/map.tsv","sha256":"e".repeat(64)}}).to_string()).unwrap();
        fs::create_dir(self.root.join("batch")).unwrap();
        let raw = fs::read_to_string(self.root.join("junctions.tsv")).unwrap();
        let mut samples = Vec::new();
        for number in 0..8 {
            let id = format!("S{number}");
            let folder = self.root.join("batch").join(&id);
            fs::create_dir(&folder).unwrap();
            let rows = raw
                .lines()
                .skip(1)
                .filter(|row| row.starts_with(&format!("{id}\t")))
                .collect::<Vec<_>>()
                .join("\n");
            fs::write(
                folder.join("junctions.tsv"),
                format!("{}\n{rows}\n", raw.lines().next().unwrap()),
            )
            .unwrap();
            let mut inputs = index.clone();
            inputs.extend([
                json!({"path":"/reference/map.tsv","sha256":"e".repeat(64)}),
                json!({"path":"/reference/splice_sites.tsv","sha256":"d".repeat(64)}),
            ]);
            let run = json!({"pipeline":"Splice Girl Cleaver/HISAT2","status":"complete","settings":{"sample_id":id,"strandedness":"RF","minimum_unique_reads":1,"minimum_anchor":8,"exclude_duplicates":false,"count_unit":"primary uniquely mapped read alignments with NH:i:1 and XS strand"},"inputs":inputs,"outputs":[{"path":"junctions.tsv","sha256":sha256_file(&folder.join("junctions.tsv")).unwrap()}],"commands":[["hisat2","--known-splicesite-infile","/reference/splice_sites.tsv"],["python","sam_junctions.py","--gene-map","/reference/map.tsv"]],"stages":[{"status":"complete"}]});
            let mut run = run;
            run["commands"] = json!([
                [
                    "hisat2",
                    "-x",
                    "/reference/genome",
                    "--rna-strandness",
                    "RF",
                    "--known-splicesite-infile",
                    "/reference/splice_sites.tsv"
                ],
                [
                    "python",
                    "sam_junctions.py",
                    "--gene-map",
                    "/reference/map.tsv",
                    "--min-anchor",
                    "8",
                    "--min-unique",
                    "1"
                ]
            ]);
            fs::write(folder.join("run_manifest.json"), run.to_string()).unwrap();
            samples.push(json!({"sample_id":id,"condition":if number < 4 {"case"} else {"control"},"status":"complete","run_manifest":{"path":format!("S{number}/run_manifest.json"),"sha256":sha256_file(&folder.join("run_manifest.json")).unwrap()}}));
        }
        fs::write(self.root.join("batch/batch_manifest.json"), json!({"pipeline":"Splice Girl sequential Cleaver/HISAT2 batch","status":"complete","samples":samples,"outputs":[{"path":"junctions.tsv","sha256":sha256_file(&self.root.join("junctions.tsv")).unwrap()}]}).to_string()).unwrap();
        self.change("analysis.json", |cfg| {
            cfg["input_kind"] = json!("pipeline");
            cfg["provenance"] = json!({"batch_manifest":"batch/batch_manifest.json","reference_manifest":"reference.json","gene_map_manifest":"map.json","node_rule":"shared-splice-site","strandedness":"RF","minimum_anchor":8,"exclude_duplicates":false});
        });
    }
    fn refresh_child(&self) {
        let digest = sha256_file(&self.root.join("batch/S0/run_manifest.json")).unwrap();
        self.change("batch/batch_manifest.json", |batch| {
            batch["samples"][0]["run_manifest"]["sha256"] = json!(digest)
        });
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

#[test]
fn synthetic_preflight_is_read_only_and_paths_follow_the_configuration() {
    let fixture = Fixture::new();
    let plan = fixture.plan().unwrap();
    assert_eq!(plan.report["cohort_samples"], 8);
    assert_eq!(plan.report["primary_eligible_genes"], 1);
    assert!(plan.inputs.junctions.is_absolute());
    assert!(!fixture.root.join("configured-results").exists());
    plan.execute(&["splice-girl".to_owned(), "run".to_owned()])
        .unwrap();
    let output = fixture.root.join("configured-results");
    let manifest: serde_json::Value = read_json(&output.join("run_manifest.json")).unwrap();
    assert_eq!(manifest["status"], "complete");
    assert!(output.join("report.html").is_file());
    assert!(
        fs::read_to_string(output.join("report.html"))
            .unwrap()
            .contains("Synthetic demonstration")
    );
    assert_eq!(
        fs::read_to_string(output.join("configuration.json")).unwrap(),
        fs::read_to_string(fixture.config()).unwrap()
    );
}

#[test]
fn configuration_rejects_unknown_duplicate_and_mistyped_fields() {
    let fixture = Fixture::new();
    fixture.change("analysis.json", |cfg| {
        cfg["shared_graph"]["minimum_edeg_count"] = json!(5)
    });
    fixture.fails("unknown field");
    fs::write(
        fixture.config(),
        "{\"schema_version\":1,\"schema_version\":2}",
    )
    .unwrap();
    fixture.fails("duplicate field");
    fs::write(
        fixture.config(),
        include_str!("../../../examples/shared_graph/analysis.json"),
    )
    .unwrap();
    fixture.change("analysis.json", |cfg| {
        cfg["shared_graph"]["minimum_gene_support"] = json!("50")
    });
    fixture.fails("invalid type");
}

#[test]
fn configuration_requires_explicit_complete_counts_and_valid_statistics() {
    let fixture = Fixture::new();
    fixture.change("analysis.json", |cfg| cfg["complete_counts"] = json!(false));
    fixture.fails("--complete-counts");
    fixture.change("analysis.json", |cfg| {
        cfg["complete_counts"] = json!(true);
        cfg["statistics"]["design"] = json!("paired");
    });
    fixture.fails("independent");
    fixture.change("analysis.json", |cfg| {
        cfg["statistics"]["design"] = json!("independent");
        cfg["statistics"]["min_per_group"] = json!(5);
    });
    fixture.fails("min_per_group");
}

#[test]
fn existing_output_and_changed_checked_inputs_are_protected() {
    let fixture = Fixture::new();
    let plan = fixture.plan().unwrap();
    fs::write(fixture.root.join("genes.tsv"), "gene_id\nOTHER\n").unwrap();
    assert!(
        plan.execute(&[])
            .unwrap_err()
            .contains("changed during analysis")
    );
    assert!(!fixture.root.join("configured-results").exists());
    fs::create_dir(fixture.root.join("configured-results")).unwrap();
    fs::write(fixture.root.join("configured-results/sentinel"), "keep").unwrap();
    let error = match fixture.plan() {
        Ok(_) => panic!("expected existing-output failure"),
        Err(error) => error,
    };
    assert!(error.contains("already exists"));
    assert_eq!(
        fs::read_to_string(fixture.root.join("configured-results/sentinel")).unwrap(),
        "keep"
    );
}

#[test]
fn complete_pipeline_checks_all_sample_counts_and_recorded_identities() {
    let fixture = Fixture::new();
    fixture.pipeline();
    let plan = fixture.plan().unwrap();
    assert_eq!(plan.report["input_kind"], "pipeline");
    assert!(plan.checked.len() >= 23);
    plan.execute(&[]).unwrap();
    assert!(
        fixture
            .root
            .join("configured-results/report.html")
            .is_file()
    );
    let resolved = fixture
        .root
        .join("configured-results/resolved_configuration.json");
    let mut saved: serde_json::Value = read_json(&resolved).unwrap();
    assert_eq!(saved["schema_version"], 1);
    saved["output"] = json!(fixture.root.join("rerun"));
    let rerun = fixture.root.join("rerun.json");
    fs::write(&rerun, saved.to_string()).unwrap();
    assert!(Configuration::load(&rerun).is_ok());
}

#[test]
fn incomplete_or_tampered_sample_provenance_fails_before_outputs() {
    let fixture = Fixture::new();
    fixture.pipeline();
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["status"] = json!("failed")
    });
    fixture.fails("checksum mismatch");
    fixture.refresh_child();
    fixture.fails("processing is incomplete");
}

#[test]
fn assembled_runs_are_accepted_and_unknown_pipelines_refused() {
    let fixture = Fixture::new();
    fixture.pipeline();
    fixture.change("batch/batch_manifest.json", |batch| {
        batch["pipeline"] = json!("Splice Girl assembled Cleaver/HISAT2 runs")
    });
    assert!(fixture.plan().is_ok());
    fixture.change("batch/batch_manifest.json", |batch| {
        batch["pipeline"] = json!("hand-merged tables")
    });
    fixture.fails("unsupported pipeline");
}

#[test]
fn sample_held_by_the_strand_check_is_refused() {
    let fixture = Fixture::new();
    fixture.pipeline();
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["strand_check"] = json!({"status": "pass", "read1_agreement": 0.99})
    });
    fixture.refresh_child();
    assert!(fixture.plan().is_ok());
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["strand_check"] = json!({"status": "hold_for_review", "read1_agreement": 0.9})
    });
    fixture.refresh_child();
    fixture.fails("failed the strand check");
}

#[test]
fn mismatched_reference_and_node_rules_are_rejected() {
    let fixture = Fixture::new();
    fixture.pipeline();
    fixture.change("map.json", |map| map["node_rule"] = json!("exon-interval"));
    fixture.fails("node rule");
    fixture.change("map.json", |map| {
        map["node_rule"] = json!("shared-splice-site");
        map["input"]["sha256"] = json!("f".repeat(64));
    });
    fixture.fails("GTF differs");
    fixture.change("map.json", |map| {
        map["input"]["sha256"] = json!("b".repeat(64))
    });
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["inputs"][0]["sha256"] = json!("f".repeat(64))
    });
    fixture.refresh_child();
    fixture.fails("reference index or splice hints differ");
}

#[test]
fn anchor_and_prefiltered_support_mismatches_are_rejected() {
    let fixture = Fixture::new();
    fixture.pipeline();
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["settings"]["minimum_anchor"] = json!(10)
    });
    fixture.refresh_child();
    fixture.fails("preprocessing settings differ");
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["settings"]["minimum_anchor"] = json!(8);
        run["settings"]["minimum_unique_reads"] = json!(3);
    });
    fixture.refresh_child();
    fixture.fails("support must be 1");
}

#[test]
fn merged_counts_cannot_omit_an_exported_edge_even_after_hashes_are_updated() {
    let fixture = Fixture::new();
    fixture.pipeline();
    let child = fixture.root.join("batch/S0/junctions.tsv");
    let text = fs::read_to_string(&child).unwrap();
    fs::write(
        &child,
        text.lines()
            .filter(|line| !line.contains("\tTREE\t"))
            .collect::<Vec<_>>()
            .join("\n")
            + "\n",
    )
    .unwrap();
    let digest = sha256_file(&child).unwrap();
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["outputs"][0]["sha256"] = json!(digest)
    });
    fixture.refresh_child();
    fixture.fails("combined and per-sample exports");
}

#[test]
fn all_excluded_primary_still_produces_a_honest_report() {
    let fixture = Fixture::new();
    fixture.change("analysis.json", |cfg| {
        cfg["shared_graph"]["minimum_edge_count"] = json!(1000000);
        cfg["sensitivity"] = json!(null);
    });
    let plan = fixture.plan().unwrap();
    assert_eq!(plan.report["primary_eligible_genes"], 0);
    plan.execute(&[]).unwrap();
    assert!(
        fixture
            .root
            .join("configured-results/report.html")
            .is_file()
    );
    assert!(
        !fixture
            .root
            .join("configured-results/analysis/primary/statistics")
            .exists()
    );
}

#[test]
fn recorded_commands_must_agree_with_claimed_processing_settings() {
    let fixture = Fixture::new();
    fixture.pipeline();
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["commands"][0][4] = json!("FR")
    });
    fixture.refresh_child();
    fixture.fails("contradicts strandedness");
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["commands"][0][4] = json!("RF");
        run["commands"][1][7] = json!("3");
    });
    fixture.refresh_child();
    fixture.fails("contradicts anchor/support");
    fixture.change("batch/S0/run_manifest.json", |run| {
        run["commands"][1][7] = json!("1");
        run["commands"][0][2] = json!("/other/genome");
    });
    fixture.refresh_child();
    fixture.fails("recorded input identity");
}
