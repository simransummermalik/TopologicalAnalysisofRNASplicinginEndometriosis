//! Argument handling for the installed shared-graph workflow.
use super::sensitivity::StatisticsSettings;
use super::*;

pub const PREPARE_HELP: &str = "Prepare the same graph for every covered sample of each gene.

usage: splice-girl prepare-shared --junctions <counts.tsv> --manifest <samples.tsv>
  [--genes <family.tsv>] --complete-counts --minimum-edge-count <n>
  --minimum-samples-per-edge <n> --minimum-gene-support <n> --output <new-folder>

All three thresholds are required positive integers. Edge selection ignores
condition labels. --complete-counts declares successful unfiltered count exports
for every manifest sample (read-counter minimum support 1). An absent edge
becomes zero only in a covered sample/gene. Genes need a cycle and coverage in
every sample to enter eligible_genes.tsv. Inspect all eligibility tables.
No Python, repository checkout or external tools are needed for this command.
Try: splice-girl example --output demo-inputs";

pub const SENSITIVITY_HELP: &str =
    "Rerun native shared-graph analysis across explicit one-factor changes.

usage: splice-girl sensitivity --junctions <counts.tsv> --manifest <samples.tsv>
  [--genes <family.tsv>] --complete-counts --primary-edge-count <n>
  --primary-sample-count <n> --primary-gene-support <n>
  --minimum-junction-counts <3,5,10> --minimum-samples-per-edge <2,4,6>
  --minimum-gene-support <20,50,100> --output <new-folder>

Optional statistics require all of:
  --run-stats --genes <family.tsv> --family <label>
  --alternative <two-sided|case-greater|case-less>
  --design independent --min-per-group <n>

Runs the explicit primary setting, then unique one-factor alternatives.
Saves shared inputs, exclusions, comparisons and optional statistics/reports.
Missing effects stay blank; all-excluded scenarios skip statistics. Anchor
sensitivity requires separate count exports and is not evaluated here.
The caller must predeclare the primary setting; no result is selected for them.
No Python, subprocesses or repository checkout are needed.";

const EXAMPLE_HELP: &str = "Write small invented eight-sample inputs for the installed CLI.

usage: splice-girl example --output <new-folder>

Writes junctions.tsv, samples.tsv, genes.tsv, analysis.json and a README.
No downloads or human data. Existing folders are protected.";

/// Dispatch three standalone commands, with errors before output creation.
pub fn run_command(name: &str, arguments: &[String], command: &[String]) -> Result<(), String> {
    let help = match name {
        "prepare-shared" => PREPARE_HELP,
        "sensitivity" => SENSITIVITY_HELP,
        "example" => EXAMPLE_HELP,
        _ => return Err("unknown shared-graph command".to_owned()),
    };
    if arguments == ["--help"] || arguments == ["-h"] {
        println!("{help}");
        return Ok(());
    }
    if arguments.is_empty() {
        return Err(help.to_owned());
    }
    let values: &[&str] = match name {
        "example" => &["--output", "--output-dir"],
        "prepare-shared" => &[
            "--junctions",
            "--manifest",
            "--genes",
            "--minimum-edge-count",
            "--minimum-samples-per-edge",
            "--minimum-gene-support",
            "--output",
            "--output-dir",
        ],
        _ => &[
            "--junctions",
            "--manifest",
            "--genes",
            "--primary-edge-count",
            "--primary-sample-count",
            "--primary-gene-support",
            "--minimum-junction-counts",
            "--minimum-samples-per-edge",
            "--minimum-gene-support",
            "--family",
            "--alternative",
            "--design",
            "--min-per-group",
            "--output",
            "--output-dir",
        ],
    };
    let flags: &[&str] = match name {
        "example" => &[],
        "prepare-shared" => &["--complete-counts"],
        _ => &["--complete-counts", "--run-stats"],
    };
    let mut options = BTreeMap::new();
    let mut arguments = arguments.iter();
    while let Some(argument) = arguments.next() {
        let canonical = if argument == "--output-dir" {
            "--output"
        } else {
            argument.as_str()
        };
        let value = if flags.contains(&argument.as_str()) {
            "true"
        } else if values.contains(&argument.as_str()) {
            arguments
                .next()
                .filter(|value| !value.starts_with("--") && !value.trim().is_empty())
                .map(String::as_str)
                .ok_or_else(|| format!("{argument} needs a value; use splice-girl {name} --help"))?
        } else {
            return Err(format!(
                "unknown {name} option {argument:?}; use splice-girl {name} --help"
            ));
        };
        if options.insert(canonical, value).is_some() {
            return Err(format!("duplicate option {canonical}"));
        }
    }
    let required = |flag: &str| {
        options
            .get(flag)
            .copied()
            .ok_or_else(|| format!("{name} requires {flag}; use splice-girl {name} --help"))
    };
    let output = Path::new(required("--output")?);
    if output.symlink_metadata().is_ok() {
        return Err(format!(
            "output path {} already exists; choose a new folder",
            output.display()
        ));
    }
    if name == "example" {
        write_example(output, command)?;
        println!(
            "Synthetic example: {}\nRead README.md inside that folder for the next commands.",
            output.display()
        );
        return Ok(());
    }
    let inputs = Inputs {
        junctions: required("--junctions")?.into(),
        manifest: required("--manifest")?.into(),
        genes: options.get("--genes").map(PathBuf::from),
        complete_counts: options.contains_key("--complete-counts"),
    };
    let number = |flag: &str| parse_number(required(flag)?, flag);
    let samples =
        |flag: &str| usize::try_from(number(flag)?).map_err(|_| format!("{flag} is too large"));
    if name == "prepare-shared" {
        let settings = Settings {
            minimum_edge_count: number("--minimum-edge-count")?,
            minimum_samples_per_edge: samples("--minimum-samples-per-edge")?,
            minimum_gene_support: number("--minimum-gene-support")?,
        };
        let eligible = prepare(&inputs, settings, output, command)?;
        println!(
            "Shared graphs: {}\nGenes eligible in every sample: {eligible}\nReview gene_eligibility.tsv and sample_eligibility.tsv.",
            output.display()
        );
    } else {
        let statistics = if options.contains_key("--run-stats") {
            if required("--design")? != "independent" {
                return Err("only --design independent is supported".to_owned());
            }
            required("--genes")?;
            Some(StatisticsSettings {
                family: required("--family")?.to_owned(),
                alternative: crate::statistics::parse_alternative(required("--alternative")?)?,
                min_per_group: samples("--min-per-group")?,
            })
        } else {
            if ["--family", "--alternative", "--design", "--min-per-group"]
                .iter()
                .any(|name| options.contains_key(name))
            {
                return Err("statistical options require --run-stats; omit them for descriptive sensitivity".to_owned());
            }
            None
        };
        let sample_grid = parse_grid(
            required("--minimum-samples-per-edge")?,
            "--minimum-samples-per-edge",
        )?
        .into_iter()
        .map(|n| usize::try_from(n).map_err(|_| "sample grid value is too large".to_owned()))
        .collect::<Result<_, _>>()?;
        let settings = SensitivityOptions {
            primary: Settings {
                minimum_edge_count: number("--primary-edge-count")?,
                minimum_samples_per_edge: samples("--primary-sample-count")?,
                minimum_gene_support: number("--primary-gene-support")?,
            },
            edge_counts: parse_grid(
                required("--minimum-junction-counts")?,
                "--minimum-junction-counts",
            )?,
            sample_counts: sample_grid,
            gene_support: parse_grid(
                required("--minimum-gene-support")?,
                "--minimum-gene-support",
            )?,
            statistics,
        };
        let count = run_sensitivity(&inputs, &settings, output, command)?;
        println!(
            "Shared-graph sensitivity: {}\nScenarios completed: {count}\nOpen sensitivity.tsv for eligibility and effects.",
            output.display()
        );
    }
    Ok(())
}

fn parse_number(text: &str, name: &str) -> Result<u64, String> {
    if text.is_empty() || !text.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(format!("{name} must contain positive integers"));
    }
    text.parse()
        .ok()
        .filter(|value| *value > 0)
        .ok_or_else(|| format!("{name} must contain positive integers"))
}

fn parse_grid(text: &str, name: &str) -> Result<Vec<u64>, String> {
    text.split(',')
        .map(|part| parse_number(part, name))
        .collect()
}

fn write_example(output: &Path, command: &[String]) -> Result<(), String> {
    // Compile these small text inputs into the binary: no installed data path.
    let files = [
        (
            "junctions.tsv",
            include_str!("../../examples/shared_graph/junctions.tsv"),
        ),
        (
            "samples.tsv",
            include_str!("../../examples/shared_graph/samples.tsv"),
        ),
        (
            "genes.tsv",
            include_str!("../../examples/shared_graph/genes.tsv"),
        ),
        (
            "README.md",
            include_str!("../../examples/shared_graph/INSTALLED.md"),
        ),
        (
            "analysis.json",
            include_str!("../../examples/shared_graph/analysis.json"),
        ),
    ];
    // The example has no biological completeness declaration or cohort inputs.
    if let Some(parent) = output.parent().filter(|path| !path.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    fs::create_dir(output).map_err(|e| format!("could not create {}: {e}", output.display()))?;
    let result: Result<Vec<String>, String> = (|| {
        let mut checksums = Vec::new();
        for (name, text) in files {
            let path = output.join(name);
            OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(&path)
                .and_then(|mut file| file.write_all(text.as_bytes()))
                .map_err(|e| e.to_string())?;
            checksums.push(file_record(&path, name)?);
        }
        Ok(checksums)
    })();
    let (status, error, outputs) = match &result {
        Ok(files) => ("complete", "null".to_owned(), files.join(",")),
        Err(error) => ("failed", json_string(error), String::new()),
    };
    fs::write(output.join("run_manifest.json"), format!("{{\"status\":{},\"error\":{},\"mode\":\"synthetic example\",\"version\":{},\"command\":[{}],\"outputs\":[{}]}}\n", json_string(status), error, json_string(env!("CARGO_PKG_VERSION")), command.iter().map(|value| json_string(value)).collect::<Vec<_>>().join(","), outputs)).map_err(|e| e.to_string())?;
    result.map(|_| ())
}
