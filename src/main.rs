//! CLI entry point for `splice-girl`. The direct command accepts a junction
//! TSV, `gene` compares one gene, `compare` summarizes every observed gene,
//! `validate` audits inputs, and `--start` opens an interactive menu. ROADMAP.md's
//! `gene` syntax is directory-based (`--case dir/ --control dir/`); this
//! implementation instead takes one junction TSV (already carrying
//! `sample_id`) plus a separate manifest TSV, matching the formats this
//! tool actually parses. `analyze` is also available as an explicit alias
//! for the direct command. Inferential statistics and ranking remain future
//! work until their rules are approved.
use std::env;
use std::io::{self, IsTerminal, Write};
use std::path::Path;
use std::process::ExitCode;

use splice_girl::analyze_junctions;
use splice_girl::audit::{ValidationSummary, summarize_validation};
use splice_girl::compare::{
    ComparisonStatus, GeneComparison, GeneComparisonSummary, compare_all_genes, compare_gene,
};
use splice_girl::manifest::parse_manifest;
use splice_girl::output::{
    ExportSummary, write_analysis_bundle, write_compare_bundle, write_gene_bundle,
};
use splice_girl::parser::parse_tsv;
use splice_girl::visualize::write_svg;

fn main() -> ExitCode {
    let command: Vec<String> = env::args().collect();
    print_banner();
    match run(&command) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn print_banner() {
    let banner = include_str!("../assets/splice_girl_banner.txt");
    if io::stdout().is_terminal() {
        print!("\x1b[38;5;217m{banner}\x1b[0m");
    } else {
        print!("{banner}");
    }
    println!();
    let _ = io::stdout().flush();
}

fn usage(program: &str) -> String {
    format!(
        "Splice Girl — turn RNA splice junctions into graphs and Hodge measurements\n\n\
usage:\n  \
  {program} --start\n  \
  {program} <junctions.tsv> [--output <folder>]\n  \
  {program} analyze <junctions.tsv> [--output <folder>]\n  \
  {program} validate <junctions.tsv> [manifest.tsv]\n  \
  {program} gene <GENE_ID> <junctions.tsv> <manifest.tsv> [--output <folder>]\n  \
  {program} compare <junctions.tsv> <manifest.tsv> [--output <folder>]\n  \
  {program} visualize <SAMPLE_ID> <GENE_ID> <junctions.tsv> <output.svg>\n\n\
commands:\n  \
  analyze     calculate graph and Hodge measurements\n  \
  validate    audit inputs, sample matching, graphs, and math checks\n  \
  gene        compare one gene across case/control samples\n  \
  compare     descriptively compare every observed gene across groups\n  \
  visualize   save one sample/gene graph as a browser-ready SVG\n  \
  --start     open the interactive menu"
    )
}

fn run(command: &[String]) -> Result<(), String> {
    let program = command
        .first()
        .cloned()
        .unwrap_or_else(|| "splice-girl".to_owned());
    let mut arguments = command.iter().skip(1).cloned();
    let first_argument = arguments.next().ok_or_else(|| usage(&program))?;

    if matches!(first_argument.as_str(), "--help" | "-h" | "help") {
        println!("{}", usage(&program));
        return Ok(());
    }

    if matches!(first_argument.as_str(), "--version" | "-V") {
        println!("splice-girl {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }

    if first_argument == "--start" {
        if arguments.next().is_some() {
            return Err(usage(&program));
        }
        return run_start_menu();
    }

    if first_argument == "gene" {
        let gene_id = arguments.next().ok_or_else(|| usage(&program))?;
        if gene_id.trim().is_empty() {
            return Err(
                "gene ID cannot be empty; provide the gene_id used in the junction TSV".to_owned(),
            );
        }
        let junctions_path = arguments.next().ok_or_else(|| usage(&program))?;
        let manifest_path = arguments.next().ok_or_else(|| usage(&program))?;
        let output_dir = parse_optional_output(arguments.collect(), &program)?;
        let (results, comparison) = run_gene_comparison(
            &gene_id,
            Path::new(&junctions_path),
            Path::new(&manifest_path),
        )?;
        if let Some(output_dir) = output_dir {
            let summary = write_gene_bundle(
                &output_dir,
                Path::new(&junctions_path),
                Path::new(&manifest_path),
                &results,
                &comparison,
                command,
            )?;
            print_export_summary(&summary);
        }
        return print_gene_comparison(&comparison);
    }

    if first_argument == "validate" {
        let junctions_path = arguments.next().ok_or_else(|| usage(&program))?;
        let remaining: Vec<String> = arguments.collect();
        let manifest_path = match remaining.as_slice() {
            [] => None,
            [path] => Some(path.as_str()),
            _ => return Err(usage(&program)),
        };
        let summary = validate_paths(Path::new(&junctions_path), manifest_path.map(Path::new))?;
        return print_validation_summary(&summary);
    }

    if first_argument == "compare" {
        let junctions_path = arguments.next().ok_or_else(|| usage(&program))?;
        let manifest_path = arguments.next().ok_or_else(|| usage(&program))?;
        let output_dir = parse_optional_output(arguments.collect(), &program)?;
        let (results, summaries) =
            run_all_gene_comparisons(Path::new(&junctions_path), Path::new(&manifest_path))?;
        if let Some(output_dir) = output_dir {
            let summary = write_compare_bundle(
                &output_dir,
                Path::new(&junctions_path),
                Path::new(&manifest_path),
                &results,
                &summaries,
                command,
            )?;
            print_export_summary(&summary);
        }
        return print_all_gene_comparisons(&summaries);
    }

    if first_argument == "visualize" {
        let sample_id = arguments.next().ok_or_else(|| usage(&program))?;
        let gene_id = arguments.next().ok_or_else(|| usage(&program))?;
        if sample_id.trim().is_empty() || gene_id.trim().is_empty() {
            return Err(
                "sample ID and gene ID cannot be empty; use values present in the junction TSV"
                    .to_owned(),
            );
        }
        let junctions_path = arguments.next().ok_or_else(|| usage(&program))?;
        let output_path = arguments.next().ok_or_else(|| usage(&program))?;
        if arguments.next().is_some() {
            return Err(usage(&program));
        }
        return create_visualization(
            &sample_id,
            &gene_id,
            Path::new(&junctions_path),
            Path::new(&output_path),
        );
    }

    if first_argument == "analyze" {
        let junctions_path = arguments.next().ok_or_else(|| usage(&program))?;
        let output_dir = parse_optional_output(arguments.collect(), &program)?;
        let results = analyze_path(Path::new(&junctions_path))?;
        if let Some(output_dir) = output_dir {
            let summary =
                write_analysis_bundle(&output_dir, Path::new(&junctions_path), &results, command)?;
            print_export_summary(&summary);
        }
        return print_results(&results);
    }

    let output_dir = parse_optional_output(arguments.collect(), &program)?;

    let results = analyze_path(Path::new(&first_argument))?;
    if let Some(output_dir) = output_dir {
        let summary =
            write_analysis_bundle(&output_dir, Path::new(&first_argument), &results, command)?;
        print_export_summary(&summary);
    }
    print_results(&results)
}

fn parse_optional_output(
    arguments: Vec<String>,
    program: &str,
) -> Result<Option<std::path::PathBuf>, String> {
    match arguments.as_slice() {
        [] => Ok(None),
        [flag, folder]
            if matches!(flag.as_str(), "--output" | "-o") && !folder.trim().is_empty() =>
        {
            Ok(Some(folder.into()))
        }
        _ => Err(usage(program)),
    }
}

fn analyze_path(path: &Path) -> Result<Vec<splice_girl::GroupResult>, String> {
    let junctions = parse_tsv(path)?;
    analyze_junctions(junctions).map_err(|error| {
        format!(
            "analysis of junction TSV {} failed: {error}",
            path.display()
        )
    })
}

fn run_gene_comparison(
    gene_id: &str,
    junctions_path: &Path,
    manifest_path: &Path,
) -> Result<(Vec<splice_girl::GroupResult>, GeneComparison), String> {
    let manifest = parse_manifest(manifest_path)?;
    let results = analyze_path(junctions_path)?;
    let comparison = compare_gene(&manifest, &results, gene_id);
    Ok((results, comparison))
}

fn run_all_gene_comparisons(
    junctions_path: &Path,
    manifest_path: &Path,
) -> Result<(Vec<splice_girl::GroupResult>, Vec<GeneComparisonSummary>), String> {
    let manifest = parse_manifest(manifest_path)?;
    let results = analyze_path(junctions_path)?;
    let summaries = compare_all_genes(&manifest, &results);
    Ok((results, summaries))
}

fn validate_paths(
    junctions_path: &Path,
    manifest_path: Option<&Path>,
) -> Result<ValidationSummary, String> {
    let junctions = parse_tsv(junctions_path)?;
    let results = analyze_junctions(junctions.clone()).map_err(|error| {
        format!(
            "analysis of junction TSV {} failed during validation: {error}",
            junctions_path.display()
        )
    })?;
    let manifest = manifest_path.map(parse_manifest).transpose()?;
    Ok(summarize_validation(
        &junctions,
        &results,
        manifest.as_deref(),
    ))
}

fn run_start_menu() -> Result<(), String> {
    loop {
        println!();
        println!("  Splice Girl");
        println!("  -----------");
        println!("  [1] Analyze a junction TSV");
        println!("  [2] Compare one gene across case/control samples");
        println!("  [3] Save one splice graph as an SVG");
        println!("  [4] Validate inputs and mathematical checks");
        println!("  [5] Compare every observed gene");
        println!("  [q] Quit");
        println!();

        let Some(choice) = prompt("  Choose an option: ")? else {
            return Ok(());
        };

        match choice.trim().to_ascii_lowercase().as_str() {
            "1" => {
                let Some(path) = prompt("  Junction TSV path: ")? else {
                    return Ok(());
                };
                let path = path.trim();
                if path.is_empty() {
                    println!("  Please enter a TSV path.");
                    continue;
                }

                let output_dir = prompt("  Output folder (leave blank for terminal only): ")?;
                let output_dir = output_dir
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(std::path::PathBuf::from);
                let menu_command = vec!["splice-girl".to_owned(), "--start".to_owned()];
                let outcome = analyze_path(Path::new(path)).and_then(|results| {
                    if let Some(output_dir) = &output_dir {
                        let summary = write_analysis_bundle(
                            output_dir,
                            Path::new(path),
                            &results,
                            &menu_command,
                        )?;
                        print_export_summary(&summary);
                    }
                    print_results(&results)
                });
                match outcome {
                    Ok(()) => {
                        if prompt("  Press Enter to return to the menu...")?.is_none() {
                            return Ok(());
                        }
                    }
                    Err(error) => println!("  error: {error}"),
                }
            }
            "2" => {
                let Some(gene_id) = prompt("  Gene ID: ")? else {
                    return Ok(());
                };
                let Some(junctions_path) = prompt("  Junction TSV path: ")? else {
                    return Ok(());
                };
                let Some(manifest_path) = prompt("  Sample manifest TSV path: ")? else {
                    return Ok(());
                };
                let gene_id = gene_id.trim();
                let junctions_path = junctions_path.trim();
                let manifest_path = manifest_path.trim();
                if gene_id.is_empty() || junctions_path.is_empty() || manifest_path.is_empty() {
                    println!("  Please enter a gene ID and both file paths.");
                    continue;
                }

                let output_dir = prompt("  Output folder (leave blank for terminal only): ")?;
                let output_dir = output_dir
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(std::path::PathBuf::from);
                let menu_command = vec!["splice-girl".to_owned(), "--start".to_owned()];
                let outcome = run_gene_comparison(
                    gene_id,
                    Path::new(junctions_path),
                    Path::new(manifest_path),
                )
                .and_then(|(results, comparison)| {
                    if let Some(output_dir) = &output_dir {
                        let summary = write_gene_bundle(
                            output_dir,
                            Path::new(junctions_path),
                            Path::new(manifest_path),
                            &results,
                            &comparison,
                            &menu_command,
                        )?;
                        print_export_summary(&summary);
                    }
                    print_gene_comparison(&comparison)
                });
                match outcome {
                    Ok(()) => {
                        if prompt("  Press Enter to return to the menu...")?.is_none() {
                            return Ok(());
                        }
                    }
                    Err(error) => println!("  error: {error}"),
                }
            }
            "3" => {
                let Some(sample_id) = prompt("  Sample ID: ")? else {
                    return Ok(());
                };
                let Some(gene_id) = prompt("  Gene ID: ")? else {
                    return Ok(());
                };
                let Some(junctions_path) = prompt("  Junction TSV path: ")? else {
                    return Ok(());
                };
                let Some(output_path) = prompt("  Output SVG path: ")? else {
                    return Ok(());
                };
                let sample_id = sample_id.trim();
                let gene_id = gene_id.trim();
                let junctions_path = junctions_path.trim();
                let output_path = output_path.trim();
                if sample_id.is_empty()
                    || gene_id.is_empty()
                    || junctions_path.is_empty()
                    || output_path.is_empty()
                {
                    println!("  Please enter a sample ID, gene ID, and both file paths.");
                    continue;
                }

                match create_visualization(
                    sample_id,
                    gene_id,
                    Path::new(junctions_path),
                    Path::new(output_path),
                ) {
                    Ok(()) => {
                        if prompt("  Press Enter to return to the menu...")?.is_none() {
                            return Ok(());
                        }
                    }
                    Err(error) => println!("  error: {error}"),
                }
            }
            "4" => {
                let Some(junctions_path) = prompt("  Junction TSV path: ")? else {
                    return Ok(());
                };
                let manifest_path = prompt("  Sample manifest path (optional): ")?;
                let junctions_path = junctions_path.trim();
                let manifest_path = manifest_path
                    .as_deref()
                    .map(str::trim)
                    .filter(|p| !p.is_empty());
                if junctions_path.is_empty() {
                    println!("  Please enter a junction TSV path.");
                    continue;
                }

                let outcome =
                    validate_paths(Path::new(junctions_path), manifest_path.map(Path::new))
                        .and_then(|summary| print_validation_summary(&summary));
                match outcome {
                    Ok(()) => {
                        if prompt("  Press Enter to return to the menu...")?.is_none() {
                            return Ok(());
                        }
                    }
                    Err(error) => println!("  error: {error}"),
                }
            }
            "5" => {
                let Some(junctions_path) = prompt("  Junction TSV path: ")? else {
                    return Ok(());
                };
                let Some(manifest_path) = prompt("  Sample manifest TSV path: ")? else {
                    return Ok(());
                };
                let junctions_path = junctions_path.trim();
                let manifest_path = manifest_path.trim();
                if junctions_path.is_empty() || manifest_path.is_empty() {
                    println!("  Please enter both file paths.");
                    continue;
                }

                let output_dir = prompt("  Output folder (leave blank for terminal only): ")?;
                let output_dir = output_dir
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(std::path::PathBuf::from);
                let menu_command = vec!["splice-girl".to_owned(), "--start".to_owned()];
                let outcome =
                    run_all_gene_comparisons(Path::new(junctions_path), Path::new(manifest_path))
                        .and_then(|(results, summaries)| {
                            if let Some(output_dir) = &output_dir {
                                let summary = write_compare_bundle(
                                    output_dir,
                                    Path::new(junctions_path),
                                    Path::new(manifest_path),
                                    &results,
                                    &summaries,
                                    &menu_command,
                                )?;
                                print_export_summary(&summary);
                            }
                            print_all_gene_comparisons(&summaries)
                        });
                match outcome {
                    Ok(()) => {
                        if prompt("  Press Enter to return to the menu...")?.is_none() {
                            return Ok(());
                        }
                    }
                    Err(error) => println!("  error: {error}"),
                }
            }
            "q" => return Ok(()),
            "" => {}
            _ => println!("  Please choose 1, 2, 3, 4, 5, or q."),
        }
    }
}

fn create_visualization(
    sample_id: &str,
    gene_id: &str,
    junctions_path: &Path,
    output_path: &Path,
) -> Result<(), String> {
    let results = analyze_path(junctions_path)?;
    let available: Vec<String> = results
        .iter()
        .map(|result| format!("{}/{}", result.sample_id, result.gene_id))
        .collect();
    let result = results
        .iter()
        .find(|result| result.sample_id == sample_id && result.gene_id == gene_id)
        .ok_or_else(|| {
            format!(
                "no graph found for sample '{sample_id}' and gene '{gene_id}'; available sample/gene pairs: {}",
                available.join(", ")
            )
        })?;

    if !result.decomposition.checks.passed {
        return Err(format!(
            "refusing to visualize sample '{sample_id}' and gene '{gene_id}' because its mathematical checks failed"
        ));
    }

    write_svg(output_path, result)?;
    println!("graph saved: {}", output_path.display());
    println!("open the SVG in a web browser to view it");
    Ok(())
}

fn prompt(message: &str) -> Result<Option<String>, String> {
    print!("{message}");
    io::stdout()
        .flush()
        .map_err(|error| format!("could not flush terminal prompt: {error}"))?;
    let mut line = String::new();
    let bytes_read = io::stdin()
        .read_line(&mut line)
        .map_err(|error| format!("could not read terminal input: {error}"))?;
    if bytes_read == 0 {
        Ok(None)
    } else {
        Ok(Some(line))
    }
}

fn print_results(results: &[splice_girl::GroupResult]) -> Result<(), String> {
    // One block of metrics per (sample, gene) group found in the input.
    for result in results {
        let decomposition = &result.decomposition;
        println!("sample: {}", result.sample_id);
        println!("gene: {}", result.gene_id);
        println!("nodes: {}", result.graph.nodes.len());
        println!("edges: {}", result.graph.edges.len());
        println!("signal_energy: {:.10}", decomposition.signal_energy);
        println!("gradient_energy: {:.10}", decomposition.gradient_energy);
        println!("cycle_energy: {:.10}", decomposition.cycle_energy);
        println!("cycle_fraction: {:.10}", decomposition.cycle_fraction);
        println!("components:");
        println!("  from\tto\tsignal\tgradient\tcycle");
        for (index, edge) in result.graph.edges.iter().enumerate() {
            println!(
                "  {}\t{}\t{:.10}\t{:.10}\t{:.10}",
                result.graph.nodes[edge.from],
                result.graph.nodes[edge.to],
                result.graph.signal[index],
                decomposition.gradient[index],
                decomposition.cycle[index]
            );
        }
        println!(
            "checks: {} (reconstruction={:.3e}, kernel={:.3e}, orthogonality={:.3e})",
            if decomposition.checks.passed {
                "PASS"
            } else {
                "FAIL"
            },
            decomposition.checks.reconstruction_error,
            decomposition.checks.kernel_error,
            decomposition.checks.orthogonality_error
        );
        println!();

        if !decomposition.checks.passed {
            return Err(format!(
                "mathematical checks failed for sample {} and gene {}",
                result.sample_id, result.gene_id
            ));
        }
    }

    Ok(())
}

fn print_gene_comparison(comparison: &GeneComparison) -> Result<(), String> {
    println!("gene: {}", comparison.gene_id);
    println!(
        "included: {} sample(s), excluded: {} sample(s)",
        comparison.included.len(),
        comparison.excluded.len()
    );

    if !comparison.included.is_empty() {
        println!();
        println!("included:");
        println!("  sample_id\tcondition\tcycle_fraction");
        for sample in &comparison.included {
            println!(
                "  {}\t{}\t{:.10}",
                sample.sample_id, sample.condition, sample.cycle_fraction
            );
        }
    }

    if !comparison.excluded.is_empty() {
        println!();
        println!("excluded:");
        println!("  sample_id\treason");
        for sample in &comparison.excluded {
            println!("  {}\t{}", sample.sample_id, sample.reason);
        }
    }
    println!();

    if comparison.included.is_empty() {
        return Err(format!(
            "no eligible samples for gene {}: every sample was excluded; check that junction sample_id values match the manifest and that this gene has junction rows",
            comparison.gene_id
        ));
    }

    let case_count = comparison
        .included
        .iter()
        .filter(|sample| sample.condition == splice_girl::manifest::Condition::Case)
        .count();
    let control_count = comparison.included.len() - case_count;
    if case_count == 0 || control_count == 0 {
        let missing_group = if case_count == 0 { "case" } else { "control" };
        return Err(format!(
            "cannot compare gene {}: no eligible {missing_group} samples remain; add a matching {missing_group} sample or review the exclusions above",
            comparison.gene_id
        ));
    }

    Ok(())
}

fn print_validation_summary(summary: &ValidationSummary) -> Result<(), String> {
    println!("validation summary");
    println!("  junction_rows: {}", summary.junction_rows);
    println!("  samples: {}", summary.samples);
    println!("  genes: {}", summary.genes);
    println!("  sample_gene_groups: {}", summary.sample_gene_groups);
    println!("  total_read_support: {:.10}", summary.total_read_support);
    println!(
        "  group_support_range: {:.10} to {:.10}",
        summary.minimum_group_support, summary.maximum_group_support
    );
    println!("  mathematical_checks_passed: {}", summary.checks_passed);
    println!("  mathematical_checks_failed: {}", summary.checks_failed);

    if let Some(manifest) = &summary.manifest {
        println!("manifest summary");
        println!("  manifest_samples: {}", manifest.manifest_samples);
        println!("  case_samples: {}", manifest.case_samples);
        println!("  control_samples: {}", manifest.control_samples);
        println!(
            "  junction_samples_without_metadata: {}",
            display_list(&manifest.junction_samples_without_metadata)
        );
        println!(
            "  manifest_samples_without_junctions: {}",
            display_list(&manifest.manifest_samples_without_junctions)
        );
    }

    let math_failed = summary.checks_failed > 0;
    let manifest_failed = summary.manifest.as_ref().is_some_and(|manifest| {
        manifest.case_samples == 0
            || manifest.control_samples == 0
            || !manifest.junction_samples_without_metadata.is_empty()
            || !manifest.manifest_samples_without_junctions.is_empty()
    });
    if math_failed || manifest_failed {
        println!("validation_status: FAIL");
        return Err(
            "validation found failed mathematical checks or sample/manifest mismatches; review the named counts and sample IDs above"
                .to_owned(),
        );
    }

    println!("validation_status: PASS");
    println!();
    Ok(())
}

fn print_all_gene_comparisons(summaries: &[GeneComparisonSummary]) -> Result<(), String> {
    println!("all-gene descriptive comparison");
    println!(
        "gene_id\tstatus\tcase_n\tcontrol_n\tcase_mean\tcontrol_mean\tmean_difference_case_minus_control\tcase_median\tcontrol_median\tcase_range\tcontrol_range\texcluded"
    );
    for summary in summaries {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            summary.comparison.gene_id,
            summary.status,
            summary.case.as_ref().map_or(0, |group| group.n),
            summary.control.as_ref().map_or(0, |group| group.n),
            terminal_number(summary.case.as_ref().map(|group| group.mean)),
            terminal_number(summary.control.as_ref().map(|group| group.mean)),
            terminal_number(summary.mean_difference),
            terminal_number(summary.case.as_ref().map(|group| group.median)),
            terminal_number(summary.control.as_ref().map(|group| group.median)),
            terminal_range(summary.case.as_ref()),
            terminal_range(summary.control.as_ref()),
            summary.comparison.excluded.len(),
        );
    }

    let complete = summaries
        .iter()
        .filter(|summary| summary.status == ComparisonStatus::Complete)
        .count();
    println!();
    println!("genes_total: {}", summaries.len());
    println!("genes_with_both_groups: {complete}");
    println!("effect_definition: case mean minus control mean");
    println!("inference: descriptive only; no p-values, q-values, or ranking calculated");
    println!();

    if complete == 0 {
        return Err(
            "no gene has at least one eligible case and one eligible control sample; review each gene status and the sample manifest"
                .to_owned(),
        );
    }
    Ok(())
}

fn terminal_number(value: Option<f64>) -> String {
    value
        .map(|number| format!("{number:.10}"))
        .unwrap_or_else(|| "NA".to_owned())
}

fn terminal_range(group: Option<&splice_girl::compare::ConditionSummary>) -> String {
    group
        .map(|group| format!("{:.10}..{:.10}", group.minimum, group.maximum))
        .unwrap_or_else(|| "NA".to_owned())
}

fn display_list(values: &[String]) -> String {
    if values.is_empty() {
        "none".to_owned()
    } else {
        values.join(", ")
    }
}

fn print_export_summary(summary: &ExportSummary) {
    println!("results saved: {}", summary.output_dir.display());
    println!("run_id: {}", summary.run_id);
    for file in &summary.files {
        println!("  {}", file.display());
    }
    println!();
}
