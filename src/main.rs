//! CLI entry point for `splice-girl`. The direct command accepts a junction
//! TSV, `gene` compares one gene across a manifest's case/control samples,
//! and `--start` opens a small interactive menu around both. ROADMAP.md's
//! `gene` syntax is directory-based (`--case dir/ --control dir/`); this
//! implementation instead takes one junction TSV (already carrying
//! `sample_id`) plus a separate manifest TSV, matching the formats this
//! tool actually parses. The `analyze`/`compare`/`rank` subcommands are
//! not built yet.
use std::env;
use std::io::{self, IsTerminal, Write};
use std::path::Path;
use std::process::ExitCode;

use splice_girl::analyze_junctions;
use splice_girl::compare::{GeneComparison, compare_gene};
use splice_girl::manifest::parse_manifest;
use splice_girl::parser::parse_tsv;

fn main() -> ExitCode {
    print_banner();
    match run() {
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
        "usage: {program} [--start | <junctions.tsv> | gene <GENE_ID> <junctions.tsv> <manifest.tsv>]"
    )
}

fn run() -> Result<(), String> {
    let mut arguments = env::args();
    let program = arguments.next().unwrap_or_else(|| "splice-girl".to_owned());
    let first_argument = arguments.next().ok_or_else(|| usage(&program))?;

    if first_argument == "--start" {
        if arguments.next().is_some() {
            return Err(usage(&program));
        }
        return run_start_menu();
    }

    if first_argument == "gene" {
        let gene_id = arguments.next().ok_or_else(|| usage(&program))?;
        let junctions_path = arguments.next().ok_or_else(|| usage(&program))?;
        let manifest_path = arguments.next().ok_or_else(|| usage(&program))?;
        if arguments.next().is_some() {
            return Err(usage(&program));
        }
        let comparison = run_gene_comparison(
            &gene_id,
            Path::new(&junctions_path),
            Path::new(&manifest_path),
        )?;
        return print_gene_comparison(comparison);
    }

    if arguments.next().is_some() {
        return Err(usage(&program));
    }

    let results = analyze_path(Path::new(&first_argument))?;
    print_results(results)
}

fn analyze_path(path: &Path) -> Result<Vec<splice_girl::GroupResult>, String> {
    let junctions = parse_tsv(path)?;
    analyze_junctions(junctions)
}

fn run_gene_comparison(
    gene_id: &str,
    junctions_path: &Path,
    manifest_path: &Path,
) -> Result<GeneComparison, String> {
    let manifest = parse_manifest(manifest_path)?;
    let results = analyze_path(junctions_path)?;
    Ok(compare_gene(&manifest, &results, gene_id))
}

fn run_start_menu() -> Result<(), String> {
    loop {
        println!();
        println!("  Splice Girl");
        println!("  -----------");
        println!("  [1] Analyze a junction TSV");
        println!("  [2] Compare one gene across case/control samples");
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

                match analyze_path(Path::new(path)).and_then(print_results) {
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

                let outcome = run_gene_comparison(
                    gene_id,
                    Path::new(junctions_path),
                    Path::new(manifest_path),
                )
                .and_then(print_gene_comparison);
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
            _ => println!("  Please choose 1, 2, or q."),
        }
    }
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

fn print_results(results: Vec<splice_girl::GroupResult>) -> Result<(), String> {
    // One block of metrics per (sample, gene) group found in the input.
    for result in results {
        let decomposition = result.decomposition;
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

fn print_gene_comparison(comparison: GeneComparison) -> Result<(), String> {
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
            "no eligible samples for gene {}: every sample was excluded",
            comparison.gene_id
        ));
    }

    Ok(())
}
