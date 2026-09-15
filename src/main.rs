//! CLI entry point for `splice-girl`. The basic direct command accepts a
//! junction TSV, while `--start` opens the small interactive menu around that
//! same analysis feature. The `analyze`/`gene`/`compare`/`rank` subcommands
//! described in ROADMAP.md are not built yet.
use std::env;
use std::io::{self, IsTerminal, Write};
use std::path::Path;
use std::process::ExitCode;

use splice_girl::analyze_junctions;
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

fn run() -> Result<(), String> {
    let mut arguments = env::args();
    let program = arguments.next().unwrap_or_else(|| "splice-girl".to_owned());
    let first_argument = arguments
        .next()
        .ok_or_else(|| format!("usage: {program} [--start | <junctions.tsv>]"))?;

    if first_argument == "--start" {
        if arguments.next().is_some() {
            return Err(format!("usage: {program} [--start | <junctions.tsv>]"));
        }
        return run_start_menu();
    }

    if arguments.next().is_some() {
        return Err(format!("usage: {program} [--start | <junctions.tsv>]"));
    }

    let results = analyze_path(Path::new(&first_argument))?;
    print_results(results)
}

fn analyze_path(path: &Path) -> Result<Vec<splice_girl::GroupResult>, String> {
    let junctions = parse_tsv(path)?;
    analyze_junctions(junctions)
}

fn run_start_menu() -> Result<(), String> {
    loop {
        println!();
        println!("  Splice Girl");
        println!("  -----------");
        println!("  [1] Analyze a junction TSV");
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
            "q" => return Ok(()),
            "" => {}
            _ => println!("  Please choose 1 or q."),
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
