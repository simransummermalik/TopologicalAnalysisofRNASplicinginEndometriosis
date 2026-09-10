use std::env;
use std::path::Path;
use std::process::ExitCode;

use splice_girl::analyze_junctions;
use splice_girl::parser::parse_tsv;

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let mut arguments = env::args();
    let program = arguments.next().unwrap_or_else(|| "splice-girl".to_owned());
    let path = arguments
        .next()
        .ok_or_else(|| format!("usage: {program} <junctions.tsv>"))?;
    if arguments.next().is_some() {
        return Err(format!("usage: {program} <junctions.tsv>"));
    }

    let junctions = parse_tsv(Path::new(&path))?;
    let results = analyze_junctions(junctions)?;

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
