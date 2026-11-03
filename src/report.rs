//! Offline, academic comparison reports built from the existing analysis.
//! Figures use individual sample values; missing measurements stay missing.

use std::collections::BTreeMap;
use std::fmt::Write;

use crate::GroupResult;
use crate::compare::{ComparisonStatus, GeneComparisonSummary};
use crate::manifest::Condition;

/// Relative filenames use numeric indices, never untrusted identifiers.
pub fn graph_filename(index: usize) -> String {
    format!("graphs/graph-{:06}.svg", index + 1)
}

/// Render the document without scripts, external fonts, or network requests.
/// `run_id` ties the report to its accompanying machine-readable bundle.
pub fn render_report(
    run_id: &str,
    results: &[GroupResult],
    summaries: &[GeneComparisonSummary],
) -> String {
    let graph_indices: BTreeMap<(&str, &str), usize> = results
        .iter()
        .enumerate()
        .filter(|(_, result)| result.decomposition.checks.passed)
        .map(|(index, result)| ((result.gene_id.as_str(), result.sample_id.as_str()), index))
        .collect();
    let complete = summaries
        .iter()
        .filter(|summary| summary.status == ComparisonStatus::Complete)
        .count();
    let failed = results
        .iter()
        .filter(|result| !result.decomposition.checks.passed)
        .count();
    let mut html = String::from(HEADER);
    writeln!(
        html,
        "<p class=metadata>Version {} · Run <code>{}</code></p>",
        env!("CARGO_PKG_VERSION"),
        escape(run_id),
    )
    .unwrap();
    writeln!(
        html,
        "<p>{} genes; {complete} with both groups; {} measured sample/gene graphs; {failed} failed mathematical checks.</p>",
        summaries.len(), results.len(),
    )
    .unwrap();
    html.push_str(
        "<p>Each point represents one eligible sample. Cycle fraction measures how much normalized junction signal remains in the graph cycle space after removing the gradient component. The effect is the case mean minus the control mean.</p>\n\
         <p class=note>Descriptive results only. Mathematical checks verify the decomposition; they do not establish biological comparability or statistical significance. Check tissue, coverage, and graph support before interpretation. Cycle-space signal alone does not establish circular RNA or causation. Demo inputs produce synthetic demonstration results.</p>\n\
         <nav aria-label=\"Result files\">Files: <a href=\"sample_metrics.tsv\">sample metrics</a> · <a href=\"junction_components.tsv\">junction components</a> · <a href=\"comparisons.tsv\">sample decisions</a> · <a href=\"gene_summaries.tsv\">gene summaries</a> · <a href=\"run_manifest.json\">provenance and input checksums</a></nav>\n\
         <h2>Gene summary</h2><div class=table-wrap><table>\n\
         <caption>Genes in identifier order; NA means the effect cannot be calculated.</caption>\n\
         <thead><tr><th scope=col>Gene</th><th scope=col>Status</th><th scope=col>Control n</th><th scope=col>Case n</th><th scope=col>Control mean</th><th scope=col>Case mean</th><th scope=col>Case − control</th><th scope=col>Excluded</th></tr></thead><tbody>\n",
    );
    for (index, summary) in summaries.iter().enumerate() {
        writeln!(
            html,
            "<tr><th scope=row><a href=\"#gene-{}\">{}</a></th><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            index + 1,
            escape(&summary.comparison.gene_id),
            summary.status,
            summary.control.as_ref().map_or(0, |group| group.n),
            summary.case.as_ref().map_or(0, |group| group.n),
            number(summary.control.as_ref().map(|group| group.mean)),
            number(summary.case.as_ref().map(|group| group.mean)),
            number(summary.mean_difference),
            summary.comparison.excluded.len(),
        )
        .unwrap();
    }
    html.push_str("</tbody></table></div>\n<h2>Individual samples and graphs</h2>\n");
    for (index, summary) in summaries.iter().enumerate() {
        let comparison = &summary.comparison;
        writeln!(
            html,
            "<section id=\"gene-{}\"><h3>{}</h3><p>Status: <code>{}</code>. Mean difference (case − control): <strong>{}</strong>.</p>",
            index + 1,
            escape(&comparison.gene_id),
            summary.status,
            number(summary.mean_difference),
        )
        .unwrap();
        html.push_str(&sample_plot(summary));
        html.push_str("<div class=table-wrap><table><caption>Eligible samples. Each graph link opens its observed, gradient, and cycle panels.</caption><thead><tr><th scope=col>Sample</th><th scope=col>Condition</th><th scope=col>Cycle fraction</th><th scope=col>Graph</th></tr></thead><tbody>\n");
        for sample in &comparison.included {
            let graph = graph_indices
                .get(&(comparison.gene_id.as_str(), sample.sample_id.as_str()))
                .map(|index| {
                    format!(
                        "<a href=\"{}\">View decomposition</a>",
                        graph_filename(*index)
                    )
                })
                .unwrap_or_else(|| "Unavailable".to_owned());
            writeln!(
                html,
                "<tr><th scope=row>{}</th><td>{}</td><td>{:.10}</td><td>{graph}</td></tr>",
                escape(&sample.sample_id),
                sample.condition,
                sample.cycle_fraction,
            )
            .unwrap();
        }
        html.push_str("</tbody></table></div>\n");
        if comparison.included.is_empty() {
            html.push_str("<p>No eligible samples; no sample values or effect were imputed.</p>\n");
        }
        if !comparison.excluded.is_empty() {
            writeln!(html, "<details><summary>Excluded samples ({})</summary><table><thead><tr><th scope=col>Sample</th><th scope=col>Condition</th><th scope=col>Reason</th></tr></thead><tbody>", comparison.excluded.len()).unwrap();
            for sample in &comparison.excluded {
                writeln!(
                    html,
                    "<tr><th scope=row>{}</th><td>{}</td><td>{}</td></tr>",
                    escape(&sample.sample_id),
                    sample
                        .condition
                        .map(|value| value.to_string())
                        .unwrap_or_else(|| "unknown".to_owned()),
                    escape(&sample.reason),
                )
                .unwrap();
            }
            html.push_str("</tbody></table></details>\n");
        }
        html.push_str("<p><a href=\"#top\">Back to summary</a></p></section>\n");
    }
    html.push_str("<footer>Generated by Splice Girl from this run’s input tables. Keep this entire folder together when sharing. Exact values, settings, exclusions, and SHA-256 input checksums are in the linked TSV and JSON files.</footer></main></body></html>\n");
    html
}

fn number(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.10}"))
        .unwrap_or_else(|| "NA".to_owned())
}

fn sample_plot(summary: &GeneComparisonSummary) -> String {
    let samples = &summary.comparison.included;
    if samples.is_empty() {
        return String::new();
    }
    let height = 95 + samples.len() * 28;
    let axis_y = height - 45;
    let mut svg = format!(
        "<figure><svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 850 {height}\" role=\"img\" aria-label=\"Individual sample cycle fractions for {}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><g font-family=\"Georgia,serif\" font-size=\"13\" fill=\"#111\">\n",
        escape(&summary.comparison.gene_id),
    );
    for (index, sample) in samples.iter().enumerate() {
        let y = 28 + index * 28;
        let x = 290.0 + sample.cycle_fraction * 400.0;
        let mut label: String = sample.sample_id.chars().take(26).collect();
        if sample.sample_id.chars().count() > 26 {
            label.push('…');
        }
        writeln!(svg, "<text x=\"10\" y=\"{}\"><title>{}</title>{} ({})</text><line x1=\"290\" y1=\"{y}\" x2=\"690\" y2=\"{y}\" stroke=\"#ddd\"/><text x=\"710\" y=\"{}\">{:.4}</text>", y + 4, escape(&sample.sample_id), escape(&label), sample.condition, y + 4, sample.cycle_fraction).unwrap();
        let title = format!(
            "{}: {:.10}",
            escape(&sample.sample_id),
            sample.cycle_fraction
        );
        if sample.condition == Condition::Case {
            writeln!(
                svg,
                "<circle cx=\"{x:.4}\" cy=\"{y}\" r=\"5\" fill=\"#111\"><title>{title}</title></circle>"
            )
            .unwrap();
        } else {
            writeln!(svg, "<rect x=\"{:.4}\" y=\"{}\" width=\"9\" height=\"9\" fill=\"white\" stroke=\"#111\"><title>{title}</title></rect>", x - 4.5, y as f64 - 4.5).unwrap();
        }
    }
    writeln!(
        svg,
        "<line x1=\"290\" y1=\"{axis_y}\" x2=\"690\" y2=\"{axis_y}\" stroke=\"#111\"/>"
    )
    .unwrap();
    for tick in 0..=4 {
        let x = 290 + tick * 100;
        writeln!(svg, "<line x1=\"{x}\" y1=\"{axis_y}\" x2=\"{x}\" y2=\"{}\" stroke=\"#111\"/><text x=\"{x}\" y=\"{}\" text-anchor=\"middle\">{:.2}</text>", axis_y + 5, axis_y + 21, tick as f64 / 4.0).unwrap();
    }
    svg.push_str("</g></svg><figcaption>Cycle fraction, on the same 0–1 scale for every gene. Filled circles: cases. Open squares: controls. Points use individual rows so overlapping values remain visible. Exact values appear below.</figcaption></figure>\n");
    svg
}

fn escape(value: &str) -> String {
    value.chars().fold(String::new(), |mut text, character| {
        text.push_str(match character {
            '&' => "&amp;",
            '<' => "&lt;",
            '>' => "&gt;",
            '"' => "&quot;",
            '\'' => "&#39;",
            _ => {
                text.push(character);
                return text;
            }
        });
        text
    })
}

const HEADER: &str = r##"<!doctype html>
<html lang=en><head><meta charset=utf-8><meta name=viewport content="width=device-width, initial-scale=1">
<title>Splice Girl — comparison report</title>
<style>
body{margin:0;background:white;color:#111;font-family:Georgia,"Times New Roman",serif;line-height:1.6}
main{max-width:1050px;margin:3rem auto;padding:0 1.2rem}h1{font-size:2rem;margin-bottom:.2rem}h2{margin-top:2rem;border-bottom:1px solid #111}h3{font-size:1.35rem}
a{color:#111;text-underline-offset:.15em}a:hover{text-decoration-thickness:2px}code{font-family:monospace;overflow-wrap:anywhere}.metadata,figcaption,footer{color:#444;font-size:.9rem}.note{border-left:2px solid #777;padding-left:1rem}
table{border-collapse:collapse;width:100%;font-size:.9rem;margin:.8rem 0}th,td{text-align:left;padding:.4rem .55rem;border-bottom:1px solid #ddd;overflow-wrap:anywhere}thead{border-top:1px solid #111;border-bottom:1px solid #111}caption{text-align:left;color:#444;margin-bottom:.4rem}.table-wrap{overflow-x:auto}
section{padding:1rem 0;border-bottom:1px solid #bbb;scroll-margin-top:1rem}figure{margin:1rem 0}svg{width:100%;height:auto}details{margin:1rem 0}summary{cursor:pointer}nav{padding:.8rem 0}footer{margin-top:2rem}
@media print{main{margin:0;max-width:none}section{break-before:page}details{display:block}details>*{display:block}a{color:inherit}figure{break-inside:avoid}}
</style></head><body><main id=top><h1>Splice Girl</h1><p>Splice graph and Hodge comparison report</p>
"##;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        analyze_junctions, compare::compare_all_genes, manifest::parse_manifest_text,
        parser::parse_tsv_text,
    };

    #[test]
    fn report_preserves_missing_values_and_escapes_input_as_text() {
        let rows =
            parse_tsv_text("sample_id\tgene_id\tfrom\tto\tcount\n<script>\t../G&\tE1\tE2\t4\n")
                .unwrap();
        let results = analyze_junctions(rows).unwrap();
        let manifest =
            parse_manifest_text("sample_id\tcondition\n<script>\tcase\nabsent\tcontrol\n").unwrap();
        let summaries = compare_all_genes(&manifest, &results);
        let report = render_report("test", &results, &summaries);
        assert!(report.contains("&lt;script&gt;"));
        assert!(!report.contains("<script>"));
        assert!(report.contains("../G&amp;"));
        assert!(report.contains("missing_control"));
        assert!(report.contains("<strong>NA</strong>"));
        assert!(report.contains("no junction data for gene"));
        assert!(report.contains("href=\"graphs/graph-000001.svg\""));
        assert!(!report.contains("href=\"../G"));
    }
}
