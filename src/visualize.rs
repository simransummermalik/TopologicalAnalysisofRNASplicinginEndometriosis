//! Dependency-free SVG rendering for one analyzed sample/gene graph.
//!
//! The renderer uses the same `GroupResult` produced by the validated
//! analysis path. It does not recalculate or reinterpret the decomposition;
//! it presents the graph and its values as a plain academic figure.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

use crate::GroupResult;

const MIN_WIDTH: f64 = 900.0;
const NODE_WIDTH: f64 = 100.0;
const NODE_HEIGHT: f64 = 42.0;
const NODE_Y: f64 = 330.0;
const NODE_SPACING: f64 = 150.0;
const SIDE_MARGIN: f64 = 105.0;

/// Render one analyzed splice graph as a standalone SVG document.
pub fn render_svg(result: &GroupResult) -> String {
    let graph = &result.graph;
    let decomposition = &result.decomposition;
    let width = MIN_WIDTH
        .max(SIDE_MARGIN * 2.0 + NODE_SPACING * graph.nodes.len().saturating_sub(1) as f64);
    let first_table_row_y = 535.0;
    let table_row_height = 26.0;
    let table_end_y = first_table_row_y + table_row_height * graph.edges.len() as f64;
    let footer_line_y = table_end_y + 22.0;
    let canvas_height = footer_line_y + 76.0;
    let max_count = graph
        .edges
        .iter()
        .map(|edge| edge.count)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let mut svg = String::new();
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{canvas_height:.0}" viewBox="0 0 {width:.0} {canvas_height:.0}" role="img" aria-labelledby="title description">
<title id="title">Splice graph for sample {} and gene {}</title>
<desc id="description">A directed splice graph. Line width represents junction read support. The table reports normalized signal and the signed Hodge cycle component.</desc>
<defs>
  <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#111111"/>
  </marker>
  <style>
    .serif {{ font-family: Georgia, "Times New Roman", serif; }}
    .mono {{ font-family: "Courier New", Courier, monospace; }}
  </style>
</defs>
<rect width="100%" height="100%" fill="#ffffff"/>

<text class="serif" x="55" y="48" font-size="24" font-weight="700" fill="#111111">Splice graph: {gene}</text>
<text class="serif" x="55" y="78" font-size="15" fill="#111111">Sample: {sample}</text>
<line x1="55" y1="96" x2="{rule_end:.1}" y2="96" stroke="#111111" stroke-width="1"/>
<text class="serif" x="55" y="121" font-size="14" fill="#111111">Nodes: {nodes}   ·   Junctions: {edges}   ·   Cycle fraction: {cycle_fraction:.4}   ·   Mathematical checks: {status}</text>

<text class="serif" x="55" y="164" font-size="17" font-weight="700" fill="#111111">Observed junction graph</text>
<text class="serif" x="55" y="188" font-size="13" font-style="italic" fill="#333333">Arrow direction follows the junction data; thicker arrows have greater read support.</text>
"##,
        xml_escape(&result.sample_id),
        xml_escape(&result.gene_id),
        gene = xml_escape(&result.gene_id),
        sample = xml_escape(&result.sample_id),
        rule_end = width - 55.0,
        nodes = graph.nodes.len(),
        edges = graph.edges.len(),
        cycle_fraction = decomposition.cycle_fraction,
        status = if decomposition.checks.passed {
            "PASS"
        } else {
            "FAIL"
        },
    ));

    // Edges are drawn first so the node boxes remain legible above them.
    for (index, edge) in graph.edges.iter().enumerate() {
        let from_x = node_x(edge.from);
        let to_x = node_x(edge.to);
        let direction = if to_x >= from_x { 1.0 } else { -1.0 };
        let start_x = from_x + direction * NODE_WIDTH / 2.0;
        let end_x = to_x - direction * NODE_WIDTH / 2.0;
        let midpoint_x = (start_x + end_x) / 2.0;
        let span = edge.from.abs_diff(edge.to).max(1) as f64;
        let arc_height = (55.0 + 24.0 * span).min(120.0);
        let control_y = if direction > 0.0 {
            NODE_Y - arc_height
        } else {
            NODE_Y + arc_height
        };
        let label_y = (NODE_Y + control_y) / 2.0 + if direction > 0.0 { -4.0 } else { 15.0 };
        let stroke_width = 0.9 + 2.1 * (edge.count / max_count).sqrt();

        svg.push_str(&format!(
            r##"<path d="M {start_x:.1} {NODE_Y:.1} Q {midpoint_x:.1} {control_y:.1} {end_x:.1} {NODE_Y:.1}" fill="none" stroke="#111111" stroke-width="{stroke_width:.2}" stroke-linecap="round" marker-end="url(#arrow)"/>
<rect x="{badge_x:.1}" y="{badge_y:.1}" width="28" height="19" fill="#ffffff" stroke="#111111" stroke-width="0.8"/>
<text class="serif" x="{midpoint_x:.1}" y="{badge_text_y:.1}" text-anchor="middle" font-size="11" fill="#111111">J{number}</text>
"##,
            badge_x = midpoint_x - 14.0,
            badge_y = label_y - 13.0,
            badge_text_y = label_y,
            number = index + 1,
        ));
    }

    for (index, node) in graph.nodes.iter().enumerate() {
        let x = node_x(index);
        let display = shorten(node, 15);
        svg.push_str(&format!(
            r##"<rect x="{left:.1}" y="{top:.1}" width="{NODE_WIDTH:.1}" height="{NODE_HEIGHT:.1}" fill="#ffffff" stroke="#111111" stroke-width="1.4"/>
<text class="serif" x="{x:.1}" y="{text_y:.1}" text-anchor="middle" font-size="14" fill="#111111">{display}<title>{full}</title></text>
"##,
            left = x - NODE_WIDTH / 2.0,
            top = NODE_Y - NODE_HEIGHT / 2.0,
            text_y = NODE_Y + 5.0,
            display = xml_escape(&display),
            full = xml_escape(node),
        ));
    }

    svg.push_str(&format!(
        r##"<text class="serif" x="55" y="482" font-size="17" font-weight="700" fill="#111111">Junction values</text>
<line x1="55" y1="497" x2="{table_end:.1}" y2="497" stroke="#111111" stroke-width="1"/>
<text class="serif" x="68" y="518" font-size="12" font-weight="700" fill="#111111">ID</text>
<text class="serif" x="120" y="518" font-size="12" font-weight="700" fill="#111111">Junction</text>
<text class="serif" x="350" y="518" font-size="12" font-weight="700" fill="#111111">Read support</text>
<text class="serif" x="500" y="518" font-size="12" font-weight="700" fill="#111111">Normalized signal</text>
<text class="serif" x="675" y="518" font-size="12" font-weight="700" fill="#111111">Cycle component</text>
<line x1="55" y1="524" x2="{table_end:.1}" y2="524" stroke="#777777" stroke-width="0.7"/>
"##,
        table_end = width - 55.0,
    ));

    for (index, edge) in graph.edges.iter().enumerate() {
        let row_y = first_table_row_y + table_row_height * index as f64;
        let from = shorten(&graph.nodes[edge.from], 24);
        let to = shorten(&graph.nodes[edge.to], 24);
        let connection = format!("{from} → {to}");
        svg.push_str(&format!(
            r##"<text class="serif" x="68" y="{row_y:.1}" font-size="12" fill="#111111">J{number}</text>
<text class="serif" x="120" y="{row_y:.1}" font-size="12" fill="#111111">{connection}</text>
<text class="mono" x="350" y="{row_y:.1}" font-size="11" fill="#111111">{support}</text>
<text class="mono" x="500" y="{row_y:.1}" font-size="11" fill="#111111">{signal:.4}</text>
<text class="mono" x="675" y="{row_y:.1}" font-size="11" fill="#111111">{cycle:+.4}</text>
<line x1="55" y1="{line_y:.1}" x2="{table_end:.1}" y2="{line_y:.1}" stroke="#cccccc" stroke-width="0.5"/>
"##,
            number = index + 1,
            connection = xml_escape(&connection),
            support = xml_escape(&format_count(edge.count)),
            signal = graph.signal[index],
            cycle = decomposition.cycle[index],
            line_y = row_y + 8.0,
            table_end = width - 55.0,
        ));
    }

    svg.push_str(&format!(
        r##"<line x1="55" y1="{footer_line_y:.1}" x2="{footer_end:.1}" y2="{footer_line_y:.1}" stroke="#111111" stroke-width="1"/>
<text class="serif" x="55" y="{footer_text_y:.1}" font-size="12" fill="#222222">Cycle fraction = cycle energy / total signal energy. Component signs follow the stored edge direction.</text>
<text class="serif" x="55" y="{caveat_y:.1}" font-size="12" font-style="italic" fill="#444444">Cycle-space structure alone does not establish circular RNA or disease causation.</text>
</svg>
"##,
        footer_end = width - 55.0,
        footer_text_y = footer_line_y + 28.0,
        caveat_y = footer_line_y + 50.0,
    ));

    svg
}

/// Write a new SVG without silently replacing an existing figure.
pub fn write_svg(path: &Path, result: &GroupResult) -> Result<(), String> {
    let svg = render_svg(result);
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| {
            if error.kind() == std::io::ErrorKind::AlreadyExists {
                format!(
                    "output {} already exists; choose another filename or remove it first",
                    path.display()
                )
            } else {
                format!("could not create {}: {error}", path.display())
            }
        })?;
    file.write_all(svg.as_bytes())
        .map_err(|error| format!("could not write {}: {error}", path.display()))
}

fn node_x(index: usize) -> f64 {
    SIDE_MARGIN + index as f64 * NODE_SPACING
}

fn format_count(count: f64) -> String {
    if count.fract().abs() <= f64::EPSILON {
        format!("{count:.0}")
    } else {
        format!("{count:.2}")
    }
}

fn shorten(value: &str, max_chars: usize) -> String {
    if value.chars().count() <= max_chars {
        return value.to_owned();
    }
    let prefix: String = value.chars().take(max_chars.saturating_sub(1)).collect();
    format!("{prefix}…")
}

fn xml_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analyze_junctions;
    use crate::parser::parse_tsv_text;

    #[test]
    fn svg_contains_graph_metrics_and_junction_values() {
        let rows = parse_tsv_text(
            "sample_id\tgene_id\tfrom\tto\tcount\n\
             case_01\tGENE_DEMO\tE1\tE2\t100\n\
             case_01\tGENE_DEMO\tE2\tE3\t100\n\
             case_01\tGENE_DEMO\tE1\tE3\t100\n",
        )
        .unwrap();
        let results = analyze_junctions(rows).unwrap();
        let svg = render_svg(&results[0]);

        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("Splice graph: GENE_DEMO"));
        assert!(svg.contains("Sample: case_01"));
        assert!(svg.contains("Cycle fraction: 0.1111"));
        assert!(svg.contains("Junction values"));
        assert!(svg.contains("E1 → E2"));
        assert!(svg.contains(">100</text>"));
        assert!(svg.contains("marker-end=\"url(#arrow)\""));
        assert!(!svg.contains("linearGradient"));
        assert!(!svg.contains("feDropShadow"));
        assert!(svg.ends_with("</svg>\n"));
    }

    #[test]
    fn svg_escapes_identifiers_and_preserves_full_node_in_title() {
        let rows = parse_tsv_text(
            "sample_id\tgene_id\tfrom\tto\tcount\n\
             sample&amp\tgene<1>\ta-very-long-node-name\tE2\t4\n",
        )
        .unwrap();
        let results = analyze_junctions(rows).unwrap();
        let svg = render_svg(&results[0]);

        assert!(svg.contains("Sample: sample&amp;amp"));
        assert!(svg.contains("gene&lt;1&gt;"));
        assert!(svg.contains("<title>a-very-long-node-name</title>"));
        assert!(!svg.contains("gene<1>"));
    }
}
