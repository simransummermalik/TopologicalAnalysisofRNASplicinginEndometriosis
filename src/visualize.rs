//! Dependency-free SVG rendering for one analyzed sample/gene graph.
//!
//! The renderer uses the same `GroupResult` produced by the validated
//! analysis path. It does not recalculate or reinterpret the decomposition;
//! it presents the graph and its values as a plain academic figure.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

use crate::GroupResult;

const MIN_PANEL_WIDTH: f64 = 320.0;
const NODE_WIDTH: f64 = 84.0;
const NODE_HEIGHT: f64 = 42.0;
const NODE_Y: f64 = 405.0;
const NODE_SPACING: f64 = 130.0;
const PANEL_PADDING: f64 = 50.0;
const PANEL_GAP: f64 = 24.0;
const PANEL_TOP: f64 = 180.0;
const PANEL_BOTTOM: f64 = 545.0;
const OUTER_MARGIN: f64 = 48.0;

/// Render one analyzed splice graph as a standalone SVG document.
pub fn render_svg(result: &GroupResult) -> String {
    let graph = &result.graph;
    let decomposition = &result.decomposition;
    let graph_span = NODE_SPACING * graph.nodes.len().saturating_sub(1) as f64;
    let panel_width = MIN_PANEL_WIDTH.max(graph_span + PANEL_PADDING * 2.0);
    let width = OUTER_MARGIN * 2.0 + panel_width * 3.0 + PANEL_GAP * 2.0;
    let first_table_row_y = 650.0;
    let table_row_height = 26.0;
    let table_end_y = first_table_row_y + table_row_height * graph.edges.len() as f64;
    let footer_line_y = table_end_y + 22.0;
    let canvas_height = footer_line_y + 76.0;
    let component_scale = graph
        .signal
        .iter()
        .chain(&decomposition.gradient)
        .chain(&decomposition.cycle)
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max)
        .max(f64::EPSILON);

    let mut svg = String::new();
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{canvas_height:.0}" viewBox="0 0 {width:.0} {canvas_height:.0}" role="img" aria-labelledby="title description">
<title id="title">Splice graph decomposition for sample {} and gene {}</title>
<desc id="description">Three aligned directed splice graphs show the observed normalized signal, gradient component, and cycle component. Every panel uses the same nodes, edges, direction, edge IDs, and line-width scale.</desc>
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

<text class="serif" x="48" y="42" font-size="24" font-weight="700" fill="#111111">Splice graph decomposition: {gene}</text>
<text class="serif" x="48" y="70" font-size="15" fill="#111111">Sample: {sample}</text>
<line x1="48" y1="88" x2="{rule_end:.1}" y2="88" stroke="#111111" stroke-width="1"/>
<text class="serif" x="48" y="113" font-size="14" fill="#111111">Nodes: {nodes}   ·   Junctions: {edges}   ·   Cycle fraction: {cycle_fraction:.4}   ·   Mathematical checks: {status}</text>
<text class="serif" x="48" y="143" font-size="13" fill="#111111">Same graph layout and edge IDs in every panel. Line width uses one shared absolute-value scale.</text>
<line x1="48" y1="163" x2="84" y2="163" stroke="#111111" stroke-width="2.2"/>
<text class="serif" x="92" y="167" font-size="12" fill="#111111">solid = positive</text>
<line x1="218" y1="163" x2="254" y2="163" stroke="#111111" stroke-width="2.2" stroke-dasharray="7 5"/>
<text class="serif" x="262" y="167" font-size="12" fill="#111111">dashed = negative</text>
<text class="serif" x="410" y="167" font-size="12" fill="#111111">Arrows always retain the stored biological junction direction.</text>
"##,
        xml_escape(&result.sample_id),
        xml_escape(&result.gene_id),
        gene = xml_escape(&result.gene_id),
        sample = xml_escape(&result.sample_id),
        rule_end = width - OUTER_MARGIN,
        nodes = graph.nodes.len(),
        edges = graph.edges.len(),
        cycle_fraction = decomposition.cycle_fraction,
        status = if decomposition.checks.passed {
            "PASS"
        } else {
            "FAIL"
        },
    ));

    let panels: [(&str, &str, &[f64]); 3] = [
        (
            "1. Observed signal",
            "Normalized junction usage (F)",
            &graph.signal,
        ),
        (
            "2. Gradient component",
            "Part explained by node values (F_grad)",
            &decomposition.gradient,
        ),
        (
            "3. Cycle component",
            "Remaining alternative-path structure (F_cycle)",
            &decomposition.cycle,
        ),
    ];
    for (panel_index, (title, subtitle, values)) in panels.into_iter().enumerate() {
        render_panel(
            &mut svg,
            result,
            panel_index,
            panel_width,
            title,
            subtitle,
            values,
            component_scale,
        );
    }

    svg.push_str(&format!(
        r##"<text class="serif" x="48" y="585" font-size="17" font-weight="700" fill="#111111">Junction values (same edge order in every panel)</text>
<line x1="48" y1="601" x2="{table_end:.1}" y2="601" stroke="#111111" stroke-width="1"/>
<text class="serif" x="60" y="628" font-size="12" font-weight="700" fill="#111111">ID</text>
<text class="serif" x="112" y="628" font-size="12" font-weight="700" fill="#111111">Junction</text>
<text class="serif" x="{support_x:.1}" y="628" font-size="12" font-weight="700" fill="#111111">Read support</text>
<text class="serif" x="{signal_x:.1}" y="628" font-size="12" font-weight="700" fill="#111111">Observed F</text>
<text class="serif" x="{gradient_x:.1}" y="628" font-size="12" font-weight="700" fill="#111111">Gradient F_grad</text>
<text class="serif" x="{cycle_x:.1}" y="628" font-size="12" font-weight="700" fill="#111111">Cycle F_cycle</text>
<line x1="48" y1="636" x2="{table_end:.1}" y2="636" stroke="#777777" stroke-width="0.7"/>
"##,
        table_end = width - OUTER_MARGIN,
        support_x = width * 0.39,
        signal_x = width * 0.55,
        gradient_x = width * 0.70,
        cycle_x = width * 0.85,
    ));

    for (index, edge) in graph.edges.iter().enumerate() {
        let row_y = first_table_row_y + table_row_height * index as f64;
        let from = shorten(&graph.nodes[edge.from], 24);
        let to = shorten(&graph.nodes[edge.to], 24);
        let connection = format!("{from} → {to}");
        svg.push_str(&format!(
            r##"<text class="serif" x="60" y="{row_y:.1}" font-size="12" fill="#111111">J{number}</text>
<text class="serif" x="112" y="{row_y:.1}" font-size="12" fill="#111111">{connection}</text>
<text class="mono" x="{support_x:.1}" y="{row_y:.1}" font-size="11" fill="#111111">{support}</text>
<text class="mono" x="{signal_x:.1}" y="{row_y:.1}" font-size="11" fill="#111111">{signal:+.4}</text>
<text class="mono" x="{gradient_x:.1}" y="{row_y:.1}" font-size="11" fill="#111111">{gradient:+.4}</text>
<text class="mono" x="{cycle_x:.1}" y="{row_y:.1}" font-size="11" fill="#111111">{cycle:+.4}</text>
<line x1="48" y1="{line_y:.1}" x2="{table_end:.1}" y2="{line_y:.1}" stroke="#cccccc" stroke-width="0.5"/>
"##,
            number = index + 1,
            connection = xml_escape(&connection),
            support = xml_escape(&format_count(edge.count)),
            signal = graph.signal[index],
            gradient = decomposition.gradient[index],
            cycle = decomposition.cycle[index],
            line_y = row_y + 8.0,
            table_end = width - OUTER_MARGIN,
            support_x = width * 0.39,
            signal_x = width * 0.55,
            gradient_x = width * 0.70,
            cycle_x = width * 0.85,
        ));
    }

    svg.push_str(&format!(
        r##"<line x1="48" y1="{footer_line_y:.1}" x2="{footer_end:.1}" y2="{footer_line_y:.1}" stroke="#111111" stroke-width="1"/>
<text class="serif" x="48" y="{footer_text_y:.1}" font-size="12" fill="#222222">Observed signal = gradient component + cycle component. Signs are relative to the stored edge direction.</text>
<text class="serif" x="48" y="{caveat_y:.1}" font-size="12" font-style="italic" fill="#444444">Cycle-space structure alone does not establish circular RNA or disease causation.</text>
</svg>
"##,
        footer_end = width - OUTER_MARGIN,
        footer_text_y = footer_line_y + 28.0,
        caveat_y = footer_line_y + 50.0,
    ));

    svg
}

#[allow(clippy::too_many_arguments)]
fn render_panel(
    svg: &mut String,
    result: &GroupResult,
    panel_index: usize,
    panel_width: f64,
    title: &str,
    subtitle: &str,
    values: &[f64],
    component_scale: f64,
) {
    let graph = &result.graph;
    let panel_x = OUTER_MARGIN + panel_index as f64 * (panel_width + PANEL_GAP);
    let panel_height = PANEL_BOTTOM - PANEL_TOP;
    svg.push_str(&format!(
        r##"<g aria-label="{title}">
<rect x="{panel_x:.1}" y="{PANEL_TOP:.1}" width="{panel_width:.1}" height="{panel_height:.1}" fill="#ffffff" stroke="#777777" stroke-width="0.8"/>
<text class="serif" x="{heading_x:.1}" y="207" font-size="17" font-weight="700" fill="#111111">{title}</text>
<text class="serif" x="{heading_x:.1}" y="230" font-size="12" fill="#333333">{subtitle}</text>
"##,
        title = xml_escape(title),
        heading_x = panel_x + 16.0,
        subtitle = xml_escape(subtitle),
    ));

    // Draw the same edges in the same order and positions in every panel.
    for (index, edge) in graph.edges.iter().enumerate() {
        let from_x = panel_node_x(panel_x, panel_width, graph.nodes.len(), edge.from);
        let to_x = panel_node_x(panel_x, panel_width, graph.nodes.len(), edge.to);
        let direction = if to_x >= from_x { 1.0 } else { -1.0 };
        let start_x = from_x + direction * NODE_WIDTH / 2.0;
        let end_x = to_x - direction * NODE_WIDTH / 2.0;
        let midpoint_x = (start_x + end_x) / 2.0;
        let span = edge.from.abs_diff(edge.to).max(1) as f64;
        let arc_height = (60.0 + 18.0 * span).min(115.0);
        let control_y = if direction > 0.0 {
            NODE_Y - arc_height
        } else {
            NODE_Y + arc_height
        };
        let label_y = control_y + if direction > 0.0 { 18.0 } else { -8.0 };
        let value = values[index];
        let relative_size = (value.abs() / component_scale).sqrt();
        let stroke_width = if value.abs() <= f64::EPSILON {
            0.7
        } else {
            1.0 + 3.0 * relative_size
        };
        let dash = if value < -f64::EPSILON { "7 5" } else { "none" };
        let label = format!("J{} {value:+.3}", index + 1);

        svg.push_str(&format!(
            r##"<path d="M {start_x:.1} {NODE_Y:.1} Q {midpoint_x:.1} {control_y:.1} {end_x:.1} {NODE_Y:.1}" fill="none" stroke="#111111" stroke-width="{stroke_width:.2}" stroke-dasharray="{dash}" stroke-linecap="round" marker-end="url(#arrow)"/>
<rect x="{badge_x:.1}" y="{badge_y:.1}" width="74" height="19" fill="#ffffff" stroke="#999999" stroke-width="0.6"/>
<text class="mono" x="{midpoint_x:.1}" y="{badge_text_y:.1}" text-anchor="middle" font-size="10" fill="#111111">{label}</text>
"##,
            badge_x = midpoint_x - 37.0,
            badge_y = label_y - 13.0,
            badge_text_y = label_y,
            label = xml_escape(&label),
        ));
    }

    for (index, node) in graph.nodes.iter().enumerate() {
        let x = panel_node_x(panel_x, panel_width, graph.nodes.len(), index);
        let display = shorten(node, 12);
        svg.push_str(&format!(
            r##"<rect x="{left:.1}" y="{top:.1}" width="{NODE_WIDTH:.1}" height="{NODE_HEIGHT:.1}" fill="#ffffff" stroke="#111111" stroke-width="1.2"/>
<text class="serif" x="{x:.1}" y="{text_y:.1}" text-anchor="middle" font-size="13" fill="#111111">{display}<title>{full}</title></text>
"##,
            left = x - NODE_WIDTH / 2.0,
            top = NODE_Y - NODE_HEIGHT / 2.0,
            text_y = NODE_Y + 5.0,
            display = xml_escape(&display),
            full = xml_escape(node),
        ));
    }
    svg.push_str("</g>\n");
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

fn panel_node_x(panel_x: f64, panel_width: f64, node_count: usize, index: usize) -> f64 {
    let graph_span = NODE_SPACING * node_count.saturating_sub(1) as f64;
    panel_x + (panel_width - graph_span) / 2.0 + index as f64 * NODE_SPACING
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
        assert!(svg.contains("Splice graph decomposition: GENE_DEMO"));
        assert!(svg.contains("Sample: case_01"));
        assert!(svg.contains("Cycle fraction: 0.1111"));
        assert!(svg.contains("1. Observed signal"));
        assert!(svg.contains("2. Gradient component"));
        assert!(svg.contains("3. Cycle component"));
        assert!(svg.contains("Junction values (same edge order in every panel)"));
        assert!(svg.contains("E1 → E2"));
        assert!(svg.contains(">100</text>"));
        assert!(svg.contains("J2 +0.444"));
        assert!(svg.contains("J2 -0.111"));
        assert!(svg.contains("Gradient F_grad"));
        assert!(svg.contains("Cycle F_cycle"));
        assert!(svg.contains("marker-end=\"url(#arrow)\""));
        assert!(svg.contains("stroke-dasharray=\"7 5\""));
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

    #[test]
    fn svg_values_change_with_the_selected_sample() {
        let rows = parse_tsv_text(
            "sample_id\tgene_id\tfrom\tto\tcount\n\
             cycle_sample\tGENE_A\tE1\tE2\t10\n\
             cycle_sample\tGENE_A\tE2\tE3\t10\n\
             cycle_sample\tGENE_A\tE3\tE1\t10\n\
             linear_sample\tGENE_A\tE1\tE2\t10\n\
             linear_sample\tGENE_A\tE2\tE3\t8\n",
        )
        .unwrap();
        let results = analyze_junctions(rows).unwrap();
        let cycle = results
            .iter()
            .find(|result| result.sample_id == "cycle_sample")
            .unwrap();
        let linear = results
            .iter()
            .find(|result| result.sample_id == "linear_sample")
            .unwrap();

        let cycle_svg = render_svg(cycle);
        let linear_svg = render_svg(linear);

        assert!(cycle_svg.contains("Sample: cycle_sample"));
        assert!(cycle_svg.contains("Cycle fraction: 1.0000"));
        assert!(linear_svg.contains("Sample: linear_sample"));
        assert!(linear_svg.contains("Cycle fraction: 0.0000"));
        assert_ne!(cycle_svg, linear_svg);
    }
}
