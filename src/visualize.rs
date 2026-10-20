//! Dependency-free SVG rendering for one analyzed sample/gene graph.
//!
//! The renderer deliberately uses the same `GroupResult` produced by the
//! validated analysis path. It does not recalculate or reinterpret the Hodge
//! decomposition; it only turns the graph, support, and cycle component into
//! a browser-readable figure.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

use crate::GroupResult;

const MIN_WIDTH: f64 = 980.0;
const NODE_WIDTH: f64 = 116.0;
const NODE_HEIGHT: f64 = 54.0;
const NODE_Y: f64 = 420.0;
const NODE_SPACING: f64 = 170.0;
const SIDE_MARGIN: f64 = 120.0;

/// Render one analyzed splice graph as a standalone SVG document.
pub fn render_svg(result: &GroupResult) -> String {
    let graph = &result.graph;
    let decomposition = &result.decomposition;
    let width = MIN_WIDTH
        .max(SIDE_MARGIN * 2.0 + NODE_SPACING * graph.nodes.len().saturating_sub(1) as f64);
    let first_table_row_y = 570.0;
    let table_row_height = 29.0;
    let table_end_y = first_table_row_y + table_row_height * graph.edges.len() as f64;
    let footer_line_y = table_end_y + 27.0;
    let canvas_height = footer_line_y + 88.0;
    let max_count = graph
        .edges
        .iter()
        .map(|edge| edge.count)
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let max_cycle = decomposition
        .cycle
        .iter()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);

    let mut svg = String::new();
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{canvas_height:.0}" viewBox="0 0 {width:.0} {canvas_height:.0}" role="img" aria-labelledby="title description">
<title id="title">Splice graph for sample {} and gene {}</title>
<desc id="description">A directed splice graph. Arrow thickness represents junction support and arrow color represents the signed Hodge cycle component.</desc>
<defs>
  <linearGradient id="background" x1="0" y1="0" x2="1" y2="1">
    <stop offset="0%" stop-color="#fff9fb"/>
    <stop offset="100%" stop-color="#f7e7ef"/>
  </linearGradient>
  <filter id="shadow" x="-20%" y="-20%" width="140%" height="140%">
    <feDropShadow dx="0" dy="3" stdDeviation="4" flood-color="#4d203b" flood-opacity="0.16"/>
  </filter>
  <marker id="arrow-positive" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#b12a76"/>
  </marker>
  <marker id="arrow-negative" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#167b86"/>
  </marker>
  <marker id="arrow-neutral" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#786b75"/>
  </marker>
  <style>
    .sans {{ font-family: ui-sans-serif, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
    .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; }}
  </style>
</defs>
<rect width="100%" height="100%" rx="22" fill="url(#background)"/>
<rect x="1.5" y="1.5" width="{border_width:.0}" height="{border_height:.0}" rx="21" fill="none" stroke="#c997b2" stroke-width="3"/>

<text class="sans" x="52" y="55" font-size="27" font-weight="750" fill="#65334f">SPLICE GIRL · GRAPH VIEW</text>
<text class="sans" x="52" y="91" font-size="18" fill="#73465e">Sample: <tspan font-weight="700">{}</tspan>   ·   Gene: <tspan font-weight="700">{}</tspan></text>
"##,
        xml_escape(&result.sample_id),
        xml_escape(&result.gene_id),
        xml_escape(&result.sample_id),
        xml_escape(&result.gene_id),
        border_width = width - 3.0,
        border_height = canvas_height - 3.0,
    ));

    let status = if decomposition.checks.passed {
        "PASS"
    } else {
        "FAIL"
    };
    svg.push_str(&summary_card(
        52.0,
        "CYCLE FRACTION",
        &format!("{:.1}%", decomposition.cycle_fraction * 100.0),
    ));
    svg.push_str(&summary_card(
        250.0,
        "NODES",
        &graph.nodes.len().to_string(),
    ));
    svg.push_str(&summary_card(
        448.0,
        "JUNCTIONS",
        &graph.edges.len().to_string(),
    ));
    svg.push_str(&summary_card(646.0, "MATH CHECKS", status));

    svg.push_str(
        r##"<text class="sans" x="52" y="243" font-size="14" font-weight="700" fill="#65334f">HOW TO READ THE GRAPH</text>
<line x1="52" y1="266" x2="102" y2="266" stroke="#786b75" stroke-width="7" stroke-linecap="round"/>
<text class="sans" x="113" y="271" font-size="14" fill="#4d3543">thicker arrow = more supporting reads</text>
<line x1="395" y1="266" x2="445" y2="266" stroke="#b12a76" stroke-width="5" stroke-linecap="round"/>
<text class="sans" x="456" y="271" font-size="14" fill="#4d3543">positive cycle contribution</text>
<line x1="695" y1="266" x2="745" y2="266" stroke="#167b86" stroke-width="5" stroke-linecap="round"/>
<text class="sans" x="756" y="271" font-size="14" fill="#4d3543">negative cycle contribution</text>
<text class="sans" x="52" y="310" font-size="18" font-weight="700" fill="#65334f">Observed splice junctions</text>
"##,
    );

    // Draw edges first so node boxes remain legible on top of them.
    for (index, edge) in graph.edges.iter().enumerate() {
        let from_x = node_x(edge.from);
        let to_x = node_x(edge.to);
        let direction = if to_x >= from_x { 1.0 } else { -1.0 };
        let start_x = from_x + direction * NODE_WIDTH / 2.0;
        let end_x = to_x - direction * NODE_WIDTH / 2.0;
        let midpoint_x = (start_x + end_x) / 2.0;
        let span = edge.from.abs_diff(edge.to).max(1) as f64;
        let arc_height = (68.0 + 28.0 * span).min(175.0);
        let control_y = if direction > 0.0 {
            NODE_Y - arc_height
        } else {
            NODE_Y + arc_height
        };
        let label_y = (NODE_Y + control_y) / 2.0 + if direction > 0.0 { -8.0 } else { 18.0 };
        let cycle = decomposition.cycle[index];
        let (color, marker) = cycle_style(cycle, max_cycle);
        let stroke_width = 2.5 + 5.5 * (edge.count / max_count).sqrt();
        let junction_label = format!("J{}", index + 1);

        svg.push_str(&format!(
            r##"<path d="M {start_x:.1} {NODE_Y:.1} Q {midpoint_x:.1} {control_y:.1} {end_x:.1} {NODE_Y:.1}" fill="none" stroke="{color}" stroke-width="{stroke_width:.2}" stroke-linecap="round" opacity="0.88" marker-end="url(#{marker})"/>
<rect x="{badge_x:.1}" y="{badge_y:.1}" width="32" height="22" rx="8" fill="#fffdfd" stroke="{color}" stroke-width="1.5"/>
<text class="mono" x="{midpoint_x:.1}" y="{badge_text_y:.1}" text-anchor="middle" font-size="11" font-weight="750" fill="{color}">{junction_label}</text>
"##,
            badge_x = midpoint_x - 16.0,
            badge_y = label_y - 15.0,
            badge_text_y = label_y,
            junction_label = xml_escape(&junction_label),
        ));
    }

    for (index, node) in graph.nodes.iter().enumerate() {
        let x = node_x(index);
        let display = shorten(node, 15);
        svg.push_str(&format!(
            r##"<g filter="url(#shadow)">
  <rect x="{left:.1}" y="{top:.1}" width="{NODE_WIDTH:.1}" height="{NODE_HEIGHT:.1}" rx="14" fill="#fffdfd" stroke="#8d5574" stroke-width="2.5"/>
</g>
<text class="mono" x="{x:.1}" y="{text_y:.1}" text-anchor="middle" font-size="14" font-weight="700" fill="#552941">{display}<title>{full}</title></text>
"##,
            left = x - NODE_WIDTH / 2.0,
            top = NODE_Y - NODE_HEIGHT / 2.0,
            text_y = NODE_Y + 5.0,
            display = xml_escape(&display),
            full = xml_escape(node),
        ));
    }

    svg.push_str(
        r##"<text class="sans" x="52" y="520" font-size="18" font-weight="700" fill="#65334f">Junction details</text>
<text class="sans" x="66" y="550" font-size="12" font-weight="750" letter-spacing="0.8" fill="#8a5b74">ID</text>
<text class="sans" x="125" y="550" font-size="12" font-weight="750" letter-spacing="0.8" fill="#8a5b74">FROM → TO</text>
<text class="sans" x="390" y="550" font-size="12" font-weight="750" letter-spacing="0.8" fill="#8a5b74">READ SUPPORT</text>
<text class="sans" x="565" y="550" font-size="12" font-weight="750" letter-spacing="0.8" fill="#8a5b74">NORMALIZED SIGNAL</text>
<text class="sans" x="785" y="550" font-size="12" font-weight="750" letter-spacing="0.8" fill="#8a5b74">CYCLE COMPONENT</text>
"##,
    );

    for (index, edge) in graph.edges.iter().enumerate() {
        let row_y = first_table_row_y + table_row_height * index as f64;
        let row_fill = if index % 2 == 0 { "#fffafd" } else { "#f8edf3" };
        let cycle = decomposition.cycle[index];
        let (cycle_color, _) = cycle_style(cycle, max_cycle);
        let from = shorten(&graph.nodes[edge.from], 24);
        let to = shorten(&graph.nodes[edge.to], 24);
        let connection = format!("{from} → {to}");
        svg.push_str(&format!(
            r##"<rect x="52" y="{row_top:.1}" width="876" height="25" rx="6" fill="{row_fill}"/>
<text class="mono" x="66" y="{row_y:.1}" font-size="12" font-weight="750" fill="{cycle_color}">J{number}</text>
<text class="mono" x="125" y="{row_y:.1}" font-size="12" fill="#4d3543">{connection}</text>
<text class="mono" x="390" y="{row_y:.1}" font-size="12" fill="#4d3543">{support}</text>
<text class="mono" x="565" y="{row_y:.1}" font-size="12" fill="#4d3543">{signal:.4}</text>
<text class="mono" x="785" y="{row_y:.1}" font-size="12" font-weight="700" fill="{cycle_color}">{cycle:+.4}</text>
"##,
            row_top = row_y - 18.0,
            number = index + 1,
            connection = xml_escape(&connection),
            support = xml_escape(&format_count(edge.count)),
            signal = graph.signal[index],
        ));
    }

    svg.push_str(&format!(
        r##"<line x1="52" y1="{footer_line_y:.1}" x2="{footer_end:.1}" y2="{footer_line_y:.1}" stroke="#d7afc3" stroke-width="1.5"/>
<text class="sans" x="52" y="{footer_text_y:.1}" font-size="15" fill="#5f4352">Cycle fraction = cycle energy ÷ total signal energy. Edge signs follow the stored arrow direction.</text>
<text class="sans" x="52" y="{caveat_y:.1}" font-size="14" fill="#765e6b">A cycle-space signal can reveal alternative-path structure; by itself, it does not prove circular RNA or disease causation.</text>
</svg>
"##,
        footer_end = width - 52.0,
        footer_text_y = footer_line_y + 34.0,
        caveat_y = footer_line_y + 61.0,
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

fn cycle_style(value: f64, max_abs: f64) -> (&'static str, &'static str) {
    let near_zero = max_abs <= f64::EPSILON || value.abs() <= max_abs * 1e-8;
    if near_zero {
        ("#786b75", "arrow-neutral")
    } else if value > 0.0 {
        ("#b12a76", "arrow-positive")
    } else {
        ("#167b86", "arrow-negative")
    }
}

fn summary_card(x: f64, label: &str, value: &str) -> String {
    format!(
        r##"<g filter="url(#shadow)">
  <rect x="{x:.1}" y="116" width="172" height="82" rx="14" fill="#fffdfd" stroke="#dab3c7"/>
</g>
<text class="sans" x="{label_x:.1}" y="143" font-size="12" font-weight="750" letter-spacing="1.1" fill="#8a5b74">{label}</text>
<text class="sans" x="{label_x:.1}" y="180" font-size="26" font-weight="750" fill="#552941">{value}</text>
"##,
        label_x = x + 17.0,
        label = xml_escape(label),
        value = xml_escape(value),
    )
}

fn format_count(count: f64) -> String {
    if count.fract().abs() <= f64::EPSILON {
        format!("{count:.0}")
    } else {
        format!("{count:.2}")
    }
}

fn shorten(value: &str, max_chars: usize) -> String {
    let count = value.chars().count();
    if count <= max_chars {
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
    fn svg_contains_graph_metrics_and_edge_encoding() {
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
        assert!(svg.contains("Sample: <tspan font-weight=\"700\">case_01</tspan>"));
        assert!(svg.contains("GENE_DEMO"));
        assert!(svg.contains("CYCLE FRACTION"));
        assert!(svg.contains("Junction details"));
        assert!(svg.contains("E1 → E2"));
        assert!(svg.contains(">100</text>"));
        assert!(svg.contains("marker-end=\"url(#arrow-positive)\""));
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

        assert!(svg.contains("sample&amp;amp"));
        assert!(svg.contains("gene&lt;1&gt;"));
        assert!(svg.contains("<title>a-very-long-node-name</title>"));
        assert!(!svg.contains("gene<1>"));
    }
}
