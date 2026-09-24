//! Reads and validates the tab-separated junction input described in the
//! README: one row per observed splice junction, with columns
//! `sample_id, gene_id, from, to, count`. This is the only place raw input
//! is trusted from disk; everything downstream assumes rows already passed
//! these checks.
use std::fs;
use std::path::Path;

/// One observed splice junction: `count` reads support the `from -> to`
/// edge for a given sample and gene.
#[derive(Clone, Debug)]
pub struct Junction {
    pub sample_id: String,
    pub gene_id: String,
    pub from: String,
    pub to: String,
    pub count: f64,
}

const REQUIRED_HEADER: [&str; 5] = ["sample_id", "gene_id", "from", "to", "count"];

/// Reads `path` from disk and parses it as junction TSV.
pub fn parse_tsv(path: &Path) -> Result<Vec<Junction>, String> {
    let contents = fs::read_to_string(path)
        .map_err(|error| format!("could not read {}: {error}", path.display()))?;
    parse_tsv_text(&contents)
}

/// Parses already-read TSV text into junction rows. Split out from
/// `parse_tsv` so tests can exercise parsing without touching the filesystem.
pub fn parse_tsv_text(contents: &str) -> Result<Vec<Junction>, String> {
    let mut lines = contents.lines();
    let header = lines.next().ok_or("input is empty")?;
    let columns: Vec<&str> = header.split('\t').collect();
    if columns != REQUIRED_HEADER {
        return Err(format!(
            "expected header {}, found {}",
            REQUIRED_HEADER.join("\\t"),
            columns.join("\\t")
        ));
    }

    let mut rows = Vec::new();
    for (offset, line) in lines.enumerate() {
        // +2: line 1 is the header, and `enumerate` is zero-based.
        let line_number = offset + 2;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != REQUIRED_HEADER.len() {
            return Err(format!(
                "line {line_number}: expected 5 tab-separated fields, found {}",
                fields.len()
            ));
        }

        // Only the first four columns are identifiers; `count` gets its own
        // numeric check below.
        for (column, value) in REQUIRED_HEADER[..4].iter().zip(&fields[..4]) {
            if value.trim().is_empty() {
                return Err(format!("line {line_number}: {column} cannot be empty"));
            }
        }

        let count: f64 = fields[4]
            .parse()
            .map_err(|_| format!("line {line_number}: count must be a number"))?;
        if !count.is_finite() || count < 0.0 {
            return Err(format!(
                "line {line_number}: count must be finite and nonnegative"
            ));
        }
        // Self-loops have no meaningful incidence-matrix column (a -1 and a
        // +1 would cancel at the same node), so the graph math can't handle
        // them yet; reject explicitly rather than silently dropping signal.
        if fields[2] == fields[3] {
            return Err(format!(
                "line {line_number}: self-loop {} -> {} is not supported in the basic prototype",
                fields[2], fields[3]
            ));
        }

        rows.push(Junction {
            sample_id: fields[0].to_owned(),
            gene_id: fields[1].to_owned(),
            from: fields[2].to_owned(),
            to: fields[3].to_owned(),
            count,
        });
    }

    if rows.is_empty() {
        return Err("input contains no junction rows".to_owned());
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::parse_tsv_text;

    #[test]
    fn rejects_negative_counts() {
        let input = "sample_id\tgene_id\tfrom\tto\tcount\ns1\tg1\tE1\tE2\t-1\n";
        assert!(parse_tsv_text(input).unwrap_err().contains("nonnegative"));
    }

    #[test]
    fn rejects_empty_sample_ids() {
        let input = "sample_id\tgene_id\tfrom\tto\tcount\n\tg1\tE1\tE2\t1\n";
        assert!(parse_tsv_text(input).unwrap_err().contains("sample_id"));
    }

    #[test]
    fn rejects_rows_with_the_wrong_number_of_fields() {
        let input = "sample_id\tgene_id\tfrom\tto\tcount\ns1\tg1\tE1\tE2\n";
        assert!(parse_tsv_text(input).unwrap_err().contains("expected 5"));
    }

    #[test]
    fn rejects_non_finite_counts() {
        let input = "sample_id\tgene_id\tfrom\tto\tcount\ns1\tg1\tE1\tE2\tNaN\n";
        assert!(parse_tsv_text(input).unwrap_err().contains("nonnegative"));
    }

    #[test]
    fn rejects_self_loops() {
        let input = "sample_id\tgene_id\tfrom\tto\tcount\ns1\tg1\tE1\tE1\t1\n";
        assert!(parse_tsv_text(input).unwrap_err().contains("self-loop"));
    }
}
