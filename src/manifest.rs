//! Parses and validates the sample manifest described in
//! `docs/data_contract.md`: a TSV assigning each `sample_id` used in the
//! junction table to a case or control group. Required columns are found
//! by name rather than a fixed header, since real metadata files carry many
//! extra columns this tool does not need.
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

/// Case/control group for one sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Condition {
    Case,
    Control,
}

impl std::fmt::Display for Condition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match self {
            Condition::Case => "case",
            Condition::Control => "control",
        };
        write!(f, "{text}")
    }
}

/// One row of the sample manifest: which group `sample_id` belongs to.
#[derive(Clone, Debug)]
pub struct SampleRecord {
    pub sample_id: String,
    pub condition: Condition,
}

/// Reads `path` from disk and parses it as a sample manifest.
pub fn parse_manifest(path: &Path) -> Result<Vec<SampleRecord>, String> {
    let contents = fs::read_to_string(path).map_err(|error| {
        format!(
            "could not read sample manifest {}: {error}; check that the path exists, is readable, and points to a UTF-8 text file",
            path.display()
        )
    })?;
    parse_manifest_text(&contents)
        .map_err(|error| format!("invalid sample manifest {}: {error}", path.display()))
}

/// Parses already-read TSV text into sample records. Split out from
/// `parse_manifest` so tests can exercise parsing without touching the
/// filesystem.
pub fn parse_manifest_text(contents: &str) -> Result<Vec<SampleRecord>, String> {
    let mut lines = contents.lines();
    let header = lines.next().ok_or(
        "file is empty; add a tab-separated header containing 'sample_id' and 'condition'",
    )?;
    let columns: Vec<&str> = header.split('\t').collect();

    // Found by name, not a fixed position — see docs/data_contract.md.
    let sample_id_index = find_required_column(&columns, "sample_id")?;
    let condition_index = find_required_column(&columns, "condition")?;

    let mut records = Vec::new();
    let mut seen_sample_ids: BTreeSet<String> = BTreeSet::new();

    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2; // the header occupies line 1
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != columns.len() {
            return Err(format!(
                "line {line_number}: expected {} tab-separated fields to match the header, found {}; check for missing values or extra tabs",
                columns.len(),
                fields.len()
            ));
        }

        let record = parse_record(&fields, sample_id_index, condition_index, line_number)?;

        if !seen_sample_ids.insert(record.sample_id.clone()) {
            return Err(format!(
                "line {line_number}, column 'sample_id': duplicate sample_id '{}'; keep exactly one manifest row per sample",
                record.sample_id
            ));
        }

        records.push(record);
    }

    if records.is_empty() {
        return Err(
            "file contains a header but no sample rows; add at least one sample with condition 'case' or 'control'"
                .to_owned(),
        );
    }

    Ok(records)
}

/// Finds `name`'s position in the header by exact name match. Must error if
/// the column is missing, and must error if it appears more than once — an
/// ambiguous header should never be silently resolved by picking whichever
/// match comes first.
fn find_required_column(columns: &[&str], name: &str) -> Result<usize, String> {
    let mut matches = columns
        .iter()
        .enumerate()
        .filter(|&(_, &column)| column == name)
        .map(|(index, _)| index);

    let Some(index) = matches.next() else {
        return Err(format!(
            "header is missing required column '{name}'; add one tab-separated '{name}' column"
        ));
    };
    if matches.next().is_some() {
        return Err(format!(
            "header has more than one '{name}' column; keep exactly one"
        ));
    }
    Ok(index)
}

/// Validates one data row and builds its `SampleRecord`.
fn parse_record(
    fields: &[&str],
    sample_id_index: usize,
    condition_index: usize,
    line_number: usize,
) -> Result<SampleRecord, String> {
    let sample_id = fields[sample_id_index].trim();
    if sample_id.is_empty() {
        return Err(format!(
            "line {line_number}, column 'sample_id': value cannot be empty; enter the identifier used in the junction TSV"
        ));
    }

    // Condition values are matched case-insensitively (unlike header column
    // names) because they repeat once per row and are far more likely to be
    // typed inconsistently by hand; see docs/data_contract.md.
    let condition_text = fields[condition_index].trim().to_ascii_lowercase();
    let condition = match condition_text.as_str() {
        "case" => Condition::Case,
        "control" => Condition::Control,
        other => {
            return Err(format!(
                "line {line_number}, column 'condition': condition must be 'case' or 'control', found '{other}'; replace it with one of the accepted labels"
            ));
        }
    };

    Ok(SampleRecord {
        sample_id: sample_id.to_owned(),
        condition,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_valid_manifest_with_extra_columns() {
        let input = "sample_id\tcondition\ttissue\n\
                     s1\tcase\tendometrium\n\
                     s2\tControl\tendometrium\n";
        let records = parse_manifest_text(input).expect("valid manifest should parse");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sample_id, "s1");
        assert_eq!(records[0].condition, Condition::Case);
        assert_eq!(records[1].sample_id, "s2");
        assert_eq!(records[1].condition, Condition::Control);
    }

    #[test]
    fn condition_is_matched_case_insensitively() {
        let input = "sample_id\tcondition\n\
                     s1\tCASE\n\
                     s2\tcOntrol\n";
        let records = parse_manifest_text(input).expect("mixed-case condition should parse");
        assert_eq!(records[0].condition, Condition::Case);
        assert_eq!(records[1].condition, Condition::Control);
    }

    #[test]
    fn rejects_invalid_condition_values() {
        let input = "sample_id\tcondition\n\
                     s1\tendometriosis\n";
        let error = parse_manifest_text(input).unwrap_err();
        assert!(error.contains("condition must be"));
        assert!(error.contains("endometriosis"));
    }

    #[test]
    fn rejects_duplicate_sample_ids() {
        let input = "sample_id\tcondition\n\
                     s1\tcase\n\
                     s1\tcontrol\n";
        let error = parse_manifest_text(input).unwrap_err();
        assert!(error.contains("duplicate sample_id"));
    }

    #[test]
    fn rejects_a_header_missing_a_required_column() {
        let input = "sample_id\ttissue\n\
                     s1\tendometrium\n";
        let error = parse_manifest_text(input).unwrap_err();
        assert!(error.contains("missing required column 'condition'"));
    }

    #[test]
    fn rejects_a_header_with_a_duplicated_required_column() {
        let input = "sample_id\tcondition\tcondition\n\
                     s1\tcase\tcase\n";
        let error = parse_manifest_text(input).unwrap_err();
        assert!(error.contains("more than one"));
    }

    #[test]
    fn rejects_a_row_with_the_wrong_number_of_fields() {
        let input = "sample_id\tcondition\ttissue\n\
                     s1\tcase\n";
        let error = parse_manifest_text(input).unwrap_err();
        assert!(error.contains("expected 3"));
    }
}
