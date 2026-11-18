//! Carry supplied candidate evidence into results by exact versioned gene ID.
//! These labels are input annotations, never inferred from Hodge measurements.
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

pub const COLUMNS: [&str; 5] = [
    "gene_symbol",
    "is_risk",
    "is_environment",
    "is_overlap",
    "category",
];

#[derive(Debug, Default)]
pub struct GeneAnnotation {
    pub gene_symbol: String,
    pub is_risk: String,
    pub is_environment: String,
    pub is_overlap: String,
    pub category: String,
}

impl GeneAnnotation {
    pub fn fields(&self) -> [&str; 5] {
        [
            &self.gene_symbol,
            &self.is_risk,
            &self.is_environment,
            &self.is_overlap,
            &self.category,
        ]
    }
}

pub fn read_annotations(path: &Path) -> Result<BTreeMap<String, GeneAnnotation>, String> {
    let text = fs::read_to_string(path).map_err(|error| {
        format!(
            "could not read gene annotations {}: {error}",
            path.display()
        )
    })?;
    parse_annotations(&text)
        .map_err(|error| format!("invalid gene annotations {}: {error}", path.display()))
}

fn parse_annotations(text: &str) -> Result<BTreeMap<String, GeneAnnotation>, String> {
    crate::statistics::parse_gene_family_text(text)?;
    let mut lines = text.lines();
    let header: Vec<_> = lines.next().unwrap().split('\t').collect();
    for name in COLUMNS {
        if header.iter().filter(|&&value| value == name).count() > 1 {
            return Err(format!(
                "duplicate annotation column '{name}'; keep one unambiguous value per gene"
            ));
        }
    }
    let gene_index = header.iter().position(|&name| name == "gene_id").unwrap();
    let indices = COLUMNS.map(|name| header.iter().position(|&value| value == name));
    let mut annotations = BTreeMap::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        let values: Vec<_> = line.split('\t').collect();
        let [gene_symbol, is_risk, is_environment, is_overlap, category] =
            indices.map(|index| index.map_or("", |index| values[index]).to_owned());
        annotations.insert(
            values[gene_index].to_owned(),
            GeneAnnotation {
                gene_symbol,
                is_risk,
                is_environment,
                is_overlap,
                category,
            },
        );
    }
    Ok(annotations)
}

#[cfg(test)]
mod tests {
    use super::parse_annotations;

    #[test]
    fn annotations_join_by_id_with_no_category_inference() {
        let rows = parse_annotations("category\tgene_id\tgene_symbol\tis_risk\tis_environment\tis_overlap\noverlap\tG2.3\tB\t1\t1\t1\nrisk_only\tG1.1\tA\t1\t0\t0\n").unwrap();
        assert_eq!(rows["G2.3"].fields(), ["B", "1", "1", "1", "overlap"]);
        assert!(!rows.contains_key("G2"));
        let plain = parse_annotations("gene_id\nG\n").unwrap();
        assert_eq!(plain["G"].fields(), [""; 5]);
    }

    #[test]
    fn ambiguous_annotation_headers_are_rejected() {
        assert!(
            parse_annotations("gene_id\tcategory\tcategory\nG\trisk_only\toverlap\n")
                .unwrap_err()
                .contains("duplicate annotation column")
        );
    }
}
