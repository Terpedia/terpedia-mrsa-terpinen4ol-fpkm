"""Descriptive audit of the archived GSE157638 processed FPKM matrix."""

from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import statistics
from itertools import combinations


ROOT = Path("research/mrsa-terpinen4ol-fpkm")
MATRIX = ROOT / "data/GSE157638_fpkm.txt.gz"
SERIES = ROOT / "data/GSE157638_series_matrix.txt.gz"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


with gzip.open(SERIES, "rt") as stream:
    series_lines = stream.read().splitlines()


def sample_field(key):
    line = next(line for line in series_lines if line.startswith(key + "\t"))
    return next(csv.reader([line], delimiter="\t"))[1:]


titles = sample_field("!Sample_title")
accessions = sample_field("!Sample_geo_accession")
conditions = sample_field("!Sample_source_name_ch1")
samples = [
    {"sample_id": title, "geo_accession": geo, "condition": condition}
    for title, geo, condition in zip(titles, accessions, conditions, strict=True)
]
assert len(samples) == 6 and len(set(titles)) == 6
controls = [item["sample_id"] for item in samples if item["condition"] == "Control"]
treated = [
    item["sample_id"]
    for item in samples
    if item["condition"] == "Treatment with terpinen-4-ol"
]
assert len(controls) == len(treated) == 3

with gzip.open(MATRIX, "rt") as stream:
    reader = csv.DictReader(stream, delimiter="\t")
    assert reader.fieldnames == ["geneID", *titles], reader.fieldnames
    matrix = {}
    for row in reader:
        gene = row["geneID"].strip()
        assert gene and gene not in matrix, "Missing or duplicate feature ID"
        values = {sample: float(row[sample]) for sample in titles}
        assert all(math.isfinite(value) and value >= 0 for value in values.values())
        matrix[gene] = values
assert matrix

sample_qc = []
for item in samples:
    sample = item["sample_id"]
    values = [row[sample] for row in matrix.values()]
    sample_qc.append(
        {
            **item,
            "features_nonzero": sum(value > 0 for value in values),
            "median_fpkm": statistics.median(values),
            "sum_fpkm": sum(values),
        }
    )
with (RESULTS / "sample_qc.tsv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=sample_qc[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(sample_qc)

pairs = []
for left, right in combinations(titles, 2):
    correlation = statistics.correlation(
        [math.log1p(row[left]) for row in matrix.values()],
        [math.log1p(row[right]) for row in matrix.values()],
    )
    pair_type = (
        "within_control"
        if left in controls and right in controls
        else "within_treated"
        if left in treated and right in treated
        else "between_arms"
    )
    pairs.append(
        {
            "left": left,
            "right": right,
            "pair_type": pair_type,
            "pearson_log1p_fpkm": correlation,
        }
    )
correlation_medians = {
    kind: statistics.median(
        row["pearson_log1p_fpkm"] for row in pairs if row["pair_type"] == kind
    )
    for kind in ("within_control", "within_treated", "between_arms")
}
with (RESULTS / "sample_correlations.tsv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=pairs[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(pairs)

effects = []
for gene, row in matrix.items():
    control_mean = statistics.mean(row[sample] for sample in controls)
    treated_mean = statistics.mean(row[sample] for sample in treated)
    if control_mean >= 1:
        effects.append(
            {
                "gene_id": gene,
                "mean_control_fpkm": control_mean,
                "mean_treated_fpkm": treated_mean,
                "log2_ratio_pc_0_1": math.log2((treated_mean + 0.1) / (control_mean + 0.1)),
                "log2_ratio_pc_1": math.log2((treated_mean + 1) / (control_mean + 1)),
            }
        )
assert effects
threshold_counts = {
    str(pc): {
        "at_least_twofold_higher": sum(row[f"log2_ratio_pc_{pc}"] >= 1 for row in effects),
        "at_least_twofold_lower": sum(row[f"log2_ratio_pc_{pc}"] <= -1 for row in effects),
    }
    for pc in ("0_1", "1")
}
with (RESULTS / "descriptive_effects.tsv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=effects[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(effects, key=lambda row: row["gene_id"]))

summary = {
    "study_accession": "GSE157638",
    "assay": "author-provided processed FPKM",
    "matrix_sha256": sha256(MATRIX),
    "geo_series_matrix_sha256": sha256(SERIES),
    "sample_count": len(samples),
    "control_count": len(controls),
    "treated_count": len(treated),
    "matrix_feature_count": len(matrix),
    "filtered_feature_rule": "mean control FPKM >= 1",
    "filtered_feature_count": len(effects),
    "median_pairwise_log1p_fpkm_correlations": correlation_medians,
    "descriptive_twofold_counts": threshold_counts,
    "limitations": [
        "FPKM is not an integer-count differential-expression input; no new p values or false-discovery estimates are calculated.",
        "The gene identifiers are not independently mapped to orthologs or pathways.",
        "Global FPKM distributions and pairwise correlations may reflect technical as well as biological differences.",
        "This is one study; cross-species replication and direct target dependence are not tested.",
    ],
}
with (RESULTS / "analysis_summary.json").open("w") as stream:
    json.dump(summary, stream, indent=2)
    stream.write("\n")
print(json.dumps(summary, indent=2))
