# Terpedia MRSA terpinen-4-ol transcriptomics audit

This forkable project audits the [GSE157638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157638) processed FPKM matrix for three untreated and three terpinen-4-ol exposed *Staphylococcus aureus* ATCC43300 libraries. The [executed Jupyter notebook](notebooks/terputer-mrsa-terpinen4ol-fpkm.ipynb), [claim-linked report](research/mrsa-terpinen4ol-fpkm/generated/white-paper.md), [claim ledger](research/mrsa-terpinen4ol-fpkm/generated/claim-ledger.tsv), and [execution receipt](research/execution-receipts/terputer-mrsa-terpinen4ol-fpkm.json) are committed.

The notebook parses 2,925 unique features and records sample QC, pairwise log1p(FPKM) correlations, and descriptive group-mean ratios for 2,456 control-expressed features. These ratios are **not** new differential-expression calls. The [source article](https://link.springer.com/article/10.1186/s12866-021-02348-2) reports 304 DEGs and a purine/pyrimidine interpretation; this project does not reproduce its count model, pathway analysis, or metabolomics.

## Reproduce

Use Python 3.10 or newer and Node.js:

```sh
python -m pip install -r requirements.txt
python restore-inputs.py research/mrsa-terpinen4ol-fpkm/report-spec.json
python notebook-runner.py notebooks/terputer-mrsa-terpinen4ol-fpkm.ipynb . 120 '["research/mrsa-terpinen4ol-fpkm/results/analysis_summary.json","research/mrsa-terpinen4ol-fpkm/results/sample_qc.tsv","research/mrsa-terpinen4ol-fpkm/results/sample_correlations.tsv","research/mrsa-terpinen4ol-fpkm/results/descriptive_effects.tsv"]' '["research/mrsa-terpinen4ol-fpkm/analysis.py","research/mrsa-terpinen4ol-fpkm/data/GSE157638_fpkm.txt.gz","research/mrsa-terpinen4ol-fpkm/data/GSE157638_series_matrix.txt.gz"]'
node report-generator.mjs research/mrsa-terpinen4ol-fpkm/report-spec.json
```

The two GEO downloads are checksum pinned and excluded from Git. `restore-inputs.py` verifies them before notebook execution. To attach this project to a `terpedia-report` workspace, fork it under the same repository name and use `add_research_project_submodule` with slug `mrsa-terpinen4ol-fpkm`; the submodule path is `projects/mrsa-terpinen4ol-fpkm`.

See the [project methods and interpretation](research/mrsa-terpinen4ol-fpkm/README.md).
