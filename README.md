# Historical drug–target discovery database backfill audit

Code and frozen data for **Historical database backfill distorts prospective
drug–target discovery benchmarks**. Project: biomedical_ai_historical_discovery.

The audit compares recovered ChEMBL releases with document-year rollbacks at
2018, 2020 and 2022 cutoffs. Database-membership mismatch and future-label
erasure are distinct from fixed-test model-performance effects. The 2022
follow-up is right-censored at the operational availability proxy 2025-09-19.

## Quick start

Use Python 3.11 or later and install `requirements.txt`. Original recorded
versions are in ENVIRONMENT; the portable presentation route requires NumPy,
pandas, matplotlib and openpyxl. From the repository root:

```sh
python CODE/verify_package.py
python CODE/reproduce_figures.py
python CODE/reproduce_v5.py
python CODE/export_source_data.py
```

Outputs are written only to REPRODUCED. The first plotting command preserves
the ten supplemental panels; the second generates the current five main
figures and S11 from frozen tables. Its Figure 3c uses MRR; top-tie any-hit is
S11. SVG/PNG output is used for portable reproduction. Fonts can change visual
appearance. Source workbook cell values, rather than ZIP bytes, define equality.

## Files and data access

- CODE: preserved scientific implementations and portable presentation/export.
- 10_phase1, 11_phase2, 14_repair: compact frozen summaries and definitions.
- SOURCE_DATA: current figure workbook, panel manifest, frozen payload and
  supplementary tables.
- FIGURES: frozen vector exports, including S1–S11.
- MANIFESTS: publisher release objects, original input provenance, public-file
  hashes and large-data archive membership/checksums.
- CONFIG, configs, ENVIRONMENT: configuration and recorded software versions.

Large candidates, scores, release-edge products, structures, descriptors and
chemical/graph audits are provided as release assets, with a linked Zenodo
dataset. Download the five data ZIPs from the v1.0.0 release and extract them
into a separate reconstruction copy of this repository, preserving paths.
Verify each archive and extracted file against MANIFESTS/DATA_ARCHIVES.tsv and
MANIFESTS/DATA_FILES.tsv. See DATA_DICTIONARY.md for keys, units and scope.

The bundled derived-table route reproduces the current presentation. It does
not independently recompute every metric from candidate scores. Scientific
implementations are supplied for inspection; a fully audited raw-SQLite-to-all-
results driver is not supplied. Historical Phase 1/2 labels can differ from the
authoritative 14_repair products and are retained only for reconciliation.
Do not replace authoritative labels with legacy labels when reusing scores.
Heavy scripts may train models or overwrite outputs: use a reconstruction copy,
HISTORICAL_DTI_ROOT and HISTORICAL_DTI_CONFIG, never a frozen original checkout.

## Interpretation

Figure 4a counts candidate compound–target pairs, with a historical ligand as
the classification witness. U2 and broader historical-active denominators
differ. Expected within-compound rank-1 placement, Recall@1, MRR and top-tie
coverage are distinct statistics. Chemical thresholds are post hoc; the
confidence-rule fallback inspected outcomes; model-panel commitment timing
is unknown. Fixed-seed graph results are descriptive. No LLM/Open Targets
benchmark was run. Missing raw archive checksums are marked UNKNOWN, not
invented. A recovered first appearance does not establish biological priority.

## Licence and citation

Code: MIT. ChEMBL-derived data: CC BY-SA 3.0; see DATA_LICENSE.md. Zenodo
creator metadata use the project name. The research manuscript has a separate
human author list. Cite the versioned Zenodo records linked from the release.

## Release verification

The portable main-figure route was executed against the packaged tables: all
five main PNGs and S11 were pixel-identical to the frozen exports. All 5,129
payload cells matched the supplied source workbook. These checks validate
presentation/export, not an independent recalculation of model metrics.
