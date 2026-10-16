# Data dictionary

Files are UTF-8 TSV/JSON, Parquet, SVG and XLSX. Names and columns are preserved
from frozen scientific outputs. Empty or NaN metrics denote unavailable values,
not zero. TRUE/FALSE or 1/0 are boolean flags. Counts are numbers of edges,
candidates or positive compounds as named. Fractions range from 0 to 1;
percentage displays multiply fractions by 100. Enrichment and score/rank
statistics are dimensionless. Release dates are ISO dates and operational
availability proxies, not definitive first-download dates.

`parent_id` is the release-local normalized ChEMBL molecule identifier;
`target_id` is the ChEMBL target identifier. Their pair is the observation unit.
`cutoff` is 2018, 2020 or 2022. `Y_future_36m` is a legacy column name for the
frozen first-observed qualifying appearance label; 2022 is censored and does
not have 36 completed months. Use AUTHORITATIVE_U2_LABELS and the label audit
to interpret older score-table labels. `M*` columns are fixed model scores;
MODEL_SPECIFICATIONS.tsv maps them to models. `subset` names the evaluated
universe; NEW/SAME/UNKNOWN refer to valid scaffold comparisons. UNKNOWN is
never pooled with NEW. `full_boundary_tie_fde` evaluates all candidates tied at
the nominal top-0.1% boundary; use actual_evaluated_n as its selected count.
`expected_recall_at_1` averages within-compound positive top-tie fractions.
`mrr_tie_expectation` is tie-exact mean reciprocal rank. `top_tie_any_hit_fraction`
is set coverage, not rank-1 discrimination. Chemical A–G class definitions,
denominators and witness selection are documented in the Fig4a audit.

Source workbook panel tabs expose year, model/group, metric, value, denominator,
selected n, source and uncertainty fields. NOT_PLOTTED marks unavailable or
unplotted bounds, not a numerical confidence interval. Figure-level sheets map
panels to source rows. SHA256 values refer to exact files. Historical archive
UNKNOWN fields describe missing provenance and must not be substituted with
materialized-edge checksums.
