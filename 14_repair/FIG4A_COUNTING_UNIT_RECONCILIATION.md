# Figure 4a counting-unit reconciliation

FINAL VERDICT:
PASS_FIG4A_UNIT_RECONCILIATION

COUNTING UNIT:
COUNTING_UNIT_CONFIRMED_CANDIDATE_TARGET

## Scope and disposition

This report resolves the Figure 4a reporting-unit discrepancy from the existing frozen classifier output. It supplies the exact replacement text and source-data metadata for the next manuscript integration. No V5 was generated, no frozen input was edited, and no classifier, model, label-building or Phase 1–4 pipeline was run. V4 DOCX and archived versions remain unchanged. The replacement text below supersedes their disputed unit descriptions for the purpose of this reconciliation; those files are not represented as already corrected publication artifacts.

## Authoritative unit check

```text
2018:
rows = 2600
unique candidate-target = 2600
unique candidate-target-witness = 2600
class sum = 2600
reported denominator = 2600
```

```text
2020:
rows = 3391
unique candidate-target = 3391
unique candidate-target-witness = 3391
class sum = 3391
reported denominator = 3391
```

```text
2022:
rows = 3478
unique candidate-target = 3478
unique candidate-target-witness = 3478
class sum = 3478
reported denominator = 3478
```

All three required equalities pass. Uniqueness of the three-field key alone would not establish the unit; uniqueness of the candidate–target key and the classifier aggregation jointly establish it.

## Multiple-ligand witness audit

| Cutoff | One matching ligand | More than one | Zero | Median | Maximum |
|---|---:|---:|---:|---:|---:|
| 2018 | 2,086 | 514 | 0 | 1 | 18 |
| 2020 | 2,768 | 623 | 0 | 1 | 18 |
| 2022 | 2,855 | 623 | 0 | 1 | 17 |

A candidate–target pair with several matching historical ligands contributes one classified row. Witness multiplicity does not multiply observational units. All recorded rows have a nonempty classification reason and selected historical-ligand identifier.

## Class-selection semantics

Read-only code inspection of `scripts/phase4_metrics_chemistry.py`, function `run_chemistry`, confirms:

1. Candidate rows are selected from the fingerprint-identical band.
2. Historical ligands are indexed by candidate target and fingerprint.
3. Every matching historical ligand is evaluated by the molecular-relationship rules.
4. One call is selected by A–G class priority. Ties within a class are broken by historical-ligand identifier using `min(calls, key=lambda x: (x[0], x[3]))`.
5. Exactly one output row is appended for the candidate–target pair. Summary counts aggregate these rows.

The selected ligand is a classification witness. It is not randomly chosen. The summary is not a census of all candidate–ligand relationships. Class-priority selection means each percentage describes candidate–target pairs assigned to a class under this deterministic rule.

## Frozen class counts and arithmetic

| Cutoff | A | B | C | D | E | F | G | Total | B (%) | C (%) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2018 | 25 | 1107 | 722 | 0 | 0 | 642 | 104 | 2600 | 42.6 | 27.8 |
| 2020 | 3 | 1516 | 947 | 0 | 0 | 803 | 122 | 3391 | 44.7 | 27.9 |
| 2022 | 4 | 1581 | 1002 | 0 | 0 | 768 | 123 | 3478 | 45.5 | 28.8 |

All seven row-level class counts match the frozen aggregate table, including D/E zeros and G counts 104, 122 and 123. Stored fractions agree with count/denominator to an absolute tolerance of 10⁻¹². One-decimal percentages remain 42.6%, 44.7%, 45.5% for defined stereo differences and 27.8%, 27.9%, 28.8% for stereo-specification differences. No reclassification was performed.

## Exact reporting-unit corrections for manuscript integration

These replacements are limited to the reconciled unit. Integration into a new Word manuscript and rendered figure is deferred to the separately authorized resumption step; they do not authorize broader V5 work in this turn.

### Results — replacement opening paragraph for Figure 4a

The fingerprint-identical band contained 2,600, 3,391 and 3,478 candidate compound–target pairs in 2018, 2020 and censored 2022, respectively (Fig. 4a). Each candidate–target pair was assigned one molecular-relationship class after evaluating its equal-fingerprint historical ligands. Defined stereochemical differences accounted for 1,107, 1,516 and 1,581 pairs (42.6%, 44.7% and 45.5%); stereo-specification differences accounted for 722, 947 and 1,002 pairs (27.8%, 27.9% and 28.8%). Exact canonical-structure counts were 25, 3 and 4, and fingerprint-collision counts were 642, 803 and 768. Complete A–G counts and classification witnesses are reported in the chemical classification tables.

### Figure 4a panel title and y-axis

Panel title: Fingerprint-identical candidate–target relationships

Full title, where space permits: Molecular relationship classes of fingerprint-identical candidate–target pairs

Y-axis: Fingerprint-identical candidate–target pairs (%)

### Figure 4 caption — panel a description

a, Seven mutually exclusive molecular-relationship classes assigned to fingerprint-identical candidate compound–target pairs. A, exact canonical structure; B, defined stereochemical difference; C, unequal or incomplete stereochemical specification; D, isotope or other InChI-layer difference; E, standardization or representation difference; F, fingerprint collision; G, unresolved.

### Figure 4 caption — replacement denominator explanation

For panel a, denominators are 2,600, 3,391 and 3,478 fingerprint-identical candidate compound–target pairs in 2018, 2020 and censored 2022, respectively. For each candidate–target pair, all equal-fingerprint historical ligands of the candidate target were evaluated. One molecular-relationship class and its corresponding historical-ligand witness were retained using A–G class priority, with ties within a class broken by historical-ligand identifier. The denominators count candidate–target pairs; the recorded historical ligand is the classification witness.

Delete the sentence “These are molecule-pair counts, separate from U2 candidate counts.” Retain the rest of panels b–d, their counts, censoring, post hoc disclosure and uncertainty wording.

### Methods — Chemical identity audit opening

Each fingerprint-identical U2 candidate compound–target pair was classified once. For a candidate target with multiple equal-fingerprint historical ligands, all such ligands were evaluated under the molecular-relationship rules, after which one class and corresponding historical-ligand witness were retained using deterministic A–G class priority, with ties within a class broken by historical-ligand identifier.

Continue with the existing molecular features and definitions of classes A–G. Replace “All U2 pairs in the fingerprint-identical band are reclassified...” with this opening. Keep the existing definitions and every numeric result.

### Supplementary Methods / table description — clarification

The chemical classification table contains one row per fingerprint-identical candidate compound–target pair. All equal-fingerprint historical ligands of the candidate target are evaluated. A–G class priority and historical-ligand identifier select one class and one classification witness. Class fractions use the number of classified candidate–target pairs at that cutoff as the denominator.

### Source Data metadata and field mapping

The current main Source Data provides Figure 4a aggregates rather than row-level witnesses. The 21 Fig4a rows and the corresponding 21 rows in the consolidated Figure4 sheet were checked: counts, denominators and class percentages all agree with the frozen summary. The existing `denominator_definition` describes candidate–historical-ligand pairs and must be replaced in both locations during controlled integration. No numeric workbook cell needs to change.

Required aggregate metadata:

```text
denominator_unit = candidate_compound_target_pair
denominator_definition = Fingerprint-identical candidate compound–target pairs; one assigned class per pair
classification_rule = A–G priority; ties within a class resolved by historical-ligand identifier
historical_ligand_role = classification_witness
row_level_source = 14_repair/TANIMOTO1_RECLASSIFIED.tsv
```

Required row-level field mapping (all seven source fields already exist; no relationships need to be recomputed):

| Publication field | Existing frozen column |
|---|---|
| cutoff | cutoff |
| candidate_parent_id | parent_id |
| candidate_target_id | target_id |
| assigned_class | chemical_class |
| matched_historical_ligand | matched_historical_ligand |
| n_equal_fingerprint_ligands | n_equal_fingerprint_ligands |
| classification_reason | classification_reason |

The aggregate `SOURCE_DATA_MANIFEST.tsv` Figure4/a entry links the correct frozen classification summary. Add the metadata above to that entry when integrating the corrected workbook. The audit TSV produced in this turn already records `denominator_unit=candidate_compound_target_pair` for every cutoff. No claim is made that the old publication workbook has already been regenerated.

## Occurrence sweep and synchronization map

The checked active V4 main-text copies describe the band as candidate pairs but need the explicit compound–target unit and deterministic-selection sentence above. Both active V4 legend copies contain the incorrect molecule-pair sentence. Both integrated supplementary-information copies and the V3 supplementary-methods source describe checking all equal-fingerprint ligands; append the one-row/one-selected-witness clarification above. V3 supplementary results and legends were inspected as historical supporting materials; archive snapshots should remain identifiable and should not supersede this reconciliation.

The Figure 4 plotting source currently uses “Molecular relationship classes A–G” and “Tanimoto = 1 pairs (%)”. Substitute the panel title and y-axis above when exporting the corrected figure; its plotted fractions are unchanged. The caption-building source also contains the incorrect sentence and must receive the same caption replacement so regeneration cannot reintroduce it.

Locations found in the text sweep:

- `13_manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:28` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:58` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:62` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:82` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:154` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:28` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:58` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:62` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:82` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/MAIN_MANUSCRIPT_V4_NCOMMS.md:154` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/FIGURE_LEGENDS_V4_NCOMMS.md:25` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/FIGURE_LEGENDS_V4_NCOMMS.md:27` — SUPERSEDED UNIT DESCRIPTION.
- `15_ncomms_submission/manuscript/FIGURE_LEGENDS_V4_NCOMMS.md:25` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/FIGURE_LEGENDS_V4_NCOMMS.md:27` — SUPERSEDED UNIT DESCRIPTION.
- `13_manuscript/SUPPLEMENTARY_INFORMATION_NCOMMS_V1.md:196` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `15_ncomms_submission/manuscript/SUPPLEMENTARY_INFORMATION_NCOMMS_V1.md:196` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/FIGURE_LEGENDS_V3.md:19` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.
- `13_manuscript/FIGURE_LEGENDS_V3.md:43` — CLARIFY / CONFIRM CANDIDATE–TARGET UNIT.

Additional checked locations: `Source_Data.xlsx`, sheets Fig4a and Figure4, column `denominator_definition`; `source_data/SOURCE_DATA_MANIFEST.tsv`, Figure4/a; `scripts/phase4_figures.py`, Figure4a labels; `scripts/prepare_ncomms_v1.py`, Figure4 caption text. The authoritative DOCX's disputed Figure 4 caption is covered by the same replacement. Archived reviewer-package and isolated-reproduction copies are preserved snapshots, not updated submission artifacts.

## Scientific interpretation

Among fingerprint-identical candidate–target pairs, defined stereochemical differences account for approximately 43–46%, and stereo-specification differences for approximately 28%. Fingerprint identity does not imply molecular identity. The near-neighbour chemotype interpretation is unchanged. These are not proportions of all molecular relationships, all candidate–ligand combinations or all historical ligands.

## Input integrity

SHA256 before and after this audit:

- `14_repair/TANIMOTO1_RECLASSIFIED.tsv`: `f932748dcff93bb98e69e48c503a92ae15b75c321724d34158a51d05fcc6d676` — UNCHANGED.
- `14_repair/CHEMICAL_CLASSIFICATION_SUMMARY.tsv`: `ba9e399c5581632bb22e31c647f99cab719f8ec0ba7e1727a4a7b5569bdc0a02` — UNCHANGED.
- `scripts/phase4_metrics_chemistry.py`: `d13b5cdeff6486b9e70a452cec203353a98ee57624c6584b51395cb618e09d65` — UNCHANGED.

## Pass conditions

- One output row equals one unique candidate–target pair: PASS.
- Class sums equal the reported denominators: PASS.
- All existing A–G counts and fractions are unchanged: PASS.
- Figure 4a Source Data numerical values agree with the frozen aggregates: PASS.
- Caption and Methods wording can be corrected without changing the analysis: PASS; exact replacements recorded above.
- Publication-file integration status: DEFERRED; no V5, Word rebuild, figure export or workbook regeneration was performed in this audit-only turn.

SCIENTIFIC RESULT CHANGED:
NO

FROZEN NUMBERS CHANGED:
NO

ERROR TYPE:
REPORTING-UNIT DESCRIPTION

CORRECT OBSERVATIONAL UNIT:
CANDIDATE COMPOUND–TARGET PAIR

HISTORICAL LIGAND ROLE:
CLASSIFICATION WITNESS AFTER EVALUATION OF EQUAL-FINGERPRINT HISTORICAL LIGANDS

NEXT STEP:
RESUME_NCOMMS_FINAL_EDITORIAL_OPTIMIZATION_V1

STOP:
Audit complete. Await user confirmation before resuming V5 editorial optimization and applying the reconciled wording and metadata to publication artifacts.
