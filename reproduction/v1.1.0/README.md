# Historical ChEMBL graph and label audit

This v1.1.0 release package reruns the existing fixed research implementation. It reconstructs eight historical qualifying graphs, three document-year rollback graphs, candidate universes, first-observed qualifying labels and record-level mechanisms. It does not train models or independently reimplement the algorithms. Release publication was authorized by the project owner on 20 November 2026. The DOI is provided by the published archive record; none is invented in this package.

## Database releases and provenance

ChEMBL 24.1, 27, 29, 31, 33, 35, 36 and 37 are required. `manifests/RAW_INPUTS.tsv` records each exact archive and extracted SQLite relative filename, byte size and SHA256. `manifests/RELEASE_SOURCES.tsv` preserves the release-specific download and checksum URLs. Its early download-state columns are historical snapshots; the current status is in its final recovery fields. Recovery hashes identify recovered files, not files proven to have been used at original study runtime. ChEMBL 24.1 contains chembl_24.db and a ChEMBL_24 VERSION row. Archive identity comes from the patch archive checksum and release notes.

## External data and licensing

Obtain the exact original releases legally from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/ using the explicit URLs in RELEASE_SOURCES.tsv. Verify publisher SHA256 before extracting. These instructions do not download anything during execution. Place each archive and extracted database under the relative paths in RAW_INPUTS.tsv; pass that parent directory explicitly with --raw-root. Current ChEMBL APIs cannot substitute for historical dumps. Raw SQLite and raw tar.gz files are excluded from this candidate by project policy. ChEMBL data and adaptations retain CC BY-SA 3.0 attribution/share-alike conditions, separately from MIT project code. See DATA_LICENSE.md, LICENSE_AUDIT.md and third-party notices.

## Installation and environment

The validated platform is macOS arm64, Python 3.12.14. Linux and other Python versions are not certified. Install an isolated environment with the pinned direct and transitive versions:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install --no-deps -r environment/requirements-installed-ssd.lock
.venv/bin/python -m pip check
```

Installation needs package-index access unless wheels are provided locally. Scientific execution performs no runtime dependency downloads. SQLite is supplied by Python; use local filesystem access supporting read-only SQLite and symbolic links. Keep at least 25 GB free for output plus the external raw inputs (230.4 GB in decimal units including archives and extracted databases; exact bytes are in RAW_INPUTS.tsv). Runtime memory should be at least 16 GB; the large candidate comparisons and canonical structure evidence can use more. Resource figures are planning estimates, not portable benchmarks.

## Reproduction commands

Run from this directory. Replace the three example paths with your actual environment, legally acquired raw root and an empty output directory. No author filesystem path is a default. Clear the two reserved optional research-script overrides shown below; the certified run used both unset. The orchestrator supplies its own canonical-mapping overrides only for the corresponding stages.

```sh
unset ENTITY_AUDIT_ROOT ENTITY_MAPPING_MODE
.venv/bin/python -I -B reproduce.py configure --raw-root /path/to/chembl_raw --run-root /path/to/empty_run
.venv/bin/python -I -B reproduce.py check --run-root /path/to/empty_run
.venv/bin/python -I -B reproduce.py smoke --run-root /path/to/empty_run
.venv/bin/python -I -B reproduce.py run --run-root /path/to/empty_run
```

Configuration copies fixed scientific sources into the run directory, materializes comparison-only expected data and links explicitly supplied SQLite inputs. `run` requires configure, then checks all raw SHA256 values before export. It stops on a discrepancy. Do not edit filters, labels or thresholds to make a check pass. A run directory cannot be overwritten by configure. Start a fresh empty directory for a repeat.

## Validation and expected outputs

`validation/CHECKS.json` records environment, raw-input, complete graph identities, qualifying activity rows, membership, candidate/label flags and record-level mechanism comparisons. `logs/clean_progress.json` reports the current stage. Per-stage logs and the rebuilt outputs are under the run root. A PASS requires all selected values after lossless Arrow storage normalization, not equal counts alone; original and observed file hashes are both retained. Parquet physical hashes may differ while typed values agree. The baseline exact mapping permits only the documented added mapping_rule annotation and still compares every original field.

Expected strict counts are 434456/499057/502711, rollback counts 459965/525328/563382, U2 candidates 9294785/10318524/12123841, U2 positives 946/1413/2104 and erased positives 633/932/1594 for 2018/2020/2022. The 2022 rollback-only denominator is 132801; eligibility-change count 97115 and molecule_type-only count 89778. The latter mechanism accounts for only 10 of 1594 erased U2 positives. Whole-graph and positive-label denominators are distinct.

## Figures and tables

`manifests/FIGURE_SOURCES.csv` links the frozen source tables for Figures 1–5 and S1–S12. The source hashes are checked; this entrypoint does not regenerate formal artwork or retrain models. Supplemental audit tables and original expected record-level data are bundled under expected/. The fixed raw chain reconstructs core graph/label/mechanism tables. Other stored model and chemistry results remain frozen; a dedicated portable plotting workflow is not certified here.

## Known limitations

Only sparse recovered release membership is observed; first appearance is not biological discovery. FTP mtimes are operational proxies. Full SQLite integrity scans are inherited checks anchored by fresh equal hashes, not rerun here; the retrospective 24.1 checkpoint has weaker provenance than directly persisted returns. Model comparisons are single-seed. Confidence and chemical decisions include post-outcome/post hoc choices. Stable entities are an all-eight-release post hoc subset. UNRESOLVED cause classifications are retained. Publication scope is the reviewed code, derived data, provenance and validation evidence. Manuscript prose and supplementary DOCX are excluded from this release; their editorial approval remains separate.

## Public payload boundary

FILE_MANIFEST.csv defines the public files. Exclude macOS ._ AppleDouble files, .DS_Store, local virtual environments and run directories when preparing an archive. Filesystem metadata were separately reviewed and are not scientific dependencies. Use the manifest to verify every intended payload file. After package finalization, run `shasum -a 256 -c SHA256SUMS` on macOS, or `sha256sum -c SHA256SUMS` on Linux, from this directory. Linux execution of the scientific chain remains uncertified.
