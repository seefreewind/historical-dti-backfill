#!/usr/bin/env python3
"""Export eligible ChEMBL activity rows from one historical SQLite release.

The export is a superset. Tier assignment happens later so confidence
thresholds can be counted before they are frozen.
"""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


STANDARD_TYPES = (
    "Ki",
    "Kd",
    "IC50",
    "EC50",
    "Potency",
    "AC50",
)


def columns(conn, table):
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}


def require(present, needed, table):
    missing = [name for name in needed if name not in present]
    if missing:
        raise SystemExit(f"{table} missing columns: {missing}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    activity_cols = columns(conn, "activities")
    assay_cols = columns(conn, "assays")
    target_cols = columns(conn, "target_dictionary")
    molecule_cols = columns(conn, "molecule_dictionary")
    doc_cols = columns(conn, "docs")
    require(activity_cols, ["activity_id", "assay_id", "molregno", "standard_type", "standard_relation", "standard_value", "standard_units"], "activities")
    require(assay_cols, ["assay_id", "doc_id", "tid", "assay_type", "confidence_score", "relationship_type"], "assays")
    require(target_cols, ["tid", "chembl_id", "target_type", "pref_name", "organism"], "target_dictionary")
    require(molecule_cols, ["molregno", "chembl_id"], "molecule_dictionary")
    require(doc_cols, ["doc_id", "year"], "docs")

    pchembl = "a.pchembl_value" if "pchembl_value" in activity_cols else "NULL AS pchembl_value"
    duplicate = "a.potential_duplicate" if "potential_duplicate" in activity_cols else "NULL AS potential_duplicate"
    validity = "a.data_validity_comment" if "data_validity_comment" in activity_cols else "NULL AS data_validity_comment"
    tax = "td.tax_id" if "tax_id" in target_cols else "NULL"
    max_phase = "md.max_phase" if "max_phase" in molecule_cols else "NULL AS max_phase"
    mol_type = "md.molecule_type" if "molecule_type" in molecule_cols else "NULL AS molecule_type"
    pref = "md.pref_name" if "pref_name" in molecule_cols else "NULL AS pref_name"
    prodrug = "md.prodrug" if "prodrug" in molecule_cols else "NULL AS prodrug"
    doc_chembl = "d.chembl_id" if "chembl_id" in doc_cols else "NULL"
    pubmed = "d.pubmed_id" if "pubmed_id" in doc_cols else "NULL"
    doc_type = "d.doc_type" if "doc_type" in doc_cols else "NULL"
    doi = "d.doi" if "doi" in doc_cols else "NULL"

    hierarchy = "molecule_hierarchy" in {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if hierarchy:
        parent_join = """
        LEFT JOIN molecule_hierarchy h ON h.molregno = a.molregno
        LEFT JOIN molecule_dictionary parent_md ON parent_md.molregno = h.parent_molregno
        """
        parent_select = "parent_md.chembl_id AS parent_molecule_chembl_id"
    else:
        parent_join = ""
        parent_select = "NULL AS parent_molecule_chembl_id"

    component_ok = {"target_components", "component_sequences"}.issubset(
        {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    )
    if component_ok:
        uniprot_join = """
        LEFT JOIN (
            SELECT tc.tid AS tid, GROUP_CONCAT(cs.accession, '|') AS uniprot_accessions
            FROM target_components tc
            JOIN component_sequences cs ON cs.component_id = tc.component_id
            GROUP BY tc.tid
        ) uni ON uni.tid = td.tid
        """
        uniprot_select = "uni.uniprot_accessions"
    else:
        uniprot_join = ""
        uniprot_select = "NULL AS uniprot_accessions"

    class_tables = {"target_components", "component_class", "protein_classification"}.issubset(
        {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    )
    class_cols = columns(conn, "protein_classification") if class_tables else set()
    if class_tables and "class_level" in class_cols:
        class_join = """
        LEFT JOIN (
            SELECT tid, GROUP_CONCAT(pref_name, '|') AS protein_class_l2
            FROM (
                SELECT tc.tid AS tid,
                       pc.pref_name AS pref_name,
                       pc.class_level AS class_level,
                       MAX(pc.class_level) OVER (PARTITION BY tc.tid) AS max_level
                FROM target_components tc
                JOIN component_class cc ON cc.component_id = tc.component_id
                JOIN protein_classification pc ON pc.protein_class_id = cc.protein_class_id
            )
            WHERE class_level = max_level
            GROUP BY tid
        ) cls ON cls.tid = td.tid
        """
        class_select = "cls.protein_class_l2"
    else:
        class_join = ""
        class_select = "NULL AS protein_class_l2"

    types = ",".join(f"'{name}'" for name in STANDARD_TYPES)
    query = f"""
    SELECT
        a.activity_id,
        a.assay_id,
        a.molregno,
        md.chembl_id AS molecule_chembl_id,
        {parent_select},
        {pref} AS molecule_pref_name,
        {max_phase} AS max_phase,
        {mol_type} AS molecule_type,
        {prodrug} AS prodrug,
        td.chembl_id AS target_chembl_id,
        td.target_type,
        td.pref_name AS target_pref_name,
        td.organism AS target_organism,
        {tax} AS tax_id,
        {uniprot_select},
        {class_select},
        ass.assay_type,
        ass.confidence_score,
        ass.relationship_type,
        a.standard_type,
        a.standard_relation,
        CAST(a.standard_value AS REAL) AS standard_value,
        a.standard_units,
        {pchembl},
        {duplicate},
        {validity},
        d.doc_id AS document_id,
        {doc_chembl} AS document_chembl_id,
        d.year AS document_year,
        {pubmed} AS pubmed_id,
        {doi} AS doi,
        {doc_type} AS doc_type
    FROM activities a
    JOIN assays ass ON ass.assay_id = a.assay_id
    JOIN target_dictionary td ON td.tid = ass.tid
    JOIN molecule_dictionary md ON md.molregno = a.molregno
    LEFT JOIN docs d ON d.doc_id = ass.doc_id
    {parent_join}
    {uniprot_join}
    {class_join}
    WHERE td.target_type = 'SINGLE PROTEIN'
      AND ass.assay_type IN ('B', 'F')
      AND ass.confidence_score >= 7
      AND a.standard_type IN ({types})
      AND a.standard_value IS NOT NULL
      AND a.standard_relation IN ('=', '<', '<=', '>', '>=', '~')
    """
    if __import__("os").environ.get("PRINT_SQL"):
        print(query)
        return
    cursor = conn.execute(query)
    names = [desc[0] for desc in cursor.description]
    float_cols = {"standard_value", "pchembl_value", "max_phase"}
    int_cols = {
        "activity_id", "assay_id", "molregno", "prodrug", "tax_id",
        "confidence_score", "potential_duplicate", "document_id",
        "document_year", "pubmed_id",
    }
    fields = []
    for name in names:
        if name in float_cols:
            fields.append(pa.field(name, pa.float64()))
        elif name in int_cols:
            fields.append(pa.field(name, pa.int64()))
        else:
            fields.append(pa.field(name, pa.string()))
    schema = pa.schema(fields)
    writer = pq.ParquetWriter(out, schema, compression="zstd")
    total = 0
    while True:
        rows = cursor.fetchmany(200_000)
        if not rows:
            break
        column_values = {name: [] for name in names}
        for row in rows:
            for name, value in zip(names, row):
                if name in int_cols and isinstance(value, float):
                    value = None if value != value else int(value)
                column_values[name].append(value)
        arrays = []
        for field in schema:
            arrays.append(pa.array(column_values[field.name], type=field.type, from_pandas=True))
        writer.write_table(pa.Table.from_arrays(arrays, schema=schema))
        total += len(rows)
        print(f"{args.version}\t{total}", flush=True)
    writer.close()
    print(f"DONE {args.version} {total} {out}", flush=True)


if __name__ == "__main__":
    main()
