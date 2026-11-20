"""Explain verified eligibility changes using actual counter-release fields and frozen predicates."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import numpy as np
import pandas as pd

V = Path(__file__).resolve().parents[1]
A = V / 'structure_equivalence'
causes = pd.read_csv(A / 'tables/EDGE_CAUSE_ATTRIBUTION.tsv', sep='\t')
summary = []
edge_summary = []
for (year, direction), group in causes.groupby(['cutoff', 'direction']):
    selected = group[group.category.eq('FILTER_ELIGIBILITY_CHANGED')]
    evidence = pd.read_parquet(A / group.activity_evidence.iloc[0])
    keys = pd.MultiIndex.from_frame(selected[['parent_id', 'target_id']])
    evidence = evidence[pd.MultiIndex.from_frame(evidence[['parent_id', 'target_id']]).isin(keys)].copy()
    if evidence.empty:
        continue
    def text(field):
        return evidence['other_' + field].astype('string')
    failures = {}
    failures['MOLECULE_TYPE_MISSING'] = evidence.other_molecule_type.isna()
    failures['MIXTURE'] = text('molecule_type').str.lower().eq('mixture').fillna(False)
    comment = text('data_validity_comment')
    failures['VALIDITY_COMMENT'] = ~(comment.isna() | comment.str.strip().isin(['', '<NA>', 'None', 'nan']))
    failures['POTENTIAL_DUPLICATE'] = pd.to_numeric(evidence.other_potential_duplicate, errors='coerce').fillna(0).ne(0)
    human = pd.to_numeric(evidence.other_tax_id, errors='coerce').eq(9606) | text('target_organism').eq('Homo sapiens')
    failures['HUMAN_FILTER'] = ~human.fillna(False)
    failures['DIRECT_RELATIONSHIP'] = ~text('relationship_type').eq('D').fillna(False)
    failures['EXACT_STANDARD_RELATION'] = ~text('standard_relation').eq('=').fillna(False)
    failures['LOOSE_CONFIDENCE_GE7'] = ~pd.to_numeric(evidence.other_confidence_score, errors='coerce').ge(7)
    failures['ASSAY_B_OR_F'] = ~text('assay_type').isin(['B', 'F'])
    failures['PCHEMBL_MISSING'] = evidence.other_pchembl_value.isna()
    units = text('standard_units').str.replace('µ', 'u', regex=False).str.lower()
    failures['STANDARD_UNITS'] = ~units.isin(['nm', 'um', 'mm', 'm'])
    failures['STANDARD_TYPE'] = ~text('standard_type').isin(['Ki', 'Kd', 'IC50', 'EC50', 'Potency', 'AC50'])
    failures['SINGLE_PROTEIN'] = ~text('target_type').eq('SINGLE PROTEIN').fillna(False)
    failures['PARENT_MISSING'] = evidence.other_parent_id.isna()
    failures['TARGET_MISSING'] = evidence.other_target_id.isna()
    flags = pd.DataFrame(failures, index=evidence.index).fillna(False).astype(bool)
    evidence['prepared_failure_count'] = flags.sum(axis=1)
    evidence['only_prepared_failure_is_missing_molecule_type'] = flags.MOLECULE_TYPE_MISSING & flags.sum(axis=1).eq(1)
    evidence['prepared_failure_reasons'] = ['|'.join(flags.columns[row]) for row in flags.to_numpy()]
    unavailable = ~evidence.other_export_present
    assert not unavailable.any()
    assert evidence.loc[~evidence.other_prepared_eligible.astype(bool), 'prepared_failure_count'].gt(0).all()
    for rule in flags:
        subset = evidence[flags[rule]]
        summary.append({'cutoff': year, 'direction': direction, 'rule': rule, 'support_records': len(subset), 'distinct_edges': len(subset[['parent_id', 'target_id']].drop_duplicates()), 'classified_filter_edges': len(selected), 'interpretation': 'Nonexclusive actual predicate failures on comparison-release records'})
    edge = evidence.groupby(['parent_id', 'target_id']).only_prepared_failure_is_missing_molecule_type.all()
    edge_summary.append({'cutoff': year, 'direction': direction, 'filter_edges': len(selected), 'all_support_records_blocked_only_by_missing_molecule_type': int(edge.sum()), 'fraction_of_filter_edges': float(edge.sum()/len(selected))})
    evidence[['activity_id', 'parent_id', 'target_id', 'other_molecule_type', 'prepared_failure_count', 'prepared_failure_reasons', 'only_prepared_failure_is_missing_molecule_type']].to_parquet(A / 'derived' / f'qualification_failure_evidence_{year}_{direction}.parquet', index=False)
pd.DataFrame(summary).to_csv(V / 'tables/QUALIFICATION_FAILURE_COUNTS.tsv', sep='\t', index=False)
pd.DataFrame(edge_summary).to_csv(V / 'tables/MOLECULE_TYPE_MISSINGNESS_EVIDENCE.tsv', sep='\t', index=False)
(V / 'logs/qualification_evidence_complete.json').write_text(json.dumps({'status': 'COMPLETE', 'scope': 'Actual comparison-release fields on all source-support records of verified FILTER_ELIGIBILITY_CHANGED edges; frozen rule unchanged'}, indent=2))
print(pd.DataFrame(edge_summary).to_string(index=False))
