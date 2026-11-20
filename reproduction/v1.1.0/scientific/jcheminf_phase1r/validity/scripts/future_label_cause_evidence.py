"""Locate verified difference mechanisms within independently rebuilt future-positive labels."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import pandas as pd
import pyarrow.parquet as pq

V = Path(__file__).resolve().parents[1]
Q = V.parent
A = V / 'structure_equivalence'
causes = pd.read_csv(A / 'tables/EDGE_CAUSE_ATTRIBUTION.tsv', sep='\t')
summary = []
missingness = []
for year in [2018, 2020, 2022]:
    c = causes[causes.cutoff.eq(year) & causes.direction.eq('rollback_only')].copy()
    failure = pd.read_parquet(A / 'derived' / f'qualification_failure_evidence_{year}_rollback_only.parquet')
    missing_edges = failure.groupby(['parent_id', 'target_id']).only_prepared_failure_is_missing_molecule_type.all()
    for domain, stem in [('U2', 'REBUILT_U2_LABELS'), ('HISTORICAL_ACTIVE_OBSERVED', 'REBUILT_HISTORICAL_ACTIVE')]:
        positives = []
        for batch in pq.ParquetFile(Q / 'derived' / f'{stem}_{year}.parquet').iter_batches(columns=['parent_id', 'target_id', 'first_discovery']):
            f = batch.to_pandas()
            positives.append(f[f.first_discovery.eq(1)])
        p = pd.concat(positives, ignore_index=True).merge(c, on=['parent_id', 'target_id'], how='left', validate='one_to_one')
        p['category'] = p.category.fillna('NOT_ERASED')
        p['all_support_blocked_only_by_missing_molecule_type'] = pd.MultiIndex.from_frame(p[['parent_id', 'target_id']]).isin(missing_edges[missing_edges].index)
        assert p.loc[p.all_support_blocked_only_by_missing_molecule_type, 'category'].eq('FILTER_ELIGIBILITY_CHANGED').all()
        erased = int(p.category.ne('NOT_ERASED').sum())
        for category, count in p.category.value_counts().items():
            summary.append({'cutoff': year, 'domain': domain, 'category': category, 'count': int(count), 'future_positives': len(p), 'erased': erased, 'fraction_of_all_positives': count/len(p), 'fraction_of_erased': count/erased if category != 'NOT_ERASED' else None})
        missingness.append({'cutoff': year, 'domain': domain, 'future_positives': len(p), 'erased': erased, 'verified_filter_change_erased': int(p.category.eq('FILTER_ELIGIBILITY_CHANGED').sum()), 'missing_molecule_type_only_erased': int(p.all_support_blocked_only_by_missing_molecule_type.sum())})
        p.to_parquet(V / 'derived' / f'FUTURE_POSITIVE_CAUSE_EVIDENCE_{domain}_{year}.parquet', index=False)
pd.DataFrame(summary).to_csv(V / 'tables/FUTURE_LABEL_CAUSE_SUMMARY.tsv', sep='\t', index=False)
pd.DataFrame(missingness).to_csv(V / 'tables/FUTURE_LABEL_MISSINGNESS_EVIDENCE.tsv', sep='\t', index=False)
(V / 'logs/future_label_cause_complete.json').write_text(json.dumps({'status': 'COMPLETE', 'scope': 'Verified source-activity explanations among the independently rebuilt future-positive labels; raw-record presence differs from qualifying-graph membership'}, indent=2))
print(pd.DataFrame(missingness).to_string(index=False))
