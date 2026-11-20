"""Record-level conservative attribution. Incomplete or mixed evidence remains UNRESOLVED."""
import os,sys,json,datetime,sqlite3
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1];Q=BASE.parent;R=Q.parent/'biomedical_ai_historical_discovery'
V=Path(os.environ.get('ENTITY_AUDIT_ROOT',str(BASE)))
for folder in ['derived','tables','logs']: (V/folder).mkdir(parents=True,exist_ok=True)
cfg=json.loads((R/'configs/repair_config.json').read_text())
def stage(name,v=''):
 d=dict(status='RUNNING',pid=os.getpid(),stage=name,version=v,updated_at=datetime.datetime.now().astimezone().isoformat())
 (V/'logs/cause_progress.json').write_text(json.dumps(d,indent=2));print(d,flush=True)
def pairs(f):return pd.MultiIndex.from_frame(f[['parent_id','target_id']])
def primary(v):return pd.read_parquet(Q/'derived'/('edges_24.1_2018.parquet' if v=='24.1' else f'primary_{v}.parquet'),columns=['parent_id','target_id'])
def external_doc_key(f):
 doi=f.doi.astype('string').str.lower().str.strip().str.replace(r'^https?://(dx\.)?doi.org/','',regex=True)
 pm=pd.to_numeric(f.pubmed_id,errors='coerce').astype('Int64').astype('string')
 return doi,pm
try:
 versions=['24.1','27','31','37'];prepared={};full={};docs={};cmaps={};tmaps={}
 for v in versions:
  stage('LOAD_RECORD_EVIDENCE',v)
  prepared[v]=pd.read_parquet(BASE/'derived'/f'qualifying_rows_{v}.parquet',columns=['activity_id','parent_id','target_id','document_year'])
  docs[v]=pd.read_parquet(BASE/'derived'/f'documents_{v}.parquet')
  m=pd.read_parquet(V/'derived'/f'compound_mapping_{v}.parquet',columns=['chembl_id','reference_id','mapping_status'])
  cmaps[v]=m.loc[m.mapping_status.eq('HIGH_CONFIDENCE_EXACT_STRUCTURE')].set_index('chembl_id').reference_id
  t=pd.read_parquet(V/'derived'/f'target_mapping_{v}.parquet')
  tmaps[v]=t.loc[t.mapping_status.eq('HIGH_CONFIDENCE_ACCESSION_AND_SEQUENCE')].set_index('local_id').reference_id
 # Identify all difference-supporting records before raw-presence scans.
 jobs=[];tracked=set()
 for year,spec in cfg['cutoffs'].items():
  strict=primary(spec['release']);retro=pd.read_parquet(Q/'derived'/('edges_37_2018.parquet' if year=='2018' else f'rollback_37_{year}.parquet'),columns=['parent_id','target_id'])
  si,ri=pairs(strict),pairs(retro)
  for direction,base,other,a,b in [('rollback_only','37',spec['release'],retro,si),('strict_only',spec['release'],'37',strict,ri)]:
   difference=a[~pairs(a).isin(b)]
   p=prepared[base]
   if direction=='rollback_only':p=p[p.document_year.le(int(year))]
   support=p[pairs(p).isin(pairs(difference))].copy()
   tracked.update(support.activity_id)
   jobs.append((year,direction,base,other,difference,support))
 for v in versions:
  stage('LOAD_TRACKED_EXPORTED_ACTIVITY_EVIDENCE',v)
  cols=['activity_id','molecule_chembl_id','parent_molecule_chembl_id','target_chembl_id','document_id','document_chembl_id','document_year','doi','pubmed_id','assay_type','confidence_score','standard_type','standard_relation','standard_units','pchembl_value','potential_duplicate','data_validity_comment','relationship_type','target_type','tax_id','target_organism','molecule_type']
  f=pd.read_parquet(Q/'derived'/f'activities_{v}.parquet',columns=cols)
  f['parent_id']=f.parent_molecule_chembl_id.fillna(f.molecule_chembl_id);f['target_id']=f.target_chembl_id
  f=f[f.activity_id.isin(tracked)].copy()
  f=f.drop(columns=['molecule_chembl_id','parent_molecule_chembl_id','target_chembl_id'])
  assert not f.activity_id.duplicated().any()
  f['prepared_eligible']=f.activity_id.isin(prepared[v].activity_id)
  full[v]=f.set_index('activity_id')
 del prepared
 raw_present={}
 cache_verified=False
 if V!=BASE and json.loads((BASE/'logs/cause_progress.json').read_text())['status']=='COMPLETE':
  cached_tracking=set()
  for year,direction,*_ in jobs:
   cached_tracking.update(pd.read_parquet(BASE/'derived'/f'edge_activity_evidence_{year}_{direction}.parquet',columns=['activity_id']).activity_id)
  assert cached_tracking==tracked, 'Raw-presence cache coverage differs from current supporting activity IDs'
  cache_verified=True
 for v in versions:
  stage('RAW_ACTIVITY_ID_PRESENCE_SCAN',v)
  cached=BASE/'derived'/f'raw_tracked_activity_presence_{v}.parquet'
  if cache_verified and cached.exists():
   present=set(pd.read_parquet(cached).activity_id)
   raw_present[v]=present
   pd.DataFrame({'activity_id':sorted(present)}).to_parquet(V/'derived'/f'raw_tracked_activity_presence_{v}.parquet',index=False)
   continue
  db=next(p for p in (Q/'raw'/('extract_'+v)).glob('*.db') if not p.name.startswith('._'))
  con=sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1',uri=True)
  con.execute('PRAGMA cache_size=-262144');present=set()
  for f in pd.read_sql_query('SELECT activity_id FROM activities',con,chunksize=500000):
   present.update(tracked.intersection(f.activity_id.to_numpy()))
  raw_present[v]=present;con.close()
  pd.DataFrame({'activity_id':sorted(present)}).to_parquet(V/'derived'/f'raw_tracked_activity_presence_{v}.parquet',index=False)
 summary=[];edge_outputs=[]
 for year,direction,base,other,difference,support in jobs:
  stage('CLASSIFY_DIRECT_EVIDENCE',year+'_'+direction)
  # Every source supporting row is retained in the evidence table. Other-export absence is distinct from raw absence.
  e=support[['activity_id','parent_id','target_id','document_year']].merge(full[base].reset_index(),on=['activity_id','parent_id','target_id','document_year'],how='left',validate='one_to_one')
  e=e.merge(full[other].reset_index().add_prefix('other_'),left_on='activity_id',right_on='other_activity_id',how='left',validate='one_to_one')
  e['other_raw_activity_present']=e.activity_id.isin(raw_present[other])
  e['other_export_present']=e.other_activity_id.notna()
  e['parent_reference']=e.parent_id.map(cmaps[base]);e['target_reference']=e.target_id.map(tmaps[base])
  e['other_parent_reference']=e.other_parent_id.map(cmaps[other]);e['other_target_reference']=e.other_target_id.map(tmaps[other])
  e['stable_semantics']=e.parent_reference.notna()&e.target_reference.notna()&e.parent_reference.eq(e.other_parent_reference)&e.target_reference.eq(e.other_target_reference)
  e['same_local_pair']=e.parent_id.eq(e.other_parent_id)&e.target_id.eq(e.other_target_id)
  doi,pm=external_doc_key(e);ddoi,dpm=external_doc_key(docs[other])
  e['document_external_identifier_present']=doi.notna()&doi.ne('') | pm.notna()
  e['document_found_in_other_release']=e.document_chembl_id.isin(docs[other].chembl_id) | doi.isin(set(ddoi.dropna())) | pm.isin(set(dpm.dropna()))
  e['source_pair_semantics_stable_across_releases']=e.parent_reference.notna()&e.target_reference.notna()&e.parent_reference.eq(e.parent_id.map(cmaps[other]))&e.target_reference.eq(e.target_id.map(tmaps[other]))
  e['late_old_document_evidence']=e.source_pair_semantics_stable_across_releases&(direction=='rollback_only') & ~e.other_raw_activity_present & ~e.document_found_in_other_release & e.document_external_identifier_present & e.document_year.le(int(year))
  e['compound_id_reassignment']=e.other_export_present&e.stable_semantics&~e.parent_id.eq(e.other_parent_id)&e.target_id.eq(e.other_target_id)
  e['target_id_reassignment']=e.other_export_present&e.stable_semantics&e.parent_id.eq(e.other_parent_id)&~e.target_id.eq(e.other_target_id)
  e['document_year_changed']=e.other_export_present&e.stable_semantics&e.same_local_pair&((e.document_year.fillna(-1)!=e.other_document_year.fillna(-1)))
  if direction=='strict_only':
   e['year_change_excludes_other']=e.document_year_changed & (e.other_document_year.isna() | e.other_document_year.gt(int(year))) & e.other_prepared_eligible.fillna(False)
  else:e['year_change_excludes_other']=False
  e['qualification_change']=e.other_export_present&e.stable_semantics&e.same_local_pair&(~e.other_prepared_eligible.fillna(False) | (e.confidence_score.ge(9)&e.other_confidence_score.lt(9)) | (e.assay_type.eq('B')&~e.other_assay_type.eq('B')))
  e['unchanged_year_excludes_other']=(direction=='strict_only')&e.other_export_present&e.stable_semantics&e.same_local_pair&e.other_prepared_eligible.fillna(False)&(e.other_document_year.isna()|e.other_document_year.gt(int(year)))&~e.document_year_changed
  # A missing comparison cannot establish a cause. pandas all() skips NA by default.
  evidence_flags=['late_old_document_evidence','compound_id_reassignment','target_id_reassignment','year_change_excludes_other','qualification_change','unchanged_year_excludes_other']
  for flag in evidence_flags:e[flag]=e[flag].fillna(False).astype(bool)
  evidence_path=V/'derived'/f'edge_activity_evidence_{year}_{direction}.parquet';e.to_parquet(evidence_path,index=False)
  grouped=e.groupby(['parent_id','target_id'],sort=False)
  aggregate=grouped.agg(support_records=('activity_id','size'),source_ids=('activity_id',lambda x:'|'.join(map(str,sorted(x)))),all_late=('late_old_document_evidence','all'),all_compound_id=('compound_id_reassignment','all'),all_target_id=('target_id_reassignment','all'),all_year_change=('year_change_excludes_other','all'),all_qualification=('qualification_change','all'),all_unchanged_year=('unchanged_year_excludes_other','all'),other_raw_present_count=('other_raw_activity_present','sum'),other_export_present_count=('other_export_present','sum')).reset_index()
  categories={'all_late':'LATE_OBSERVED_OLD_DOCUMENT_EVIDENCE','all_compound_id':'COMPOUND_IDENTIFIER_REASSIGNMENT','all_target_id':'TARGET_IDENTIFIER_REASSIGNMENT','all_year_change':'DOCUMENT_YEAR_CHANGED_OR_MISSING','all_qualification':'FILTER_ELIGIBILITY_CHANGED','all_unchanged_year':'OTHER_VERIFIED_DOCUMENT_YEAR_ROLLBACK_EXCLUSION'}
  aggregate['category']='UNRESOLVED';aggregate['evidence_rule']='INCOMPLETE_OR_MIXED_RECORD_EVIDENCE'
  # A single fully-supported explanation is required. Mixed explanations remain unresolved.
  mask=aggregate[list(categories)].sum(axis=1).eq(1)
  for col,category in categories.items():
   hit=mask&aggregate[col];aggregate.loc[hit,'category']=category;aggregate.loc[hit,'evidence_rule']='ALL_SOURCE_SUPPORTING_RECORDS_SATISFY_'+col
  out=difference.merge(aggregate,on=['parent_id','target_id'],how='left',validate='one_to_one')
  assert out.support_records.notna().all() and len(out)==len(difference)
  out['cutoff']=year;out['direction']=direction;out['source_release']=base;out['comparison_release']=other;out['activity_evidence']=str(evidence_path.relative_to(V))
  edge_outputs.append(out)
  for category,count in out.category.value_counts().items():summary.append(dict(cutoff=year,direction=direction,category=category,count=int(count),denominator=len(out),fraction=count/len(out)))
  pd.DataFrame(summary).to_csv(V/'tables/EDGE_CAUSE_SUMMARY.tsv',sep='\t',index=False)
 pd.concat(edge_outputs,ignore_index=True).to_csv(V/'tables/EDGE_CAUSE_ATTRIBUTION.tsv',sep='\t',index=False)
 stage('COMPLETE')
 p=V/'logs/cause_progress.json';d=json.loads(p.read_text());d['status']='COMPLETE';p.write_text(json.dumps(d,indent=2))
except Exception as exc:
 p=V/'logs/cause_progress.json';d=json.loads(p.read_text()) if p.exists() else {};d.update(status='FAILED',error=repr(exc));p.write_text(json.dumps(d,indent=2));raise
