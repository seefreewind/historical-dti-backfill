"""Conservative one-to-one structure and protein mapping, with global stable-subset sensitivity."""
import sys,os,json,datetime
sys.dont_write_bytecode=True
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
BASE=Path(__file__).resolve().parents[1];Q=BASE.parent;R=Q.parent/'biomedical_ai_historical_discovery'
V=Path(os.environ.get('ENTITY_AUDIT_ROOT',str(BASE)))
mode=os.environ.get('ENTITY_MAPPING_MODE','EXACT_PUBLISHER_STRING')
for folder in ['derived','tables','logs']: (V/folder).mkdir(parents=True,exist_ok=True)
cfg=json.loads((R/'configs/repair_config.json').read_text());versions=list(cfg['releases'])
def stage(name,v=''):
 d=dict(status='RUNNING',pid=os.getpid(),stage=name,version=v,updated_at=datetime.datetime.now().astimezone().isoformat())
 (V/'logs/entity_mapping_progress.json').write_text(json.dumps(d,indent=2));print(d,flush=True)
def path(v):return Q/'derived'/('edges_24.1_2018.parquet' if v=='24.1' else f'primary_{v}.parquet')
def read_molecules(v,ids):
 return pq.read_table(BASE/'derived'/f'molecules_{v}.parquet',filters=[('chembl_id','in',list(ids))]).to_pandas()
def target_summary(v):
 f=pd.read_parquet(BASE/'derived'/f'target_components_{v}.parquet')
 f=f[f.target_type.eq('SINGLE PROTEIN')].copy()
 out=[]
 for tid,g in f.groupby('chembl_id',sort=False):
  valid=(len(g)==1 and g.accession.notna().all() and g.sequence_md5sum.notna().all() and g.component_type.eq('PROTEIN').all() and g.component_tax_id.eq(9606).all() and g.homologue.eq(0).all())
  signature=str(g.accession.iloc[0])+'|'+str(g.sequence_md5sum.iloc[0]) if valid else ''
  out.append(dict(local_id=tid,tid=g.tid.iloc[0],accessions='|'.join(sorted(g.accession.dropna().unique())),sequence_md5='|'.join(sorted(g.sequence_md5sum.dropna().unique())),signature=signature,component_count=len(g),eligible=valid))
 return pd.DataFrame(out)
try:
 for v in versions:assert (BASE/'logs'/f'entities_{v}_complete.json').exists()
 frames={v:pd.read_parquet(path(v),columns=['parent_id','target_id']) for v in versions}
 allparents=set().union(*(set(f.parent_id) for f in frames.values()))
 # Compare all graph-used local entities, also if inactive in an individual release.
 stage('REFERENCE_STRUCTURES')
 ref=read_molecules('37',allparents)
 ref=ref[ref.standard_inchi_key.notna() & ref.canonical_smiles.notna()].copy()
 assert not ref.standard_inchi_key.duplicated().any()
 ref=ref[['chembl_id','standard_inchi_key','canonical_smiles']].rename(columns={'chembl_id':'reference_id','canonical_smiles':'reference_smiles'})
 reft=target_summary('37');validt=reft[reft.eligible].copy()
 duplicate_signatures=set(validt.loc[validt.signature.duplicated(keep=False),'signature'])
 validt[validt.signature.isin(duplicate_signatures)].to_csv(V/'tables/TARGET_REFERENCE_ONE_TO_MANY.tsv',sep='\t',index=False)
 validt=validt[~validt.signature.isin(duplicate_signatures)]
 canonical={}
 if mode=='CANONICAL_STRUCTURE_EQUIVALENCE':
  proof=pd.read_parquet(BASE/'derived/SMILES_REPRESENTATION_EQUIVALENCE.parquet').dropna()
  canonical=proof.set_index('publisher_smiles').rdkit_canonical_isomeric_smiles.to_dict()
 cmaps={};tmaps={};mapping_summary=[];common_c=None;common_t=None
 for v in versions:
  stage('BUILD_HIGH_CONFIDENCE_MAPPING',v)
  f=read_molecules(v,allparents)
  full=pd.DataFrame({'chembl_id':sorted(allparents)}).merge(f,on='chembl_id',how='left',validate='one_to_one')
  m=full.merge(ref,on='standard_inchi_key',how='left',validate='many_to_one')
  m['mapping_status']='UNMAPPED_STRUCTURE_OR_REFERENCE'
  m.loc[m.reference_id.notna(),'mapping_status']='SMILES_CHANGED_WITH_SAME_INCHIKEY'
  same=m.reference_id.notna() & m.canonical_smiles.eq(m.reference_smiles).fillna(False)
  m['mapping_rule']=np.where(same,'IDENTICAL_PUBLISHER_SMILES','UNCONFIRMED')
  if mode=='CANONICAL_STRUCTURE_EQUIVALENCE':
   source_canonical=m.canonical_smiles.map(canonical);reference_canonical=m.reference_smiles.map(canonical)
   equivalent=m.reference_id.notna() & source_canonical.notna() & source_canonical.eq(reference_canonical)
   m.loc[equivalent & ~same,'mapping_rule']='SAME_FULL_INCHIKEY_AND_RDKIT_ISOMERIC_GRAPH'
   same=same | equivalent
  m.loc[same,'mapping_status']='HIGH_CONFIDENCE_EXACT_STRUCTURE'
  duplicates=m.loc[same,'reference_id'].duplicated(keep=False)
  ambiguous=set(m.loc[same].loc[duplicates,'chembl_id'])
  m.loc[m.chembl_id.isin(ambiguous),'mapping_status']='AMBIGUOUS_MANY_TO_ONE'
  m.loc[m.molregno.isna(),'mapping_status']='LOCAL_ENTITY_ABSENT'
  m['same_chembl_id']=m.chembl_id.eq(m.reference_id)
  m.to_parquet(V/'derived'/f'compound_mapping_{v}.parquet',index=False)
  good=m.mapping_status.eq('HIGH_CONFIDENCE_EXACT_STRUCTURE')
  cmaps[v]=m.loc[good].set_index('chembl_id').reference_id.to_dict()
  for status,count in m.mapping_status.value_counts().items():mapping_summary.append(dict(version=v,entity='compound',scope='union_of_graph_used_local_ids',status=status,count=int(count),denominator=len(m)))
  active=m[m.chembl_id.isin(set(frames[v].parent_id))]
  for status,count in active.mapping_status.value_counts().items():mapping_summary.append(dict(version=v,entity='compound',scope='release_graph_used_local_ids',status=status,count=int(count),denominator=len(active)))
  t=target_summary(v).merge(validt[['local_id','signature']].rename(columns={'local_id':'reference_id'}),on='signature',how='left',validate='many_to_one')
  t['mapping_status']=np.where(t.eligible & t.reference_id.notna(),'HIGH_CONFIDENCE_ACCESSION_AND_SEQUENCE','UNMAPPED_OR_NONUNIQUE_COMPONENT')
  t.loc[t.eligible & t.signature.isin(duplicate_signatures),'mapping_status']='AMBIGUOUS_ONE_TO_MANY'
  dup=t.loc[t.mapping_status.eq('HIGH_CONFIDENCE_ACCESSION_AND_SEQUENCE'),'reference_id'].duplicated(keep=False)
  bad=set(t.loc[t.mapping_status.eq('HIGH_CONFIDENCE_ACCESSION_AND_SEQUENCE')].loc[dup,'local_id'])
  t.loc[t.local_id.isin(bad),'mapping_status']='AMBIGUOUS_MANY_TO_ONE'
  t.to_parquet(V/'derived'/f'target_mapping_{v}.parquet',index=False)
  tmaps[v]=t.loc[t.mapping_status.eq('HIGH_CONFIDENCE_ACCESSION_AND_SEQUENCE')].set_index('local_id').reference_id.to_dict()
  for status,count in t[t.local_id.isin(set(frames[v].target_id))].mapping_status.value_counts().items():mapping_summary.append(dict(version=v,entity='target',scope='release_graph_used_local_ids',status=status,count=int(count),denominator=t.local_id.isin(set(frames[v].target_id)).sum()))
  common_c=set(cmaps[v].values()) if common_c is None else common_c & set(cmaps[v].values())
  common_t=set(tmaps[v].values()) if common_t is None else common_t & set(tmaps[v].values())
  pd.DataFrame(mapping_summary).to_csv(V/'tables/ENTITY_MAPPING_COVERAGE.tsv',sep='\t',index=False)
  del m,full,f,t
 stage('GLOBAL_STABLE_SUBSET')
 pd.DataFrame({'reference_parent_id':sorted(common_c)}).to_csv(V/'tables/STABLE_COMPOUND_IDS.tsv',sep='\t',index=False)
 pd.DataFrame({'reference_target_id':sorted(common_t)}).to_csv(V/'tables/STABLE_TARGET_IDS.tsv',sep='\t',index=False)
 def mapped(f,v):
  out=f[['parent_id','target_id']].copy();out.parent_id=out.parent_id.map(cmaps[v]);out.target_id=out.target_id.map(tmaps[v])
  return out[out.parent_id.isin(common_c)&out.target_id.isin(common_t)].drop_duplicates()
 mapped_frames={v:mapped(f,v) for v,f in frames.items()}
 membership=pd.concat([f.assign(bit=1<<j) for j,f in enumerate(mapped_frames.values())],ignore_index=True).groupby(['parent_id','target_id'],sort=False).bit.sum().reset_index()
 membership['first_index']=[next(j for j in range(8) if int(b)&(1<<j)) for b in membership.bit]
 membership.to_parquet(V/'derived/STABLE_MAPPED_RELEASE_MEMBERSHIP.parquet',index=False)
 times=np.array(list(cfg['releases'].values()),dtype='datetime64[D]');firstdate=times[membership.first_index.to_numpy()]
 results=[]
 baseline=pd.read_csv(Q/'tables/REBUILT_LABEL_SUMMARY.tsv',sep='\t').set_index('cutoff')
 for year,spec in cfg['cutoffs'].items():
  stage('STABLE_MAPPING_GRAPH_AND_LABEL_SENSITIVITY',year)
  v=spec['release'];s=mapped_frames[v]
  rpath=Q/'derived'/('edges_37_2018.parquet' if year=='2018' else f'rollback_37_{year}.parquet')
  rawretro=pd.read_parquet(rpath,columns=['parent_id','target_id']);t=mapped(rawretro,'37')
  sset=set(map(tuple,s.to_numpy()));tset=set(map(tuple,t.to_numpy()))
  stable_backfill=len(tset-sset)/len(tset)
  local_s=frames[v][frames[v].parent_id.map(cmaps[v]).isin(common_c)&frames[v].target_id.map(tmaps[v]).isin(common_t)]
  local_t=rawretro[rawretro.parent_id.map(cmaps['37']).isin(common_c)&rawretro.target_id.map(tmaps['37']).isin(common_t)]
  local_sset=set(map(tuple,local_s.to_numpy()));local_tset=set(map(tuple,local_t.to_numpy()))
  subset_original_backfill=len(local_tset-local_sset)/len(local_tset)
  event=(firstdate>np.datetime64(spec['cutoff']))&(firstdate<=np.datetime64(spec['followup_end']))
  positive_set=set(map(tuple,membership.loc[event,['parent_id','target_id']].to_numpy()))
  original_retro=set(map(tuple,rawretro.to_numpy()))
  n=oldpos=old_erased=newpos=new_erased=leak_changed=0
  for batch in pq.ParquetFile(Q/'derived'/f'REBUILT_U2_LABELS_{year}.parquet').iter_batches(columns=['parent_id','target_id','first_discovery'],batch_size=250000):
   f=batch.to_pandas();mappedp=f.parent_id.map(cmaps[v]);mappedt=f.target_id.map(tmaps[v]);keep=mappedp.isin(common_c)&mappedt.isin(common_t)
   f=f[keep];rp=mappedp[keep];rt=mappedt[keep];n+=len(f)
   pairs=list(zip(rp,rt));new=np.array([x in positive_set for x in pairs]);old=f.first_discovery.to_numpy(dtype=bool)
   oldpos+=int(old.sum());old_erased+=sum(bool(e) and pair in original_retro for e,pair in zip(old,f[['parent_id','target_id']].itertuples(index=False,name=None)))
   newpos+=int(new.sum());new_erased+=sum(bool(e) and pair in tset for e,pair in zip(new,pairs));leak_changed+=int((old!=new).sum())
  original_graph=set(map(tuple,frames[v][['parent_id','target_id']].to_numpy()))
  results.append(dict(cutoff=year,stable_compounds_all_eight=len(common_c),stable_targets_all_eight=len(common_t),strict_edges=len(sset),strict_edge_coverage=len(sset)/len(original_graph),rollback_edges=len(tset),rollback_edge_coverage=len(tset)/len(rawretro),strict_only=len(sset-tset),rollback_only=len(tset-sset),stable_backfill_fraction=stable_backfill,subset_original_backfill_fraction=subset_original_backfill,mapping_delta_backfill_fraction=stable_backfill-subset_original_backfill,subset_original_strict_only=len(local_sset-local_tset),original_backfill_fraction=len(original_retro-original_graph)/len(original_retro),retained_U2_candidates=n,candidate_coverage=n/int(baseline.loc[int(year),'candidates']),subset_original_future_positives=oldpos,subset_mapped_future_positives=newpos,subset_label_changes=leak_changed,subset_original_erased=old_erased,subset_mapped_erased=new_erased,subset_original_erasure_fraction=old_erased/oldpos if oldpos else np.nan,subset_mapped_erasure_fraction=new_erased/newpos if newpos else np.nan,mapping_delta_erasure_fraction=(new_erased/newpos-old_erased/oldpos) if newpos and oldpos else np.nan))
  pd.DataFrame(results).to_csv(V/'tables/STABLE_MAPPING_SENSITIVITY.tsv',sep='\t',index=False)
 stage('COMPLETE')
 p=V/'logs/entity_mapping_progress.json';d=json.loads(p.read_text());d['status']='COMPLETE';p.write_text(json.dumps(d,indent=2))
except Exception as e:
 p=V/'logs/entity_mapping_progress.json';d=json.loads(p.read_text()) if p.exists() else {};d.update(status='FAILED',error=repr(e));p.write_text(json.dumps(d,indent=2));raise
