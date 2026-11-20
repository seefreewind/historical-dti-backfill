"""Authoritative first-seen labels from materialized release edge sequence.

No Phase 1 summary or old outcome label enters construction. Parent/target IDs
are the release-local frozen identifiers, with no later canonical remapping.
"""
from repair_common import *
import pyarrow.parquet as pq
from datetime import date
import re

def build():
 versions=list(CFG['releases']);dates=np.array([CFG['releases'][v] for v in versions],dtype='datetime64[D]')
 frames={v:pd.read_parquet(configured('primary_dir')/f'primary_{v}.parquet',columns=['parent_id','target_id']).drop_duplicates() for v in versions}
 compounds=pd.Index(pd.unique(np.concatenate([x.parent_id.to_numpy() for x in frames.values()])))
 targets=pd.Index(pd.unique(np.concatenate([x.target_id.to_numpy() for x in frames.values()])))
 nt=len(targets);keyframes=[]
 for j,(v,f) in enumerate(frames.items()):
  k=compounds.get_indexer(f.parent_id).astype(np.int64)*nt+targets.get_indexer(f.target_id)
  keyframes.append(pd.DataFrame({'key':k,'membership':np.uint16(1<<j)}))
 union=pd.concat(keyframes,ignore_index=True).groupby('key',sort=True).membership.sum().astype(np.uint16)
 keys=union.index.to_numpy();bits=union.to_numpy();first=np.full(len(keys),-1,np.int16)
 for j in range(len(versions)):first[(first<0)&((bits&(1<<j))!=0)]=j
 membership=np.array(['|'.join(v for j,v in enumerate(versions) if int(x)&(1<<j)) for x in np.unique(bits)],dtype=object)
 membership_map=dict(zip(np.unique(bits),membership))
 master=pd.DataFrame({'parent_id':compounds.to_numpy()[keys//nt],'target_id':targets.to_numpy()[keys%nt],
 'release_membership_bitmask':bits,'all_release_membership':[membership_map[x] for x in bits],
 'first_public_release_present':[versions[j] for j in first],'first_operational_date':dates[first].astype(str)})
 master.to_parquet(OUT/'AUTHORITATIVE_RELEASE_MEMBERSHIP.parquet',index=False)
 migrations=[];summary=[]
 def label(frame,year):
  spec=CFG['cutoffs'][str(year)];cut=np.datetime64(spec['cutoff']);end=np.datetime64(spec['followup_end'])
  c=compounds.get_indexer(frame.parent_id);t=targets.get_indexer(frame.target_id);k=c.astype(np.int64)*nt+t
  loc=np.searchsorted(keys,k);safe=np.minimum(loc,len(keys)-1);valid=(c>=0)&(t>=0)&(loc<len(keys))&(keys[safe]==k)
  b=np.where(valid,bits[safe],0).astype(np.uint16);f=np.where(valid,first[safe],-1)
  before_mask=sum(1<<j for j,d in enumerate(dates) if d<=cut)
  strict_index=versions.index(spec['release']);presentbefore=(b&before_mask)!=0;at=(b&(1<<strict_index))!=0
  after=np.full(len(frame),-1,np.int16)
  for j,d in enumerate(dates):
   if d>cut:after[(after<0)&((b&(1<<j))!=0)]=j
  firstdate=np.full(len(frame),np.datetime64('NaT'),dtype='datetime64[D]');firstdate[valid]=dates[f[valid]]
  afterdate=np.full(len(frame),np.datetime64('NaT'),dtype='datetime64[D]');afterdate[after>=0]=dates[after[after>=0]]
  event=valid&~presentbefore&(firstdate>cut)&(firstdate<=end)
  gap=presentbefore&~at&(after>=0)&(afterdate<=end)
  result=frame[['parent_id','target_id']].copy();result['cutoff']=year
  result['first_public_release_present']=np.where(valid,np.array(versions,dtype=object)[np.maximum(f,0)],'NOT_OBSERVED')
  result['all_release_membership']=[membership_map.get(x,'NONE') for x in b]
  result['release_membership_bitmask']=b
  result['present_before_cutoff']=presentbefore;result['present_at_cutoff']=at
  result['first_present_after_cutoff']=np.where(after>=0,np.array(versions,dtype=object)[np.maximum(after,0)],'NOT_OBSERVED')
  result['gap_reappearance']=gap;result['first_discovery']=event.astype(np.int8);result['event_observed']=event
  result['event_release_date']=np.where(event,firstdate.astype(str),'')
  result['days_to_event']=np.where(event,(firstdate-cut).astype('timedelta64[D]').astype(float),np.nan)
  result['followup_end']=str(end);result['right_censored']=~event
  result['later_first_appearance_after_followup']=valid&~presentbefore&(firstdate>end)
  result['censor_status']=np.where(event,'EVENT_OBSERVED',np.where(result.later_first_appearance_after_followup,'FIRST_SEEN_AFTER_FOLLOWUP','NO_FIRST_DISCOVERY_OBSERVED'))
  result['days_observed']=np.where(event,(firstdate-cut).astype('timedelta64[D]').astype(float),float((end-cut).astype(int)))
  return result
 def compare(old,new,year,name,requires_first_seen=True):
  oy=old.Y_future_36m.to_numpy().astype(bool);ny=new.first_discovery.to_numpy().astype(bool)
  migrations.append(dict(cutoff=year,source=name,requires_first_seen_reconciliation=requires_first_seen,label_definition='FIRST_SEEN' if requires_first_seen else 'LEGACY_DOCUMENT_YEAR_ENDPOINT',n_pairs=len(old),old_positive_new_positive=int((oy&ny).sum()),old_positive_new_negative=int((oy&~ny).sum()),old_negative_new_positive=int((~oy&ny).sum()),old_negative_new_negative=int((~oy&~ny).sum()),gap_corrected=int((new.gap_reappearance&~ny).sum()),censor_changed=0,mapping_changed=0,censor_note='EXPLICIT_METADATA_ADDED; NO_OLD_CENSOR_MASK_TO_COMPARE',mapping_note='SAME_FROZEN_RELEASE_LOCAL_KEYS'))
  diff=new[oy!=ny].copy()
  if len(diff):diff['old_Y']=oy[oy!=ny].astype(int);diff['source']=name;diff.to_csv(OUT/f'LABEL_DIFFERENCES_{year}_{name}.tsv',sep='\t',index=False)
 for year in (2018,2020,2022):
  old=pd.read_parquet(ROOT/f'10_phase1/DRUG_TARGET_LABELS_{year}.parquet');new=label(old,year)
  compare(old,new,year,'PHASE1_CORRECTED');new.to_parquet(OUT/f'LABELS_PHASE1_DOMAIN_{year}.parquet',index=False)
  strict=frames[CFG['cutoffs'][str(year)]['release']]
  possible=master[master.parent_id.isin(strict.parent_id.unique())&master.target_id.isin(strict.target_id.unique())].copy()
  broad=label(possible,year);broad=broad[~broad.present_at_cutoff];broad.to_parquet(OUT/f'LABELS_HISTORICAL_ACTIVE_OBSERVED_{year}.parquet',index=False)
  positive=broad[broad.first_discovery.eq(1)]
  oldpos=set(map(tuple,old.loc[old.Y_future_36m.eq(1),['parent_id','target_id']].to_numpy()));newpos=set(map(tuple,positive[['parent_id','target_id']].to_numpy()))
  migrations.append(dict(cutoff=year,source='HISTORICAL_ACTIVE_COMPLETE_POSITIVE_SET',requires_first_seen_reconciliation=True,label_definition='FIRST_SEEN',n_pairs=len(broad),old_positive_new_positive=len(oldpos&newpos),old_positive_new_negative=len(oldpos-newpos),old_negative_new_positive=len(newpos-oldpos),gap_corrected=int(broad.gap_reappearance.sum()),mapping_changed=0,censor_changed=0))
  if oldpos!=newpos:save([dict(cutoff=year,parent_id=p,target_id=t,change='NEW_ONLY' if (p,t) in newpos else 'OLD_ONLY') for p,t in oldpos^newpos],f'LABEL_DIFFERENCES_BROAD_{year}.tsv')
  score_ids=pd.read_parquet(ROOT/f'11_phase2/MODEL_SCORES_STRICT_{year}.parquet',columns=['parent_id','target_id','Y_future_36m'])
  n=label(score_ids,year);compare(score_ids,n,year,'PHASE2_FULL_U2');n.to_parquet(OUT/f'AUTHORITATIVE_U2_LABELS_{year}.parquet',index=False)
  for prefix in ('STRICT_BENCHMARK','COMMON_CLEAN','MODEL_SCORES_RETRO'):
   additional=pd.read_parquet(ROOT/f'11_phase2/{prefix}_{year}.parquet',columns=['parent_id','target_id','Y_future_36m'])
   compare(additional,label(additional,year),year,'PHASE2_'+prefix,requires_first_seen=prefix!='MODEL_SCORES_RETRO')
  table=(ROOT/'13_manuscript/MAIN_MANUSCRIPT_V1.md').read_text()
  tableline=next(x for x in table.splitlines() if x.startswith(f'| {year}-12-31 |'))
  cells=[x.strip() for x in tableline.strip('|').split('|')]
  manuscript_positive=int(re.sub(r'[^0-9]','',cells[7]))
  manuscript_erased=int(cells[9].split('/')[0].replace(',',''))
  # Phase 3 identity, first-positive strata, and unscored panels are migration inputs only.
  for filename in ('CHEMICAL_IDENTITY_AUDIT.tsv','CHEMICAL_NOVELTY_STRATA.tsv','LLM_PAIR_PANEL.tsv','LLM_COMPOUND_RANKING_PANEL.tsv'):
   p=pd.read_csv(ROOT/'12_phase3'/filename,sep='\t');p=p[p.cutoff.eq(year)]
   if len(p):compare(p,label(p,year),year,'PHASE3_'+filename.removesuffix('.tsv'))
  retro=pd.read_parquet(ROOT/f'03_historical_snapshots/CHEMBL_RETRO_{year}/edges.parquet',columns=['parent_id','target_id'])
  ri=pd.MultiIndex.from_frame(retro).unique();erased=int(n.first_discovery.to_numpy()[pd.MultiIndex.from_frame(n[['parent_id','target_id']]).isin(ri)].sum())
  broad_erased=int(pd.MultiIndex.from_frame(positive[['parent_id','target_id']]).isin(ri).sum())
  expected={2018:(946,633,5611,4508),2020:(1413,932,6775,5287),2022:(2104,1594,8924,6750)}[year]
  observed=(int(n.first_discovery.sum()),erased,len(positive),broad_erased)
  summary.append(dict(cutoff=year,u2_positive=observed[0],u2_erased=erased,historical_active_positive=len(positive),historical_active_erased=broad_erased,expected_frozen=repr(expected),observed=repr(observed),frozen_reconciled=observed==expected,manuscript_u2_positive=manuscript_positive,manuscript_u2_erased=manuscript_erased,manuscript_counts_reconciled=(observed[0]==manuscript_positive and observed[1]==manuscript_erased),followup_end=CFG['cutoffs'][str(year)]['followup_end']))
  save(migrations,'LABEL_MIGRATION_AUDIT.tsv');save(summary,'AUTHORITATIVE_LABEL_SUMMARY.tsv');print('LABEL',year,observed,flush=True)
def reconcile_migration():
 migration=pd.read_csv(OUT/'LABEL_MIGRATION_AUDIT.tsv',sep='\t')
 summary=pd.read_csv(OUT/'AUTHORITATIVE_LABEL_SUMMARY.tsv',sep='\t')
 master=pd.read_parquet(OUT/'AUTHORITATIVE_RELEASE_MEMBERSHIP.parquet')
 evidence=[]
 for year,end_version in ((2018,'29'),(2020,'33'),(2022,'36')):
  old=pd.read_parquet(ROOT/f'10_phase1/DRUG_TARGET_LABELS_{year}.parquet',columns=['parent_id','target_id','Y_future_36m'])
  broad=pd.read_parquet(OUT/f'LABELS_HISTORICAL_ACTIVE_OBSERVED_{year}.parquet')
  positive=broad[broad.first_discovery.eq(1)]
  at_endpoint=positive.all_release_membership.str.split('|',regex=False).apply(lambda x:end_version in x)
  endpoint_set=set(map(tuple,positive.loc[at_endpoint,['parent_id','target_id']].to_numpy()))
  old_set=set(map(tuple,old.loc[old.Y_future_36m.eq(1),['parent_id','target_id']].to_numpy()))
  endpoint_reconstructs_frozen=endpoint_set==old_set
  year_pass=endpoint_reconstructs_frozen
  for i,r in migration[migration.cutoff.eq(year)].iterrows():
   if not r.requires_first_seen_reconciliation:
    migration.loc[i,'reconciliation_status']='EXPECTED_DIFFERENT_ENDPOINT'
    migration.loc[i,'difference_explanation']='Legacy document-year endpoint, not a first-seen outcome'
    continue
   if r.old_positive_new_negative==0 and r.old_negative_new_positive==0:
    migration.loc[i,'reconciliation_status']='PASS_UNCHANGED'
    continue
   if r.source=='HISTORICAL_ACTIVE_COMPLETE_POSITIVE_SET':
    newonly=positive[~pd.MultiIndex.from_frame(positive[['parent_id','target_id']]).isin(pd.MultiIndex.from_tuples(list(old_set)))]
    d=newonly[['parent_id','target_id','first_public_release_present','all_release_membership']].copy()
   else:
    d=pd.read_csv(OUT/f'LABEL_DIFFERENCES_{year}_{r.source}.tsv',sep='\t')[['parent_id','target_id']].merge(master,on=['parent_id','target_id'],validate='one_to_one')
   absent_endpoint=d.all_release_membership.str.split('|',regex=False).apply(lambda x:end_version not in x)
   supported=endpoint_reconstructs_frozen and r.old_positive_new_negative==0 and len(d)==r.old_negative_new_positive and absent_endpoint.all()
   migration.loc[i,'reconciliation_status']='PASS_EXPLAINED_CORRECTION' if supported else 'FAIL_UNEXPLAINED'
   migration.loc[i,'difference_explanation']='First appeared within follow-up, absent at old endpoint release; old endpoint-conditioned construction omitted the event'
   year_pass=year_pass and supported
   d['cutoff']=year;d['source']=r.source;d['old_endpoint_release']=end_version;d['absent_old_endpoint']=absent_endpoint;d['evidence_status']='DIRECT_RELEASE_MEMBERSHIP'
   evidence.extend(d.to_dict('records'))
  strict=migration[migration.cutoff.eq(year)&migration.source.isin(['PHASE2_FULL_U2','PHASE2_STRICT_BENCHMARK','PHASE2_COMMON_CLEAN'])]
  year_pass=year_pass and strict.old_positive_new_negative.eq(0).all() and strict.old_negative_new_positive.eq(0).all()
  summary.loc[summary.cutoff.eq(year),'migration_reconciled']=bool(year_pass)
  summary.loc[summary.cutoff.eq(year),'legacy_positive_set_reconstructed_by_endpoint']=bool(endpoint_reconstructs_frozen)
 save(migration,'LABEL_MIGRATION_AUDIT.tsv');save(summary,'AUTHORITATIVE_LABEL_SUMMARY.tsv');save(evidence,'LABEL_RECONCILIATION_EVIDENCE.tsv')
 print('MIGRATION_RECONCILED',summary.migration_reconciled.tolist(),flush=True)
if __name__=='__main__':
 import sys
 if '--reconcile-only' not in sys.argv:build()
 reconcile_migration()
