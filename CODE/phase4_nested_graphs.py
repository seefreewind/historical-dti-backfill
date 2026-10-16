"""QC-04 pure addition/removal on a fixed common-clean test set, frozen training."""
from repair_common import *
import phase2_benchmark as p2
import pyarrow.parquet as pq
import gc,time

def main():
 frames={v:pd.read_parquet(configured('primary_dir')/f'primary_{v}.parquet') for v in CFG['releases']}
 compounds=pd.Index(pd.unique(np.concatenate([f.parent_id.to_numpy() for f in frames.values()])))
 targets=pd.Index(pd.unique(np.concatenate([f.target_id.to_numpy() for f in frames.values()])))
 nc,nt=len(compounds),len(targets)
 h=json.loads((ROOT/'11_phase2/HYPERPARAMETER_FREEZE.json').read_text());results=[];edges=[]
 for year,ver in ((2018,'24.1'),(2020,'27')):
  strict=frames[ver][['parent_id','target_id']].drop_duplicates()
  retro=pd.read_parquet(ROOT/f'03_historical_snapshots/CHEMBL_RETRO_{year}/edges.parquet',columns=['parent_id','target_id']).drop_duplicates()
  rindex=pd.MultiIndex.from_frame(retro)
  common=strict[pd.MultiIndex.from_frame(strict).isin(rindex)]
  union=pd.concat([strict,retro],ignore_index=True).drop_duplicates()
  test=pd.read_parquet(ROOT/f'11_phase2/COMMON_CLEAN_{year}.parquet')
  cc=compounds.get_indexer(test.parent_id).astype(np.int32);ct=targets.get_indexer(test.target_id).astype(np.int32)
  assert (cc>=0).all() and (ct>=0).all()
  y=test.Y_future_36m.to_numpy();parents=test.parent_id.to_numpy()
  for arm,graph in [('G0_STRICT',strict),('G1_STRICT_PLUS_ROLLBACK_ONLY',union),('G2_COMMON_EDGES',common),('G3_FULL_ROLLBACK',retro)]:
   t0=time.time();print('START',year,arm,len(graph),flush=True)
   c=compounds.get_indexer(graph.parent_id).astype(np.int32);t=targets.get_indexer(graph.target_id).astype(np.int32)
   dc=np.bincount(c,minlength=nc).astype(np.float32);dt=np.bincount(t,minlength=nt).astype(np.float32)
   score={'M3_product':dc[cc]*dt[ct]}
   if arm.startswith(('G0','G3')):
    prefix='strict_' if arm.startswith('G0') else 'retro_'
    for model in ('M7_svd','M9_lightgcn'):score[model]=test[prefix+model].to_numpy()
   else:
    u,s,vt=p2.fit_svd(c,t,nc,nt,h['svd_dim'],p2.SEED)
    # Chunk score construction to avoid a candidates x embedding memory spike.
    values=np.empty(len(cc),dtype=np.float32)
    for a in range(0,len(cc),100000):
     b=min(a+100000,len(cc));values[a:b]=np.sum(u[cc[a:b]]*s*vt[:,ct[a:b]].T,axis=1)
    score['M7_svd']=values;del u,s,vt;gc.collect()
    eu,ei=p2.train_lightgcn(c,t,nc,nt,h['lightgcn_dim'],h['epochs'],h['lightgcn_lr'],h['reg'],h['negatives'],p2.SEED,h['lightgcn_layers'])
    values=np.empty(len(cc),dtype=np.float32)
    for a in range(0,len(cc),100000):
     b=min(a+100000,len(cc));values[a:b]=np.sum(eu[cc[a:b]]*ei[ct[a:b]],axis=1)
    score['M9_lightgcn']=values;del eu,ei;gc.collect()
    frame=test[['parent_id','target_id','Y_future_36m']].copy()
    for model,values in score.items():frame[model]=values
    frame.to_parquet(OUT/f'NESTED_SCORES_{year}_{arm[:2]}.parquet',index=False)
   for model,values in score.items():
    results.append(dict(cutoff=year,arm=arm,model=model,**global_metrics(y,values),**compound_metrics(parents,y,values),
        expected_random_tie_fde_at_1000=global_metrics(y,values,k=1000)['expected_random_tie_fde'],
        training_source='FROZEN_COMMON_CLEAN_SCORES' if arm.startswith(('G0','G3')) else 'FROZEN_HYPERPARAMETERS_RETRAIN',seed=p2.SEED,n_nodes_compound=nc,n_nodes_target=nt,n_graph_edges=len(graph)))
   edges.append(dict(cutoff=year,arm=arm,n_edges=len(graph),n_test=len(test),n_positive_test=int(y.sum()),training_seconds=time.time()-t0))
   save(results,'NESTED_GRAPH_ABLATION.tsv');save(edges,'NESTED_GRAPH_ARM_MANIFEST.tsv')
   print('DONE',year,arm,'seconds',round(time.time()-t0,1),flush=True)
  del test;gc.collect()
 df=pd.DataFrame(results);deltas=[]
 for (year,model),g in df.groupby(['cutoff','model']):
  base=g[g.arm.eq('G0_STRICT')].iloc[0]
  for effect,arm in [('ADD','G1_STRICT_PLUS_ROLLBACK_ONLY'),('REMOVE','G2_COMMON_EDGES'),('FULL','G3_FULL_ROLLBACK')]:
   row=g[g.arm.eq(arm)].iloc[0]
   deltas.append(dict(cutoff=year,model=model,effect=effect,**{'delta_'+metric:row[metric]-base[metric] for metric in ('full_boundary_tie_fde','expected_random_tie_fde','expected_recall_at_1','mrr_tie_expectation','expected_random_tie_fde_at_1000')}))
 save(deltas,'NESTED_GRAPH_EFFECTS.tsv')
def augment_top1000():
 table=pd.read_csv(OUT/'NESTED_GRAPH_ABLATION.tsv',sep='\t')
 for year in (2018,2020):
  frozen=pd.read_parquet(ROOT/f'11_phase2/COMMON_CLEAN_{year}.parquet')
  for arm in table.arm.unique():
   if arm.startswith(('G0','G3')):
    frame=frozen;prefix='strict_' if arm.startswith('G0') else 'retro_'
   else:
    frame=pd.read_parquet(OUT/f'NESTED_SCORES_{year}_{arm[:2]}.parquet');prefix=''
   for model in ('M3_product','M7_svd','M9_lightgcn'):
    g=global_metrics(frame.Y_future_36m,frame[prefix+model],k=1000)
    ix=table.cutoff.eq(year)&table.arm.eq(arm)&table.model.eq(model)
    for key in ('full_boundary_tie_fde','actual_evaluated_n','actual_evaluated_fraction','expected_random_tie_fde','boundary_score','n_boundary_tie'):
     table.loc[ix,key+'_at_1000']=g[key]
 save(table,'NESTED_GRAPH_ABLATION.tsv')
 effects=pd.read_csv(OUT/'NESTED_GRAPH_EFFECTS.tsv',sep='\t')
 for i,r in effects.iterrows():
  q=table[table.cutoff.eq(r.cutoff)&table.model.eq(r.model)]
  arm={'ADD':'G1_STRICT_PLUS_ROLLBACK_ONLY','REMOVE':'G2_COMMON_EDGES','FULL':'G3_FULL_ROLLBACK'}[r.effect]
  effects.loc[i,'delta_full_boundary_tie_fde_at_1000']=q[q.arm.eq(arm)].full_boundary_tie_fde_at_1000.iloc[0]-q[q.arm.eq('G0_STRICT')].full_boundary_tie_fde_at_1000.iloc[0]
 save(effects,'NESTED_GRAPH_EFFECTS.tsv')
if __name__=='__main__':
 import sys
 if '--augment-top1000' not in sys.argv:main()
 augment_top1000()
