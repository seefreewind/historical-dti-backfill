"""P2 controls. Global ranks and exact simple bipartite degree preservation."""
from repair_common import *
import phase2_benchmark as p2
import gc,sys,time

def gaps():
 rows=[]
 for year in (2018,2020,2022):
  labels=pd.read_parquet(OUT/f'AUTHORITATIVE_U2_LABELS_{year}.parquet',columns=['parent_id','target_id','first_discovery','gap_reappearance'])
  df=pd.read_parquet(ROOT/f'11_phase2/MODEL_SCORES_STRICT_{year}.parquet')
  assert labels[['parent_id','target_id']].equals(df[['parent_id','target_id']])
  masks={'FUTURE_POSITIVES':labels.first_discovery.eq(1).to_numpy(),'GAP_REAPPEARANCE':labels.gap_reappearance.to_numpy(),'ALL_CANDIDATES':np.ones(len(df),bool)}
  for model in [c for c in df if c.startswith('M')]:
   rank=df[model].rank(method='average',pct=True).to_numpy()
   for name,mask in masks.items():
    vals=rank[mask];qs=np.quantile(vals,[.25,.5,.75]) if len(vals) else [np.nan]*3
    rows.append(dict(cutoff=year,model=model,group=name,n=len(vals),mean=float(vals.mean()) if len(vals) else np.nan,median=qs[1],q25=qs[0],q75=qs[2],iqr=qs[2]-qs[0],reference='FULL_U2_AVERAGE_RANK_ONCE'))
  save(rows,'GAP_CONTROL_GLOBAL_PERCENTILES.tsv');print('GLOBAL_GAP_RANKS',year,flush=True)

def rewire():
 frames={v:pd.read_parquet(configured('primary_dir')/f'primary_{v}.parquet',columns=['parent_id','target_id']) for v in CFG['releases']}
 compounds=pd.Index(pd.unique(np.concatenate([f.parent_id.to_numpy() for f in frames.values()])))
 targets=pd.Index(pd.unique(np.concatenate([f.target_id.to_numpy() for f in frames.values()])))
 nc,nt=len(compounds),len(targets);h=json.loads((ROOT/'11_phase2/HYPERPARAMETER_FREEZE.json').read_text());audit=[];metrics=[]
 for year,v in ((2018,'24.1'),(2020,'27')):
  f=frames[v];c=compounds.get_indexer(f.parent_id).astype(np.int32);t=targets.get_indexer(f.target_id).astype(np.int32)
  basec=np.bincount(c,minlength=nc);baset=np.bincount(t,minlength=nt)
  def record(method,cc,tt,dup,swaps,attempts):
   dc=np.bincount(cc,minlength=nc)-basec;dt=np.bincount(tt,minlength=nt)-baset
   audit.append(dict(cutoff=year,method=method,n_original_edges=len(c),n_rewired_edges=len(cc),compound_changed_nodes=int(np.count_nonzero(dc)),target_changed_nodes=int(np.count_nonzero(dt)),compound_max_absolute_degree_difference=int(abs(dc).max()),target_max_absolute_degree_difference=int(abs(dt).max()),compound_mean_absolute_degree_difference=float(abs(dc).mean()),target_mean_absolute_degree_difference=float(abs(dt).mean()),number_duplicate_edges_removed=dup,successful_swaps=swaps,attempts=attempts,seed=20261015,note='Legacy algorithm reenactment, not recovery of unsaved frozen RNG realization' if method=='LEGACY_SHUFFLE_DEDUP_REENACTED' else 'Simple graph edge swaps; exact degrees verified'))
  rng=np.random.default_rng(20261015);oldt=t.copy();rng.shuffle(oldt);keep=np.unique(p2.pack(c,oldt,nt),return_index=True)[1]
  record('LEGACY_SHUFFLE_DEDUP_REENACTED',c[keep],oldt[keep],len(c)-len(keep),0,1)
  rng=np.random.default_rng(20261015);newt=t.copy();edge_set=set(p2.pack(c,t,nt).tolist());done=0;attempts=0;goal=5*len(c)
  while done<goal and attempts<50*len(c):
   # Candidate edge positions are independent of outcome labels.
   pairs=rng.integers(0,len(c),size=(10000,2))
   for i,j in pairs:
    attempts+=1
    a,b=int(c[i]),int(c[j]);x,y=int(newt[i]),int(newt[j])
    if a==b or x==y:continue
    k1=a*nt+y;k2=b*nt+x
    if k1 in edge_set or k2 in edge_set:continue
    edge_set.remove(a*nt+x);edge_set.remove(b*nt+y);edge_set.add(k1);edge_set.add(k2)
    newt[i]=y;newt[j]=x;done+=1
    if done>=goal:break
  assert len(edge_set)==len(c) and np.array_equal(np.bincount(c,minlength=nc),basec) and np.array_equal(np.bincount(newt,minlength=nt),baset)
  record('EXACT_SIMPLE_BIPARTITE_DOUBLE_EDGE_SWAP',c,newt,0,done,attempts)
  save(audit,'REWIRE_DEGREE_AUDIT.tsv')
  pd.DataFrame({'parent_id':compounds.to_numpy()[c],'target_id':targets.to_numpy()[newt]}).to_parquet(OUT/f'REWIRED_EDGES_{year}.parquet',index=False)
  eu,ei=p2.train_lightgcn(c,newt,nc,nt,h['lightgcn_dim'],h['epochs'],h['lightgcn_lr'],h['reg'],h['negatives'],20261015,h['lightgcn_layers'])
  df=pd.read_parquet(ROOT/f'11_phase2/MODEL_SCORES_STRICT_{year}.parquet',columns=['parent_id','target_id','Y_future_36m','M9_lightgcn'])
  cc=compounds.get_indexer(df.parent_id);tt=targets.get_indexer(df.target_id);s=np.empty(len(df),np.float32)
  for a in range(0,len(df),100000):
   b=min(a+100000,len(df));s[a:b]=np.sum(eu[cc[a:b]]*ei[tt[a:b]],axis=1)
  for condition,score in [('STRICT',df.M9_lightgcn.to_numpy()),('EXACT_REWIRE',s)]:metrics.append(dict(cutoff=year,model='M9_lightgcn',condition=condition,**global_metrics(df.Y_future_36m.to_numpy(),score)))
  save(metrics,'REWIRE_CONTROL_METRICS.tsv');print('REWIRE',year,done,flush=True)
  del eu,ei,df;gc.collect()
if __name__=='__main__':{'gaps':gaps,'rewire':rewire}[sys.argv[1]]()
