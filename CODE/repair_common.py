"""Shared Phase 4 definitions. Input paths are project relative and configurable."""
from pathlib import Path
import os,json,math
import numpy as np
import pandas as pd
ROOT=Path(os.environ.get('HISTORICAL_DTI_ROOT',Path(__file__).resolve().parents[1]))
CONFIG_PATH=Path(os.environ.get('HISTORICAL_DTI_CONFIG',Path(__file__).resolve().parents[1]/'CONFIG/repair_config.json'))
CFG=json.loads(CONFIG_PATH.read_text())
OUT=ROOT/'14_repair'
def configured(key):
 p=Path(os.environ.get('HISTORICAL_DTI_'+key.upper(),CFG[key])); return p if p.is_absolute() else ROOT/p

def save(rows,name):
 df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
 df.to_csv(OUT/name,sep='\t',index=False);return df

def expected_first_reciprocal(before,b,p):
 # Position of first positive in a uniformly permuted block with p positives.
 prob=p/b; result=0.
 for j in range(1,b-p+2):
  result+=prob/(before+j)
  if j<=b-p:prob*= (b-j-p+1)/(b-j)
 return result

def compound_metrics(parent,y,score):
 y=np.asarray(y,dtype=bool);s=np.asarray(score,dtype=float)
 if len(y)==0 or not y.any():
  return dict(n_positive_compounds=0,old_metric=np.nan,old_phase2_mean_recall_at_1=np.nan,expected_recall_at_1=np.nan,top_tie_any_hit_fraction=np.nan,mrr_tie_expectation=np.nan,n_top_tied_compounds=0,median_top_tie_size=np.nan)
 codes,_=pd.factorize(parent,sort=False)
 positive_codes=np.unique(codes[y]); keep=np.isin(codes,positive_codes)
 c=codes[keep];s=s[keep];yy=y[keep]
 order=np.lexsort((-s,c)); c=c[order];s=s[order];yy=yy[order]
 cuts=np.r_[0,np.flatnonzero(np.diff(c))+1,len(c)]
 ex=[];anyhit=[];mrr=[];old=[];phase2=[];ties=[]
 for a,b in zip(cuts[:-1],cuts[1:]):
  ys=yy[a:b];ss=s[a:b];t=int(np.sum(ss==ss[0]));p=int(ys[:t].sum());ties.append(t)
  ex.append(p/t);anyhit.append(float(p>0));old.append(float(ys[0]));phase2.append(float(ys[0])/ys.sum())
  first=int(np.flatnonzero(ys)[0]);val=ss[first];before=int(np.sum(ss>val));block=ss==val
  mrr.append(expected_first_reciprocal(before,int(block.sum()),int(ys[block].sum())))
 return dict(n_positive_compounds=len(ex),old_metric=float(np.mean(old)) if ex else np.nan,
  old_phase2_mean_recall_at_1=float(np.mean(phase2)) if ex else np.nan,
  expected_recall_at_1=float(np.mean(ex)) if ex else np.nan,
  top_tie_any_hit_fraction=float(np.mean(anyhit)) if ex else np.nan,
  mrr_tie_expectation=float(np.mean(mrr)) if ex else np.nan,
  n_top_tied_compounds=int(np.sum(np.array(ties)>1)),median_top_tie_size=float(np.median(ties)) if ties else np.nan)

def global_metrics(y,s,fraction=.001,k=None):
 y=np.asarray(y,dtype=bool);s=np.asarray(s,dtype=float);n=len(y);p=int(y.sum())
 if n==0:return dict(n_candidates=0,n_positives=0,base_rate=np.nan,nominal_k=0,actual_evaluated_n=0,full_boundary_tie_fde=np.nan,expected_random_tie_fde=np.nan)
 k=min(n,max(1,round(n*fraction))) if k is None else min(k,n)
 th=np.partition(s,n-k)[n-k];above=s>th;tie=s==th;na=int(above.sum());nt=int(tie.sum());ha=int(y[above].sum());ht=int(y[tie].sum());slots=k-na
 factor=(n/p) if p else np.nan
 return dict(n_candidates=n,n_positives=p,base_rate=p/n,nominal_k=k,boundary_score=float(th),n_above=na,n_boundary_tie=nt,n_positive_boundary_tie=ht,
 actual_evaluated_n=na+nt,actual_evaluated_fraction=(na+nt)/n,
 full_boundary_tie_hits=ha+ht,full_boundary_tie_fde=(ha+ht)/(na+nt)*factor,
 expected_random_tie_hits=ha+slots*ht/nt,expected_random_tie_fde=(ha+slots*ht/nt)/k*factor,
 order_fde_min=(ha+max(0,slots-(nt-ht)))/k*factor,order_fde_max=(ha+min(slots,ht))/k*factor)
