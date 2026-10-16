"""QC-01, QC-02, QC-03 in priority order. Never overwrite frozen results."""
from repair_common import *
import sys,base64
from functools import lru_cache
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
VERS={2018:'24.1',2020:'27',2022:'31'}
YEARS=(2018,2020,2022)
def scores(year):return pd.read_parquet(ROOT/f'11_phase2/MODEL_SCORES_STRICT_{year}.parquet')
def model_columns(df):return [c for c in df if c.startswith('M')]
def scaffold_data(year,df):
 desc=pd.read_parquet(configured('descriptor_dir')/f'desc_{VERS[year]}.parquet')
 smiles=pd.read_parquet(configured('structure_dir')/f'smiles_{VERS[year]}.parquet')
 smap=smiles.set_index('molecule_chembl_id').canonical_smiles
 # Rebuild from actual molecules, not the old cached empty-to-ACYCLIC mapping.
 rows=[]
 for mid in pd.unique(df.parent_id):
  raw=smap.get(mid);mol=Chem.MolFromSmiles(raw) if isinstance(raw,str) else None
  if mol is None:token=None;reason='CANDIDATE_STRUCTURE_UNAVAILABLE_OR_INVALID'
  else:
   token=MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
   token=token if token else 'ACYCLIC';reason='VALID_STRUCTURE'
  rows.append((mid,token,reason))
 # Include all historical ligands, including those without candidate rows.
 have={x[0] for x in rows}
 edges=pd.read_parquet(ROOT/f'03_historical_snapshots/CHEMBL_STRICT_{year}/edges.parquet',columns=['parent_id','target_id'])
 for mid in pd.unique(edges.parent_id):
  if mid in have:continue
  raw=smap.get(mid);mol=Chem.MolFromSmiles(raw) if isinstance(raw,str) else None
  token=MurckoScaffold.MurckoScaffoldSmiles(mol=mol) if mol is not None else None
  if mol is not None and not token:token='ACYCLIC'
  rows.append((mid,token,'VALID_STRUCTURE' if mol is not None else 'STRUCTURE_UNAVAILABLE_OR_INVALID'))
 tokens=pd.DataFrame(rows,columns=['parent_id','scaffold_token','structure_status']).set_index('parent_id')
 tokens.reset_index().to_parquet(OUT/f'SCAFFOLD_TOKENS_{year}.parquet',index=False)
 lig=edges.join(tokens,on='parent_id');valid=lig[lig.scaffold_token.notna()]
 valid_pairs=pd.MultiIndex.from_frame(valid[['target_id','scaffold_token']].drop_duplicates())
 ct=df.parent_id.map(tokens.scaffold_token);has=df.target_id.isin(valid.target_id.unique())
 good=ct.notna() & has
 same=pd.MultiIndex.from_arrays([df.target_id,ct]).isin(valid_pairs)&good.to_numpy()
 state=np.where(~good,'UNKNOWN',np.where(same,'SAME','NEW'))
 reason=np.where(ct.isna(),'CANDIDATE_STRUCTURE_UNAVAILABLE_OR_INVALID',np.where(~has,'NO_VALID_HISTORICAL_LIGAND_SCAFFOLD','VALID_COMPARISON'))
 return state,ct,reason

def run_ties():
 rows=[]
 for year in YEARS:
  df=scores(year);y=df.Y_future_36m.to_numpy();parents=df.parent_id.to_numpy()
  for m in model_columns(df):rows.append(dict(cutoff=year,model=m,subset='ALL_U2',**compound_metrics(parents,y,df[m].to_numpy())))
  save(rows,'COMPOUND_CENTRIC_TIE_REPAIR.tsv');print('QC01',year,flush=True)

def run_scaffolds():
 allmetrics=[];toprows=[];classcounts=[]
 ties=pd.read_csv(OUT/'COMPOUND_CENTRIC_TIE_REPAIR.tsv',sep='\t').to_dict('records')
 for year in YEARS:
  df=scores(year)
  existing=OUT/f'SCAFFOLD_CLASSIFICATION_{year}.parquet'
  if existing.exists():
   cached=pd.read_parquet(existing)
   assert cached[['parent_id','target_id']].equals(df[['parent_id','target_id']])
   state=cached.scaffold_state.to_numpy();tokens=cached.candidate_scaffold;reasons=cached.reason
  else:state,tokens,reasons=scaffold_data(year,df)
  cls=df[['parent_id','target_id','Y_future_36m']].copy();cls['cutoff']=year;cls['scaffold_state']=state;cls['candidate_scaffold']=tokens;cls['reason']=reasons
  cls.to_parquet(OUT/f'SCAFFOLD_CLASSIFICATION_{year}.parquet',index=False)
  # Full row-level TSV is explicitly requested; append every candidate, not just positives.
  cls.to_csv(OUT/'SCAFFOLD_CLASSIFICATION_REPAIRED.tsv',sep='\t',index=False,mode='w' if year==2018 else 'a',header=year==2018)
  y=df.Y_future_36m.to_numpy();par=df.parent_id.to_numpy()
  masks={'ALL_U2':np.ones(len(df),bool),'SAME':state=='SAME','NEW':state=='NEW','UNKNOWN':state=='UNKNOWN',
    'TANIMOTO_LT_0.95':df.M5_ligand_max.to_numpy()<.95,'TANIMOTO_LT_0.80':df.M5_ligand_max.to_numpy()<.80}
  for sub,mask in masks.items():
   classcounts.append(dict(cutoff=year,subset=sub,n_candidates=int(mask.sum()),n_positives=int(y[mask].sum())))
   for m in model_columns(df):
    g=global_metrics(y[mask],df[m].to_numpy()[mask]);cm=compound_metrics(par[mask],y[mask],df[m].to_numpy()[mask])
    row=dict(cutoff=year,subset=sub,model=m,**g,**cm)
    toprows.append(dict(cutoff=year,subset=sub,model=m,**g))
    if sub in ('SAME','NEW','UNKNOWN'):allmetrics.append(row)
    if sub!='ALL_U2':ties.append(dict(cutoff=year,subset=sub,model=m,**cm))
  save(allmetrics,'SCAFFOLD_METRICS_REPAIRED.tsv');save(toprows,'TOP_PERCENTILE_TIE_AUDIT.tsv');save(ties,'COMPOUND_CENTRIC_TIE_REPAIR.tsv');save(classcounts,'SCAFFOLD_STATE_COUNTS.tsv')
  print('QC02',year,flush=True)

TAUT=rdMolStandardize.TautomerEnumerator()
@lru_cache(maxsize=60000)
def molecular(raw):
 mol=Chem.MolFromSmiles(raw) if isinstance(raw,str) else None
 if mol is None:return dict(valid=False)
 Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
 centres=Chem.FindMolChiralCenters(mol,includeUnassigned=True,includeCIP=True)
 info=Chem.FindPotentialStereo(mol)
 undef=sum(str(x.specified)=='Unspecified' for x in info)
 defined=sum(str(x.specified)=='Specified' for x in info)
 iso=Chem.MolToSmiles(mol,isomericSmiles=True,canonical=True)
 no=Chem.MolToSmiles(mol,isomericSmiles=False,canonical=True)
 ik=inchi.MolToInchiKey(mol);ii=inchi.MolToInchi(mol)
 isotopes=sorted((a.GetSymbol(),a.GetIsotope()) for a in mol.GetAtoms() if a.GetIsotope())
 return dict(valid=True,isomeric_smiles=iso,nonisomeric_smiles=no,inchikey=ik,connectivity=ik[:14],inchi=ii,
  chiral_tags=sum(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()),
  assigned_centres=sum(v!='?' for _,v in centres),undefined_centres=sum(v=='?' for _,v in centres),
  defined_stereo=defined,undefined_stereo=undef,
  defined_double_bonds=sum(b.GetStereo() not in (Chem.BondStereo.STEREONONE,Chem.BondStereo.STEREOANY) for b in mol.GetBonds()),
  isotopes=json.dumps(isotopes),formal_charge=Chem.GetFormalCharge(mol),
  tautomer_key=inchi.MolToInchiKey(TAUT.Canonicalize(mol)))

def relation(a,b):
 if not a['valid'] or not b['valid']:return 'G_UNRESOLVED','INVALID_STRUCTURE'
 if a['isomeric_smiles']==b['isomeric_smiles']:return 'A_EXACT_STRUCTURE','CANONICAL_ISOMERIC_SMILES_EQUAL'
 if a['connectivity']==b['connectivity']:
  if a['isotopes']!=b['isotopes'] or a['formal_charge']!=b['formal_charge']:return 'D_ISOTOPE_OR_OTHER_INCHI_LAYER_DIFFERENCE','ISOTOPE_OR_CHARGE_DIFFERENCE'
  if a['nonisomeric_smiles']==b['nonisomeric_smiles']:
   if a['defined_stereo']>0 and b['defined_stereo']>0 and a['undefined_stereo']==0 and b['undefined_stereo']==0:
    return 'B_DEFINED_STEREO_DIFFERENCE','BOTH_FULLY_SPECIFIED_STEREO_AND_DIFFERENT_ISOMERIC_SMILES'
   if a['defined_stereo']!=b['defined_stereo'] or a['undefined_stereo']!=b['undefined_stereo'] or min(a['defined_stereo'],b['defined_stereo'])==0:
    return 'C_STEREO_SPECIFICATION_DIFFERENCE','UNEQUAL_OR_ABSENT_STEREO_SPECIFICATION'
   return 'G_UNRESOLVED','BOTH_PARTIAL_STEREO_NO_PROVEN_DEFINED_PAIR'
  return 'E_STANDARDIZATION_OR_REPRESENTATION_DIFFERENCE','CONNECTIVITY_KEY_EQUAL_BUT_REPRESENTATION_DIFFERS'
 if a['tautomer_key']==b['tautomer_key']:return 'E_STANDARDIZATION_OR_REPRESENTATION_DIFFERENCE','CANONICAL_TAUTOMER_EQUAL'
 return 'F_FINGERPRINT_COLLISION','MORGAN_EQUAL_CONNECTIVITY_DIFFERENT'

def run_chemistry():
 old=pd.read_csv(ROOT/'12_phase3/CHEMICAL_IDENTITY_AUDIT.tsv',sep='\t')
 rows=[];sums=[]
 priority={c:i for i,c in enumerate('ABCDEFG')}
 for year in YEARS:
  desc=pd.read_parquet(configured('descriptor_dir')/f'desc_{VERS[year]}.parquet').set_index('molecule_chembl_id')
  ed=pd.read_parquet(ROOT/f'03_historical_snapshots/CHEMBL_STRICT_{year}/edges.parquet',columns=['parent_id','target_id'])
  by={}
  for mid,t in ed.itertuples(index=False,name=None):
   fp=desc.at[mid,'morgan'] if mid in desc.index else None
   if isinstance(fp,str):by.setdefault((t,fp),[]).append(mid)
  part=old[old.cutoff.eq(year)&old.m5_max_tanimoto.ge(1-1e-6)]
  for row in part.itertuples(index=False):
   mid=row.parent_id;fp=desc.at[mid,'morgan'];a=molecular(desc.at[mid,'raw_smiles']);matches=by.get((row.target_id,fp),[])
   calls=[]
   for lig in matches:
    b=molecular(desc.at[lig,'raw_smiles']);cl,why=relation(a,b);calls.append((priority[cl[0]],cl,why,lig,b))
   if calls:
    _,cl,why,lig,b=min(calls,key=lambda x:(x[0],x[3]))
   else:cl,why,lig,b='G_UNRESOLVED','NO_EQUAL_MORGAN_LIGAND_FOUND','',dict(valid=False)
   rows.append(dict(cutoff=year,parent_id=mid,target_id=row.target_id,Y_future_36m=row.Y_future_36m,old_class=row.relation_class,
     chemical_class=cl,classification_reason=why,matched_historical_ligand=lig,n_equal_fingerprint_ligands=len(matches),
     **{'candidate_'+k:v for k,v in a.items()},**{'ligand_'+k:v for k,v in b.items()}))
  frame=pd.DataFrame([r for r in rows if r['cutoff']==year])
  for cls in ('A_EXACT_STRUCTURE','B_DEFINED_STEREO_DIFFERENCE','C_STEREO_SPECIFICATION_DIFFERENCE','D_ISOTOPE_OR_OTHER_INCHI_LAYER_DIFFERENCE','E_STANDARDIZATION_OR_REPRESENTATION_DIFFERENCE','F_FINGERPRINT_COLLISION','G_UNRESOLVED'):
   q=frame[frame.chemical_class.eq(cls)]
   sums.append(dict(cutoff=year,chemical_class=cls,n_pairs=len(q),fraction=len(q)/len(frame),n_future_positives=int(q.Y_future_36m.sum()),denominator=len(frame)))
  save(rows,'TANIMOTO1_RECLASSIFIED.tsv');save(sums,'CHEMICAL_CLASSIFICATION_SUMMARY.tsv');print('QC03',year,flush=True)
if __name__=='__main__':
 mode=sys.argv[1]
 {'ties':run_ties,'scaffolds':run_scaffolds,'chemistry':run_chemistry}[mode]()
