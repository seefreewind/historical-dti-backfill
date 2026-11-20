"""Orchestrate byte-identical scientific functions in isolated outputs; no new algorithms."""
import os,sys,json,csv,hashlib,sqlite3,subprocess,datetime,time,traceback,fcntl,shutil
from pathlib import Path
sys.dont_write_bytecode=True
S=Path(os.environ['REPRODUCTION_RUN_ROOT']).resolve();W=Path(__file__).resolve().parent/'expected';O=S/'reproduction/work';Q=O/'jcheminf_phase1r';R=O/'biomedical_ai_historical_discovery';V=Q/'validity';HQ=W/'jcheminf_phase1r'
os.environ.update(MPLCONFIGDIR=str(S/'environment/matplotlib-cache'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
for _key in list(os.environ):
 if _key.startswith('HISTORICAL_DTI_'):os.environ.pop(_key)
os.environ.update(HISTORICAL_DTI_ROOT=str(R),HISTORICAL_DTI_CONFIG=str(R/'configs/repair_config.json'))
state={'status':'RUNNING','pid':os.getpid(),'started_at':datetime.datetime.now().astimezone().isoformat()};results=[]
lock=(S/'logs/clean_run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def stage(name,**extra):
 state.update(stage=name,updated_at=datetime.datetime.now().astimezone().isoformat(),**extra);(S/'logs/clean_progress.json').write_text(json.dumps(state,indent=2));print(state,flush=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def record(name,passed,expected='',observed='',**extra):
 results.append(dict(check=name,status='PASS' if passed else 'FAIL',expected=str(expected),observed=str(observed),**extra));(S/'validation/CHECKS.json').write_text(json.dumps(results,indent=2));assert passed,('RESULT_DISCREPANCY',name,expected,observed)
def run(script,tag,args=(),env=None):
 stage(tag)
 with (S/'logs'/f'{tag}.log').open('x') as f:subprocess.run([sys.executable,'-I','-B',str(script),*map(str,args)],cwd=O,stdout=f,stderr=subprocess.STDOUT,env=env or dict(os.environ),check=True)
def canonical_table(p,columns=None,keys=None):
 import pyarrow as pa,pyarrow.compute as pc,pyarrow.parquet as pq
 t=pq.read_table(p,columns=columns)
 for i in range(t.num_columns):
  a=t.column(i);typ=a.type
  if pa.types.is_dictionary(typ):a=pc.cast(a,typ.value_type);typ=a.type
  if pa.types.is_string(typ) or pa.types.is_large_string(typ):a=pc.cast(a,pa.large_string())
  t=t.set_column(i,t.schema.names[i],a)
 t=t.replace_schema_metadata(None)
 if keys:t=t.take(pc.sort_indices(t,sort_keys=[(k,'ascending')for k in keys]))
 return t.combine_chunks()
def compare(new,old,name,keys,cols=None):
 import pyarrow as pa
 if cols is None:
  import pyarrow.parquet as pq
  ns=pq.read_schema(new).names;oldcols=pq.read_schema(old).names
  added=set(ns)-set(oldcols);removed=set(oldcols)-set(ns)
  if added or removed:
   record(name+'/schema_annotation',added=={'mapping_rule'} and not removed,str(oldcols),str(ns))
   cols=oldcols
 a=canonical_table(new,cols,keys);b=canonical_table(old,cols,keys)
 record(name,a.equals(b),str(old),f'{new}; rows={a.num_rows}',new_sha256=sha(new),expected_sha256=sha(old),comparison='ALL_SELECTED_VALUES_AFTER_LOSSLESS_STORAGE_NORMALIZATION',columns=a.column_names)
 del a,b
try:
 import importlib.metadata as md,platform
 import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq
 historical=json.loads((S/'environment/HISTORICAL_DEPENDENCY_METADATA.json').read_text())
 for d in historical['dependencies']:record('ENV/'+d['name'],md.version(d['name'])==d['version'],d['version'],md.version(d['name']))
 record('ENV/no_system_site_packages',not any('codex-primary-runtime/dependencies/python/lib/python3.12/site-packages' in p for p in sys.path),'isolated venv',sys.path)
 (S/'environment/RUN_ENVIRONMENT.json').write_text(json.dumps(dict(python=sys.version,executable=sys.executable,platform=platform.platform(),prefix=sys.prefix,base_prefix=sys.base_prefix,sys_path=sys.path,sqlite_version=sqlite3.sqlite_version),indent=2))
 sys.path.insert(0,str(R/'scripts'));import phase1_build_benchmark as frozen
 cfg=json.loads((R/'configs/repair_config.json').read_text());versions=list(cfg['releases'])
 # Tiny fixture exercises the existing functions only. No fabricated study result.
 stage('SMOKE_TEST')
 f=pd.read_parquet(Path(__file__).resolve().parent/'fixtures/activity_smoke.parquet',columns=frozen.READ_COLS)
 p=frozen.prepare_activities(f);e=frozen.primary_edges(frozen.reduce_edges(p),9)
 record('SMOKE/existing_functions',len(e)>0, 'nonempty qualifying fixture',len(e))
 (S/'validation/SMOKE_TEST.json').write_text(json.dumps({'status':'PASS','fixture':'first1000 existing exported24.1 rows; fixture only, not raw reproduction','activities':len(f),'prepared':len(p),'edges':len(e)},indent=2));del f,p,e
 raw=list(csv.DictReader((S/'manifests/RAW_INPUT_SHA256.tsv').open(),delimiter='\t'))
 for version in versions:
  stage('A_RAW_INPUT_INTEGRITY',version=version)
  needed=[r for r in raw if r['version']==version]
  for r in needed:
   path=Path(r['path']);observed=sha(path)
   record('RAW/'+version+'/'+r['kind'],observed==r['SHA256'],r['SHA256'],observed)
  db=next(Path(r['path'])for r in needed if r['kind']=='SQLite')
  con=sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1',uri=True);con.execute('PRAGMA cache_size=-1048576');con.execute('PRAGMA query_only=ON')
  tables={r[0] for r in con.execute("select name from sqlite_master where type='table'")}
  record('RAW/schema/'+version,{'activities','assays','docs','molecule_hierarchy','molecule_dictionary','target_dictionary','version'}<=tables)
  ver=con.execute('SELECT * FROM version').fetchall()
  record('RAW/version/'+version,('ChEMBL_'+version.replace('.','_') in str(ver).replace('.','_')) or (version=='24.1' and ver==[('ChEMBL_24','2018-04-23 00:00:00.000000','ChEMBL Release 24')] and next(r for r in needed if r['kind']=='archive')['SHA256']=='6bb1030408c68b26ad8e9e9ae34ca7226e55d4214a6a057d93be987bbba5ea8c'),version,ver)
  # Fresh full-file SHA256 establishes identical bytes to the prior checked DB.
  # Do not repeat hours of PRAGMA work; report the inherited check as inherited.
  inherited=json.loads((HQ/'logs'/f'integrity_checkpoint_{version}.json').read_text())
  record('RAW/prior_integrity_checkpoint/'+version,inherited['integrity_check']==[['ok']], 'historical ok; same fresh DB SHA256',inherited.get('evidence','recorded checkpoint'))
  (S/'validation'/f'SQLITE_{version}.json').write_text(json.dumps(dict(version=ver,full_integrity_check='NOT_RERUN; prior direct or24.1reconstructed checkpoint retained and fresh byte hash matched',checkpoint=inherited,schema=list(con.execute("select name,sql from sqlite_master where type='table'"))),indent=2));con.close()
  run(Path(__file__).resolve().parent/'cached_export.py','B_EXPORT_'+version,[R/'scripts/phase1_export_chembl_release.py','--db',db,'--version',version,'--out',Q/'derived'/f'activities_{version}.parquet'])
  stage('B_PREPARE_C_BUILD_PRIMARY',version=version)
  prepared=frozen.prepare_activities(pd.read_parquet(Q/'derived'/f'activities_{version}.parquet',columns=frozen.READ_COLS))
  prepared.to_parquet(V/'derived'/f'qualifying_rows_{version}.parquet',index=False)
  primary=frozen.primary_edges(frozen.reduce_edges(prepared),9)
  target=Q/'derived'/('edges_24.1_2018.parquet' if version=='24.1' else f'primary_{version}.parquet');primary.to_parquet(target,index=False)
  compare(target,W/'biomedical_ai_historical_discovery/14_repair/inputs/primary'/f'primary_{version}.parquet','EDGE/primary/'+version,['parent_id','target_id'],['parent_id','target_id'])
  compare(V/'derived'/f'qualifying_rows_{version}.parquet',HQ/'validity/derived'/f'qualifying_rows_{version}.parquet','ACTIVITY/qualifying/'+version,['activity_id'])
  if version=='37':
   for year in [2018,2020,2022]:
    stage('D_BUILD_ROLLBACK',year=year)
    retro=frozen.primary_edges(frozen.reduce_edges(prepared[prepared.document_year.le(year)]),9)
    dest=Q/'derived'/('edges_37_2018.parquet' if year==2018 else f'rollback_37_{year}.parquet');retro.to_parquet(dest,index=False)
    compare(dest,W/f'biomedical_ai_historical_discovery/03_historical_snapshots/CHEMBL_RETRO_{year}/edges.parquet','EDGE/rollback/'+str(year),['parent_id','target_id'],['parent_id','target_id'])
  del prepared,primary
 # All files below are actual PASS gate summaries generated by this run, not copies
 # of historical completion flags.
 stage('E_CORE_COUNTS')
 edgechecks=[dict(label=x['check'],status=x['status'])for x in results if x['check'].startswith('EDGE/')]
 pd.DataFrame(edgechecks).to_csv(Q/'tables/EDGE_MEMBERSHIP_COMPARISON.tsv',sep='\t',index=False);pd.DataFrame(edgechecks).to_csv(Q/'tables/EXPANSION_EDGE_MEMBERSHIP.tsv',sep='\t',index=False)
 core=[]
 for year,spec in cfg['cutoffs'].items():
  a=Q/'derived'/('edges_24.1_2018.parquet'if year=='2018'else f"primary_{spec['release']}.parquet");b=Q/'derived'/('edges_37_2018.parquet'if year=='2018'else f'rollback_37_{year}.parquet')
  si=set(map(tuple,pd.read_parquet(a,columns=['parent_id','target_id']).to_numpy()));ri=set(map(tuple,pd.read_parquet(b,columns=['parent_id','target_id']).to_numpy()))
  obs=dict(strict_edges=len(si),rollback_edges=len(ri),intersection=len(si&ri),strict_only=len(si-ri),rollback_only=len(ri-si),backfill_fraction=len(ri-si)/len(ri),Jaccard=len(si&ri)/len(si|ri))
  ex=pd.read_csv(HQ/'tables/CORE_RESULT_RECALCULATION.tsv',sep='\t');ex=ex[ex.cutoff.eq(int(year))].set_index('metric')
  for n,val in obs.items():
   exp=float(ex.loc[n,'expected']);ok=abs(val-exp)<1e-12;record('CORE/'+year+'/'+n,ok,exp,val);core.append(dict(cutoff=int(year),metric=n,expected=exp,observed=val,difference=val-exp,status='PASS'if ok else'FAIL'))
 pd.DataFrame(core).to_csv(Q/'tables/CORE_RESULT_RECALCULATION.tsv',sep='\t',index=False)
 (Q/'logs/expansion_progress.json').write_text(json.dumps({'status':'RAW_RESCUE_PASS','evidence':'Current clean run all8archive/SQLite fresh hashes and primary/rollback identities passed'}))
 (Q/'logs/pilot_run_result.json').write_text(json.dumps({'status':'PILOT_REPRODUCIBILITY_PASS','evidence':'Current clean run raw exports/graphs passed; this gate does not claim labels passed'}))
 run(Q/'scripts/rebuild_label_chain.py','F_TO_J_LABEL_CHAIN')
 compare(Q/'derived/REBUILT_RELEASE_MEMBERSHIP.parquet',HQ/'derived/REBUILT_RELEASE_MEMBERSHIP.parquet','LABEL/release_membership',['parent_id','target_id'])
 for year in [2018,2020,2022]:
  for stem in ['REBUILT_U2_LABELS','REBUILT_HISTORICAL_ACTIVE']:compare(Q/'derived'/f'{stem}_{year}.parquet',HQ/'derived'/f'{stem}_{year}.parquet','LABEL/'+stem+'/'+str(year),['parent_id','target_id'])
 run(V/'scripts/extract_entities.py','K_RAW_ENTITIES')
 run(V/'scripts/entity_mapping.py','K_EXACT_ENTITY_MAPPING')
 run(V/'scripts/canonical_structure_evidence.py','K_CANONICAL_EVIDENCE')
 env=dict(os.environ,ENTITY_AUDIT_ROOT=str(V/'structure_equivalence'),ENTITY_MAPPING_MODE='CANONICAL_STRUCTURE_EQUIVALENCE')
 run(V/'scripts/entity_mapping.py','K_STRUCTURE_ENTITY_MAPPING',env=env)
 run(V/'scripts/cause_attribution.py','K_BASELINE_CAUSE_ATTRIBUTION')
 run(V/'scripts/cause_attribution.py','K_CAUSE_ATTRIBUTION',env=env)
 run(V/'scripts/qualification_failure_evidence.py','K_QUALIFICATION_FAILURES')
 run(V/'scripts/future_label_cause_evidence.py','K_FUTURE_LABEL_CAUSES')
 for root,expected in [(V,HQ/'validity'),(V/'structure_equivalence',HQ/'validity/structure_equivalence')]:
  for p in sorted((root/'derived').glob('*.parquet')):
   if p.name.startswith('._'):continue
   if p.name.startswith(('qualifying_rows_','molecules_','target_components_','documents_')):continue
   old=expected/'derived'/p.name
   if not old.exists():continue
   names=pq.read_schema(p).names;keys=['parent_id','target_id','activity_id']if all(n in names for n in ['parent_id','target_id','activity_id'])else['chembl_id']if 'chembl_id'in names else['local_id']if 'local_id'in names else['publisher_smiles']if 'publisher_smiles'in names else['parent_id','target_id']if all(n in names for n in ['parent_id','target_id'])else['activity_id']if'activity_id'in names else names[:1]
   compare(p,old,'RECORD/'+str(p.relative_to(V)),keys)
 # Entire cause classes, membership, rule/evidence flags must agree, not counts alone.
 for name in ['EDGE_CAUSE_ATTRIBUTION.tsv','EDGE_CAUSE_SUMMARY.tsv','STABLE_MAPPING_SENSITIVITY.tsv']:
  a=pd.read_csv(V/'structure_equivalence/tables'/name,sep='\t');b=pd.read_csv(HQ/'validity/structure_equivalence/tables'/name,sep='\t')
  keys=[c for c in ['cutoff','direction','parent_id','target_id','category']if c in a.columns];a=a.sort_values(keys).reset_index(drop=True);b=b.sort_values(keys).reset_index(drop=True);record('RECORD/'+name,a.equals(b),len(b),len(a),new_sha256=sha(V/'structure_equivalence/tables'/name),expected_sha256=sha(HQ/'validity/structure_equivalence/tables'/name))
 for name in ['MOLECULE_TYPE_MISSINGNESS_EVIDENCE.tsv','FUTURE_LABEL_MISSINGNESS_EVIDENCE.tsv','FUTURE_LABEL_CAUSE_SUMMARY.tsv']:
  a=pd.read_csv(V/'tables'/name,sep='\t');b=pd.read_csv(HQ/'validity/tables'/name,sep='\t');keys=[c for c in ['cutoff','domain','direction','category']if c in a.columns];record('RECORD/'+name,a.sort_values(keys).reset_index(drop=True).equals(b.sort_values(keys).reset_index(drop=True)),len(b),len(a))
 state['status']='FULL_REPRODUCTION_PASS';stage('COMPLETE')
except Exception as e:
 state.update(status='FAIL',error_type='RESULT_DISCREPANCY'if'RESULT_DISCREPANCY'in str(e)else'ENVIRONMENT_OR_CODE_OR_INPUT_ERROR',error=repr(e));stage('STOPPED');traceback.print_exc();raise
