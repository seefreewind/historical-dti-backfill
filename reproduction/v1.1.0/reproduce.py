"""Portable launcher for existing fixed implementations; raw inputs remain external."""
from pathlib import Path
import argparse,json,csv,shutil,os,subprocess,sys,hashlib
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('command',choices=['configure','check','smoke','run'])
ap.add_argument('--run-root',type=Path,required=True)
ap.add_argument('--raw-root',type=Path)
a=ap.parse_args();S=a.run_root.resolve()
def digest(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
if a.command=='configure':
 if not a.raw_root:ap.error('--raw-root required')
 if S.exists()and any(S.iterdir()):ap.error('run-root must be empty; previous runs cannot be overwritten')
 for n in ['logs','validation','manifests','environment','scripts']: (S/n).mkdir(parents=True,exist_ok=True)
 shutil.copyfile(P/'environment/HISTORICAL_DEPENDENCY_METADATA.json',S/'environment/HISTORICAL_DEPENDENCY_METADATA.json')
 O=S/'reproduction/work';shutil.copytree(P/'scientific',O)
 Q=O/'jcheminf_phase1r';R=O/'biomedical_ai_historical_discovery'
 for root in [Q,Q/'validity',Q/'validity/structure_equivalence']:
  for n in ['derived','logs','tables','manifests']:(root/n).mkdir(parents=True,exist_ok=True)
 for folder in ['14_repair','03_historical_snapshots']:
  shutil.copytree(P/'expected/biomedical_ai_historical_discovery'/folder,R/folder)
 raw=list(csv.DictReader((P/'manifests/RAW_INPUTS.tsv').open(),delimiter='\t'))
 for r in raw:
  r['path']=str((a.raw_root.resolve()/r['relative_path']).resolve())
  if not Path(r['path']).is_file():raise FileNotFoundError(r['relative_path'])
  if r['kind']=='SQLite':
   dst=Q/'raw'/('extract_'+r['version']);dst.mkdir(parents=True,exist_ok=True);(dst/Path(r['path']).name).symlink_to(r['path'])
 with(S/'manifests/RAW_INPUT_SHA256.tsv').open('x',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(raw[0]),delimiter='\t');w.writeheader();w.writerows(raw)
 (S/'configuration.json').write_text(json.dumps({'package':str(P),'raw_root':str(a.raw_root.resolve()),'run_root':str(S),'network_required_during_run':False},indent=2))
 print('CONFIGURATION_INITIALIZED')
elif a.command=='check':
 raw=list(csv.DictReader((S/'manifests/RAW_INPUT_SHA256.tsv').open(),delimiter='\t'))
 for r in raw:
  assert Path(r['path']).is_file(),r['relative_path']
 print('INPUT_LOCATIONS_PASS',len(raw),'hashes verified by run command')
elif a.command=='smoke':
 import pandas as pd
 O=S/'reproduction/work';os.environ['HISTORICAL_DTI_ROOT']=str(O/'biomedical_ai_historical_discovery');os.environ['HISTORICAL_DTI_CONFIG']=str(O/'biomedical_ai_historical_discovery/configs/repair_config.json')
 sys.path.insert(0,str(O/'biomedical_ai_historical_discovery/scripts'));import phase1_build_benchmark as f
 x=pd.read_parquet(P/'fixtures/activity_smoke.parquet',columns=f.READ_COLS);y=f.primary_edges(f.reduce_edges(f.prepare_activities(x)),9);assert len(y)>0
 (S/'validation/ENTRY_SMOKE.json').write_text(json.dumps({'status':'PASS','fixture_rows':len(x),'edges':len(y)}));print('SMOKE_PASS',len(y))
else:
 env=dict(os.environ,REPRODUCTION_RUN_ROOT=str(S),PYTHONDONTWRITEBYTECODE='1')
 subprocess.run([sys.executable,'-I','-B',str(P/'run_chain.py')],env=env,check=True)
