"""Raw release-local entity evidence for cross-version mapping."""
import sys, os
sys.dont_write_bytecode=True
from pathlib import Path
import sqlite3, json, datetime, time
import pandas as pd
import pyarrow.parquet as pq
V=Path(__file__).resolve().parents[1]; Q=V.parent
versions=list(json.loads((Q.parent/'biomedical_ai_historical_discovery/configs/repair_config.json').read_text())['releases'])
def progress(stage,version):
 d=dict(status='RUNNING',pid=os.getpid(),stage=stage,version=version,updated_at=datetime.datetime.now().astimezone().isoformat())
 (V/'logs/entity_extraction_progress.json').write_text(json.dumps(d,indent=2));print(d,flush=True)
for v in versions:
 done=V/'logs'/f'entities_{v}_complete.json'
 if done.exists():continue
 progress('RAW_ENTITY_EXTRACTION',v)
 db=next(p for p in (Q/'raw'/('extract_'+v)).glob('*.db') if not p.name.startswith('._'))
 con=sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1',uri=True)
 con.execute('PRAGMA cache_size=-1048576')
 # Full raw tables allow ambiguity checks beyond the qualifying graph subset.
 query='''SELECT md.molregno,md.chembl_id, h.parent_molregno, pm.chembl_id AS hierarchy_parent_id,
 cs.canonical_smiles,cs.standard_inchi_key
 FROM molecule_dictionary md
 LEFT JOIN molecule_hierarchy h ON h.molregno=md.molregno
 LEFT JOIN molecule_dictionary pm ON pm.molregno=h.parent_molregno
 LEFT JOIN compound_structures cs ON cs.molregno=md.molregno'''
 writer=None;n=0
 try:
  for frame in pd.read_sql_query(query,con,chunksize=100000):
   import pyarrow as pa
   # Explicit stable types for batches with missing structure values.
   for c in ['chembl_id','hierarchy_parent_id','canonical_smiles','standard_inchi_key']:frame[c]=frame[c].astype('string')
   for c in ['molregno','parent_molregno']:frame[c]=pd.to_numeric(frame[c]).astype('Int64')
   table=pa.Table.from_pandas(frame,preserve_index=False)
   if writer is None:writer=pq.ParquetWriter(V/'derived'/f'molecules_{v}.parquet',table.schema,compression='zstd')
   writer.write_table(table);n+=len(frame)
   print(v,'molecules',n,flush=True)
 finally:
  if writer:writer.close()
 progress('RAW_TARGET_EXTRACTION',v)
 target=pd.read_sql_query('''SELECT td.tid,td.chembl_id,td.target_type,td.tax_id,td.organism,
 tc.component_id,tc.homologue,cs.accession,cs.sequence_md5sum,cs.component_type,cs.tax_id AS component_tax_id
 FROM target_dictionary td LEFT JOIN target_components tc ON tc.tid=td.tid
 LEFT JOIN component_sequences cs ON cs.component_id=tc.component_id''',con)
 target.to_parquet(V/'derived'/f'target_components_{v}.parquet',index=False)
 progress('RAW_DOCUMENT_EXTRACTION',v)
 docs=pd.read_sql_query('SELECT doc_id,chembl_id,year,doi,pubmed_id,doc_type FROM docs',con)
 docs.to_parquet(V/'derived'/f'documents_{v}.parquet',index=False)
 con.close()
 done.write_text(json.dumps(dict(version=v,molecules=n,target_component_rows=len(target),documents=len(docs),db_size=db.stat().st_size,db_mtime_ns=db.stat().st_mtime_ns,raw_validation_file=str(Q/'logs'/f'sqlite_validation_{v}.json')),indent=2))
progress('COMPLETE','')
p=V/'logs/entity_extraction_progress.json';d=json.loads(p.read_text());d['status']='COMPLETE';p.write_text(json.dumps(d,indent=2))
