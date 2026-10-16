from pathlib import Path
import hashlib,csv
base=Path(__file__).resolve().parents[1]
with (base/'MANIFESTS/PACKAGE_SHA256.tsv').open() as f:
 rows=list(csv.DictReader(f,delimiter='\t'))
for r in rows:
 p=base/r['path'];actual=hashlib.sha256(p.read_bytes()).hexdigest()
 assert actual==r['sha256'],r['path']
print('PASS:',len(rows),'packaged file hashes')
