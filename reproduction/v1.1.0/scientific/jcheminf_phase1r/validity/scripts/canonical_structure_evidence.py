"""Compare publisher SMILES under one RDKit version without tautomer or parent transformations."""
import sys
sys.dont_write_bytecode = True
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import datetime
import json
import os
import pandas as pd
from rdkit import Chem, RDLogger, rdBase

V = Path(__file__).resolve().parents[1]

def canonicalize(smiles):
    RDLogger.DisableLog('rdApp.*')
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return smiles, None
    return smiles, Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)

def main():
    versions = list(json.loads((V.parent.parent / 'biomedical_ai_historical_discovery/configs/repair_config.json').read_text())['releases'])
    strings = set()
    for version in versions:
        f = pd.read_parquet(V / 'derived' / f'compound_mapping_{version}.parquet', columns=['canonical_smiles', 'reference_smiles', 'mapping_status'])
        f = f[f.mapping_status.eq('SMILES_CHANGED_WITH_SAME_INCHIKEY')]
        strings.update(f.canonical_smiles.dropna())
        strings.update(f.reference_smiles.dropna())
    records = []
    progress = V / 'logs/structure_equivalence_progress.json'
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i, (source, canonical) in enumerate(pool.map(canonicalize, sorted(strings), chunksize=500), 1):
            records.append({'publisher_smiles': source, 'rdkit_canonical_isomeric_smiles': canonical})
            if i % 25000 == 0:
                status = {'status': 'RUNNING', 'pid': os.getpid(), 'processed': i, 'total': len(strings), 'updated_at': datetime.datetime.now().astimezone().isoformat()}
                progress.write_text(json.dumps(status, indent=2))
                print(status, flush=True)
    result = pd.DataFrame(records)
    result.to_parquet(V / 'derived/SMILES_REPRESENTATION_EQUIVALENCE.parquet', index=False)
    status = {'status': 'COMPLETE', 'processed': len(result), 'parse_failures': int(result.rdkit_canonical_isomeric_smiles.isna().sum()), 'rdkit_version': rdBase.rdkitVersion, 'rule': 'Same complete publisher InChIKey plus matching RDKit canonical isomeric SMILES; no tautomer, charge, fragment or parent transformation'}
    progress.write_text(json.dumps(status, indent=2))
    print(status, flush=True)

if __name__ == '__main__':
    main()
