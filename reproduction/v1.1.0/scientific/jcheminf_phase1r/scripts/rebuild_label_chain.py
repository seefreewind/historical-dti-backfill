"""Independent inputs for candidates and first-observed labels, frozen rules reused."""
import os, sys
sys.dont_write_bytecode = True
from pathlib import Path
import ast, json, datetime, hashlib, traceback, fcntl
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
Q = Path(__file__).resolve().parents[1]
R = Q.parent / 'biomedical_ai_historical_discovery'
os.environ['MPLCONFIGDIR'] = str(Q/'logs/matplotlib_cache')
sys.path.insert(0, str(R/'scripts'))
import phase1_build_benchmark as frozen
CFG = json.loads((R/'configs/repair_config.json').read_text())
lock = (Q/'logs/label_chain.lock').open('a')
fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
state = dict(status='RUNNING', pid=os.getpid(), started_at=datetime.datetime.now().astimezone().isoformat())
checks = []
def stage(name, year=''):
    state.update(stage=name, cutoff=year, updated_at=datetime.datetime.now().astimezone().isoformat())
    (Q/'logs/label_chain_progress.json').write_text(json.dumps(state, indent=2))
    print(state['updated_at'], name, year, flush=True)
def check(name, expected, observed, year=''):
    ok = expected == observed
    checks.append(dict(cutoff=year, check=name, expected=str(expected), observed=str(observed), status='PASS' if ok else 'FAIL'))
    pd.DataFrame(checks).to_csv(Q/'tables/LABEL_RECONSTRUCTION_AUDIT.tsv', sep='\t', index=False)
    if not ok: raise RuntimeError('HOLD_RESULT_DISCREPANCY '+str(year)+' '+name)
def encode(frame):
    c = compounds.get_indexer(frame.parent_id)
    t = targets.get_indexer(frame.target_id)
    assert (c>=0).all() and (t>=0).all()
    return c.astype(np.int64)*nt+t
try:
    stage('VERIFY_RAW_GRAPH_GATE')
    assert json.loads((Q/'logs/expansion_progress.json').read_text())['status']=='RAW_RESCUE_PASS'
    assert json.loads((Q/'logs/pilot_run_result.json').read_text())['status']=='PILOT_REPRODUCIBILITY_PASS'
    for file in ['CORE_RESULT_RECALCULATION.tsv','EDGE_MEMBERSHIP_COMPARISON.tsv','EXPANSION_EDGE_MEMBERSHIP.tsv']:
        assert pd.read_csv(Q/'tables'/file,sep='\t').status.eq('PASS').all()
    versions = list(CFG['releases'])
    dates = np.array([CFG['releases'][v] for v in versions], dtype='datetime64[D]')
    # Pilot materialized only the 2018 rollback for 37. Reconstruct its full primary graph now.
    stage('RECONSTRUCT_FULL_PRIMARY_37')
    prepared = frozen.prepare_activities(pd.read_parquet(Q/'derived/activities_37.parquet', columns=frozen.READ_COLS))
    full37 = frozen.primary_edges(frozen.reduce_edges(prepared), 9)
    full37.to_parquet(Q/'derived/primary_37.parquet', index=False)
    del prepared, full37
    paths = {v:Q/'derived'/('edges_24.1_2018.parquet' if v=='24.1' else f'primary_{v}.parquet') for v in versions}
    frames = {v:pd.read_parquet(p, columns=['parent_id','target_id']).drop_duplicates() for v,p in paths.items()}
    stage('BUILD_RELEASE_MEMBERSHIP')
    compounds = pd.Index(pd.unique(np.concatenate([x.parent_id.to_numpy() for x in frames.values()])))
    targets = pd.Index(pd.unique(np.concatenate([x.target_id.to_numpy() for x in frames.values()])))
    nt = len(targets)
    parts = [pd.DataFrame(dict(key=encode(f), membership=np.uint16(1<<j))) for j,f in enumerate(frames.values())]
    union = pd.concat(parts,ignore_index=True).groupby('key',sort=True).membership.sum().astype(np.uint16)
    del parts
    keys, bits = union.index.to_numpy(), union.to_numpy()
    first = np.full(len(keys),-1,np.int16)
    for j in range(len(versions)): first[(first<0)&((bits&(1<<j))!=0)] = j
    membership_map = {x:'|'.join(v for j,v in enumerate(versions) if int(x)&(1<<j)) for x in np.unique(bits)}
    master = pd.DataFrame(dict(parent_id=compounds.to_numpy()[keys//nt],target_id=targets.to_numpy()[keys%nt],release_membership_bitmask=bits,all_release_membership=[membership_map[x] for x in bits],first_public_release_present=[versions[j] for j in first],first_operational_date=dates[first].astype(str)))
    master.to_parquet(Q/'derived/REBUILT_RELEASE_MEMBERSHIP.parquet',index=False)
    # Extract only the frozen label function. Its inputs all come from the new graphs.
    tree = ast.parse((R/'scripts/build_frozen_first_seen_labels.py').read_text())
    build = next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='build')
    label_node = next(x for x in build.body if isinstance(x,ast.FunctionDef) and x.name=='label')
    exec(compile(ast.Module(body=[label_node],type_ignores=[]), '<frozen_label_function>', 'exec'),globals())
    expected_master = pd.read_parquet(R/'14_repair/AUTHORITATIVE_RELEASE_MEMBERSHIP.parquet')
    oldkeys = encode(expected_master)
    order = np.argsort(oldkeys)
    check('release_membership_pair_identity',True,np.array_equal(keys,oldkeys[order]))
    check('release_membership_bitmasks',True,np.array_equal(bits,expected_master.release_membership_bitmask.to_numpy()[order]))
    del expected_master, oldkeys, order
    expected_summary = pd.read_csv(R/'14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv',sep='\t').set_index('cutoff')
    summaries = []
    compare_columns = ['first_discovery','gap_reappearance','present_before_cutoff','present_at_cutoff','right_censored','later_first_appearance_after_followup']
    for year in [2018,2020,2022]:
        stage('REBUILD_CANDIDATE_UNIVERSE',year)
        spec = CFG['cutoffs'][str(year)]
        strict = pd.read_parquet(paths[spec['release']])
        known, target_classes, usable = frozen.u2_structures(strict)
        # Frozen candidate/label files are opened only for expected identity and metadata comparison.
        expected_frame = pd.read_parquet(R/f'14_repair/AUTHORITATIVE_U2_LABELS_{year}.parquet',columns=['parent_id','target_id']+compare_columns)
        expected_keys = encode(expected_frame)
        order = np.argsort(expected_keys)
        sorted_expected = expected_keys[order]
        expected_values = expected_frame[compare_columns].to_numpy(dtype=np.int8)[order]
        del expected_frame, expected_keys, order
        check('expected_candidate_uniqueness',len(sorted_expected),len(np.unique(sorted_expected)),year)
        retro_path = Q/'derived'/('edges_37_2018.parquet' if year==2018 else f'rollback_37_{year}.parquet')
        retro_keys = np.sort(encode(pd.read_parquet(retro_path,columns=['parent_id','target_id'])))
        writer = None
        total = positives = erased = gap = prior = mismatch = unknown = duplicates = 0
        parent_names = sorted(known)
        try:
            for start in range(0,len(parent_names),1000):
                pairs = [(p,t) for p in parent_names[start:start+1000] for t in sorted(frozen.u2_pool(p,known,target_classes,usable))]
                if not pairs: continue
                candidate = pd.DataFrame(pairs,columns=['parent_id','target_id'])
                new = label(candidate,year)
                k = encode(candidate)
                duplicates += len(k)-len(np.unique(k))
                loc = np.searchsorted(sorted_expected,k)
                safe = np.minimum(loc,len(sorted_expected)-1)
                valid = (loc<len(sorted_expected))&(sorted_expected[safe]==k)
                unknown += int((~valid).sum())
                mismatch += int(np.any(new[compare_columns].to_numpy(dtype=np.int8)[valid]!=expected_values[safe[valid]],axis=1).sum())
                event = new.first_discovery.to_numpy(dtype=bool)
                total += len(new); positives += int(event.sum())
                erased += int((event & np.isin(k,retro_keys)).sum())
                gap += int(new.gap_reappearance.sum()); prior += int(new.present_before_cutoff.sum())
                assert not (event & (new.present_before_cutoff | new.gap_reappearance).to_numpy()).any()
                assert not new.present_at_cutoff.any()
                if event.any():
                    d=np.array(new.loc[event,'event_release_date'],dtype='datetime64[D]')
                    assert ((d>np.datetime64(spec['cutoff']))&(d<=np.datetime64(spec['followup_end']))).all()
                table = pa.Table.from_pandas(new,preserve_index=False)
                if writer is None: writer=pq.ParquetWriter(Q/'derived'/f'REBUILT_U2_LABELS_{year}.parquet',table.schema)
                writer.write_table(table)
        finally:
            if writer: writer.close()
        check('candidate_count',len(sorted_expected),total,year)
        check('candidate_identity_unknown',0,unknown,year)
        check('candidate_duplicates',0,duplicates,year)
        check('all_six_label_flags_pairwise_differences',0,mismatch,year)
        check('frozen_candidate_count',{2018:9294785,2020:10318524,2022:12123841}[year],total,year)
        del sorted_expected, expected_values
        stage('REBUILD_HISTORICAL_ACTIVE_FIRST_APPEARANCES',year)
        possible = master[master.parent_id.isin(strict.parent_id.unique()) & master.target_id.isin(strict.target_id.unique())]
        broad = label(possible,year)
        broad = broad[~broad.present_at_cutoff].copy()
        broad.to_parquet(Q/'derived'/f'REBUILT_HISTORICAL_ACTIVE_{year}.parquet',index=False)
        bpos = broad[broad.first_discovery.eq(1)]
        broad_erased = int(np.isin(encode(bpos),retro_keys).sum())
        expected_broad = pd.read_parquet(R/f'14_repair/LABELS_HISTORICAL_ACTIVE_OBSERVED_{year}.parquet',columns=['parent_id','target_id']+compare_columns)
        left=broad.set_index(['parent_id','target_id'])[compare_columns].sort_index()
        right=expected_broad.set_index(['parent_id','target_id'])[compare_columns].sort_index()
        check('historical_active_pair_and_label_identity',True,left.equals(right),year)
        observed = dict(u2_positive=positives,u2_erased=erased,historical_active_positive=len(bpos),historical_active_erased=broad_erased)
        for name,value in observed.items(): check(name,int(expected_summary.loc[year,name]),value,year)
        summaries.append(dict(cutoff=year,candidates=total,**observed,gap_reappearances=gap,prior_appearances=prior,followup_end=spec['followup_end'],erasure_fraction=erased/positives))
        pd.DataFrame(summaries).to_csv(Q/'tables/REBUILT_LABEL_SUMMARY.tsv',sep='\t',index=False)
        del broad, expected_broad, left, right, strict
    state['status']='CORE_LABEL_CHAIN_REPRODUCIBILITY_PASS'
    state['verdict_note']='Full record-level core candidate and label comparison passed. Final manifests and reports still require completion.'
    stage('CORE_CHAIN_COMPLETE')
except Exception as exc:
    state['status']='HOLD_RESULT_DISCREPANCY' if 'HOLD_RESULT_DISCREPANCY' in str(exc) else 'HOLD_REPRODUCIBILITY'
    state['error']=repr(exc)
    stage('STOPPED')
    traceback.print_exc()
    raise
