"""Export traceable plotting values from frozen summary tables, no model execution.

Run inside reviewer package CODE after copying this file there. JSON and TSV
outputs use Python only; optional XLSX authoring is a separate JavaScript step.
"""
from pathlib import Path
import json,hashlib
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
DATA=BASE/'DATA_DERIVED'; OUT=BASE/'REPRODUCED'; OUT.mkdir(exist_ok=True)
def rd(p):return pd.read_csv(DATA/p,sep='\t')
def h(p):return hashlib.sha256((DATA/p).read_bytes()).hexdigest()
back=rd('10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv');arms=rd('11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv');labs=rd('14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv')
top=rd('14_repair/TOP_PERCENTILE_TIE_AUDIT.tsv');tie=rd('14_repair/COMPOUND_CENTRIC_TIE_REPAIR.tsv');chem=rd('14_repair/CHEMICAL_CLASSIFICATION_SUMMARY.tsv');bins=rd('14_repair/FINGERPRINT_BIN_COUNTS.tsv');effect=rd('14_repair/NESTED_GRAPH_EFFECTS.tsv');nested=rd('14_repair/NESTED_GRAPH_ABLATION.tsv')
payload=[];manifest=[]
def add(name,rows,sources,columns,figure=None,panel=None):
 df=pd.DataFrame(rows).fillna('NOT_APPLICABLE')
 payload.append({'name':name,'matrix':[df.columns.tolist()]+df.values.tolist()})
 for source in sources:manifest.append(dict(figure=figure or 'Figure'+name[3],panel=panel or name[4:],source_data_sheet=name,original_artifact=source,original_columns=columns,transformation_script='CODE/ncomms_source_data.py',checksum=h(source)))
def row(year,group,metric,value,denom,n,src,**kw):
 return dict(year=year,group_model=group,metric=metric,value=value,denominator=denom,n=n,uncertainty_lower='NOT_PLOTTED',uncertainty_upper='NOT_PLOTTED',source_artifact=src,**kw)
for p in 'abcd':add('Fig1'+p,[dict(status='SCHEMATIC_NO_NUMERICAL_SOURCE_DATA',source_manifest='MANIFESTS/FIGURE_CONTENT_MANIFEST.tsv',cutoffs='2018;2020;2022 censored',uncertainty='NOT_APPLICABLE')],['14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv','14_repair/CHEMBL_RELEASE_DATE_PROVENANCE.tsv'],'design and date provenance')
src='10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv'
add('Fig2a',[row(int(r.cutoff),c,'edge_count',int(r[c]),int(r.N_overlap+r.N_strict_only+r.N_retro_only),int(r[c]),src,denominator_definition='Union of strict and rollback edges') for _,r in back.iterrows() for c in ['N_overlap','N_strict_only','N_retro_only']],[src],'cutoff;N_overlap;N_strict_only;N_retro_only')
add('Fig2b',[row(int(r.cutoff),'rollback','rollback_only_percent',100*r.backfill_fraction,int(r.N_edges_retro),int(r.N_retro_only),src) for _,r in back.iterrows()],[src],'cutoff;backfill_fraction;N_edges_retro;N_retro_only')
src='14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv'
ers=[row(int(r.cutoff),g,'future_positive_erasure_percent',100*r[n]/r[d],int(r[d]),int(r[n]),src,denominator_definition=g+' future-positive pairs') for _,r in labs.iterrows() for g,n,d in [('U2','u2_erased','u2_positive'),('historical_active','historical_active_erased','historical_active_positive')]]
add('Fig2c',ers,[src,'11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv'],'u2_positive;u2_erased;historical_active_positive;historical_active_erased;label_erasure_fraction_of_u2_positives')
src='11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv'
add('Fig2d',[row(int(r.cutoff),'U2','candidates_millions',r.u2_n/1e6,int(r.u2_n),int(r.u2_n),src) for _,r in arms.iterrows()],[src],'cutoff;u2_n')
lig='M5_ligand_max';deg='M2_compound_degree';tar='M1_target_degree'
def globalrows(d):
 return [row(int(r.cutoff),r.model,'full_boundary_tie_fde',r.full_boundary_tie_fde,int(r.n_candidates),int(r.actual_evaluated_n),'14_repair/TOP_PERCENTILE_TIE_AUDIT.tsv',subset=r.subset,n_positives=int(r.n_positives),selected_hits=r.full_boundary_tie_hits,base_rate=r.base_rate,nominal_k=int(r.nominal_k),denominator_definition='All subset candidates; n is actual selected boundary-tie set',post_hoc=True) for _,r in d.iterrows()]
def compoundrows(d,metric):
 return [row(int(r.cutoff),r.model,metric,r[metric],int(r.n_positive_compounds),int(r.n_positive_compounds),'14_repair/COMPOUND_CENTRIC_TIE_REPAIR.tsv',subset=r.subset,denominator_definition='Positive compounds within named subset') for _,r in d.iterrows()]
gsrc='14_repair/TOP_PERCENTILE_TIE_AUDIT.tsv';tsrc='14_repair/COMPOUND_CENTRIC_TIE_REPAIR.tsv'
add('Fig3a',globalrows(top[top.subset.eq('ALL_U2')&top.model.isin([lig,deg])]),[gsrc],'cutoff;model;subset;full_boundary_tie_fde;n_candidates;n_positives;actual_evaluated_n;full_boundary_tie_hits;base_rate;nominal_k')
add('Fig3b',compoundrows(tie[tie.subset.eq('ALL_U2')&tie.model.isin([lig,deg,tar])],'expected_recall_at_1'),[tsrc],'cutoff;model;subset;expected_recall_at_1;n_positive_compounds')
add('Fig3c',compoundrows(tie[tie.subset.eq('ALL_U2')&tie.model.isin([lig,deg])],'top_tie_any_hit_fraction'),[tsrc],'cutoff;model;subset;top_tie_any_hit_fraction;n_positive_compounds')
add('Fig3d',globalrows(top[top.cutoff.eq(2022)&top.subset.eq('ALL_U2')&top.model.isin([lig,deg,'M7_svd','M9_lightgcn'])]),[gsrc],'2022;ALL_U2;full_boundary_tie_fde;n_candidates;n_positives;actual_evaluated_n')
src='14_repair/CHEMICAL_CLASSIFICATION_SUMMARY.tsv'
add('Fig4a',[row(int(r.cutoff),r.chemical_class,'fingerprint_identity_class_percent',r.fraction*100,int(r.denominator),int(r.n_pairs),src,n_future_positives=int(r.n_future_positives),denominator_definition='Candidate–historical-ligand pairs with Tanimoto=1') for _,r in chem.iterrows()],[src],'cutoff;chemical_class;fraction;denominator;n_pairs;n_future_positives')
src='14_repair/FINGERPRINT_BIN_COUNTS.tsv'
add('Fig4b',[row(int(r.cutoff),r.stratum,'future_positive_count',int(r.n_positive),int(r.n_candidates),int(r.n_positive),src) for _,r in bins[bins.cutoff.eq(2018)].iterrows()],[src],'cutoff;stratum;n_positive;n_candidates')
add('Fig4c',globalrows(top[top.cutoff.isin([2018,2020])&top.model.eq(lig)&top.subset.isin(['ALL_U2','TANIMOTO_LT_0.95','TANIMOTO_LT_0.80','NEW'])]),[gsrc],'cutoff;model;subset;full_boundary_tie_fde;n_candidates;n_positives;actual_evaluated_n')
add('Fig4d',globalrows(top[top.subset.eq('NEW')&top.model.isin([lig,deg])]),[gsrc],'cutoff;model;subset;full_boundary_tie_fde;n_candidates;n_positives;actual_evaluated_n')
add('Fig5a',[r for r in ers if r['group_model']=='U2'],['14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv', '11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv'],'u2_positive;u2_erased;label_erasure_fraction_of_u2_positives')
for p,ename,arm in [('b','ADD','G1_ADD_ONLY'),('c','REMOVE','G2_INTERSECTION'),('d','FULL','G3_ROLLBACK')]:
 rows=[]
 for _,r in effect[effect.effect.eq(ename)].iterrows():
  d=nested[nested.cutoff.eq(r.cutoff)&nested.model.eq(r.model)]
  # Match frozen arm identifiers semantically, never reconstruct scores.
  a=d[d.arm.str.startswith('G'+str({'ADD':1,'REMOVE':2,'FULL':3}[ename]))].iloc[0];b=d[d.arm.str.startswith('G0')].iloc[0]
  rows.append(row(int(r.cutoff),r.model,'delta_full_boundary_tie_fde',r.delta_full_boundary_tie_fde,int(a.n_candidates),int(a.actual_evaluated_n),'14_repair/NESTED_GRAPH_EFFECTS.tsv;14_repair/NESTED_GRAPH_ABLATION.tsv',effect=ename,arm=a.arm,baseline_arm=b.arm,n_positives=int(a.n_positives),arm_fde=a.full_boundary_tie_fde,baseline_fde=b.full_boundary_tie_fde,baseline_selected_n=int(b.actual_evaluated_n),base_rate=a.base_rate,denominator_definition='Common-clean candidate universe; n is arm selected count',post_hoc=True))
 add('Fig5'+p,rows,['14_repair/NESTED_GRAPH_EFFECTS.tsv','14_repair/NESTED_GRAPH_ABLATION.tsv'],'cutoff;model;effect;delta_full_boundary_tie_fde;arm;n_candidates;n_positives;actual_evaluated_n;full_boundary_tie_fde')
dates=rd('14_repair/CHEMBL_RELEASE_DATE_PROVENANCE.tsv');table=[]
summary=json.loads((DATA/'10_phase1/phase1_summary.json').read_text())['summary_counts']
for _,r in back.iterrows():
 y=int(r.cutoff);a=arms[arms.cutoff.eq(y)].iloc[0];l=labs[labs.cutoff.eq(y)].iloc[0];v=str(r.strict_release).removesuffix('.0');date=dates[dates.chembl_version.astype(str).str.removesuffix('.0').eq(v)].iloc[0]
 vals={'cutoff':str(y)+'-12-31','historically_active_compounds':summary[str(y)]['compounds'],'historically_active_targets':summary[str(y)]['targets'],'operational_availability_proxy':date.operational_date_used,'strict_release':v,'strict_edges':int(r.N_edges_strict),'rollback_edges':int(r.N_edges_retro),'overlap_edges':int(r.N_overlap),'rollback_only':int(r.N_retro_only),'rollback_only_percent':100*r.backfill_fraction,'strict_only':int(r.N_strict_only),'U2_candidates':int(a.u2_n),'U2_positives':int(l.u2_positive),'U2_erasure_count':int(l.u2_erased),'U2_erasure_percent':100*l.u2_erased/l.u2_positive,'broader_positives':int(l.historical_active_positive),'broader_erasure_count':int(l.historical_active_erased),'broader_erasure_percent':100*l.historical_active_erased/l.historical_active_positive,'followup_end':l.followup_end}
 for metric,value in vals.items():
  denom=int(l.u2_positive) if metric.startswith('U2_erasure') else int(l.historical_active_positive) if metric.startswith('broader_erasure') else int(r.N_edges_retro) if metric=='rollback_only_percent' else 'NOT_APPLICABLE'
  table.append(row(y,'Table1',metric,value,denom,value if isinstance(value,int) else 'NOT_APPLICABLE',';'.join(['10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv','11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv','14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv','14_repair/CHEMBL_RELEASE_DATE_PROVENANCE.tsv'])))
add('Table1',table,['10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv','11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv','14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv','14_repair/CHEMBL_RELEASE_DATE_PROVENANCE.tsv','10_phase1/phase1_summary.json'],'Table 1 fields; summary_counts.compounds/targets only (legacy outcomes not used)',figure='Table1',panel='whole')
# Supplementary values are explicit too: the unchanged V3 workbook is supporting,
# while these tabs close the bin-count and schematic-provenance coverage gaps.
for i,source in [(1,'14_repair/CHEMBL_RELEASE_DATE_PROVENANCE.tsv'),(3,'14_repair/ANALYSIS_DECISION_TIMELINE.tsv'),(4,'11_phase2/U2_SENSITIVITY.tsv'),(5,'14_repair/AUTHORITATIVE_LABEL_SUMMARY.tsv'),(7,'14_repair/NESTED_GRAPH_ARM_MANIFEST.tsv')]:
 add('FigS'+str(i),[dict(status='SCHEMATIC_NO_NUMERICAL_SOURCE_DATA',source_artifact=source,source_manifest='MANIFESTS/FIGURE_CONTENT_MANIFEST.tsv',description='Release provenance / rule / scope / uncertainty schematic; no measured numerical ordinate')],[source],'schematic content and provenance',figure='FigureS'+str(i),panel='a')
add('FigS2',[row(int(r.cutoff),'U2','candidates_millions',r.u2_n/1e6,int(r.u2_n),int(r.u2_n),'11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv') for _,r in arms.iterrows()],['11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv'],'cutoff;u2_n',figure='FigureS2',panel='a')
rew=rd('14_repair/REWIRE_CONTROL_METRICS.tsv');src='14_repair/REWIRE_CONTROL_METRICS.tsv'
add('FigS6',[row(int(r.cutoff),r.condition,'full_boundary_tie_fde',r.full_boundary_tie_fde,int(r.n_candidates),int(r.actual_evaluated_n),src,n_positives=int(r.n_positives),base_rate=r.base_rate,selected_hits=r.full_boundary_tie_hits) for _,r in rew.iterrows()],[src,'14_repair/REWIRE_DEGREE_AUDIT.tsv'],'cutoff;condition;full_boundary_tie_fde;n_candidates;actual_evaluated_n;n_positives;base_rate;full_boundary_tie_hits',figure='FigureS6',panel='a')
add('FigS8',globalrows(top[top.cutoff.eq(2018)&top.subset.eq('ALL_U2')]),[gsrc],'cutoff;model;full_boundary_tie_fde;n_candidates;n_positives;actual_evaluated_n',figure='FigureS8',panel='a')
src='14_repair/FINGERPRINT_BIN_COUNTS.tsv'
add('FigS9',[row(int(r.cutoff),r.stratum,'future_positive_count',int(r.n_positive),int(r.n_candidates),int(r.n_positive),src) for _,r in bins.iterrows()],[src],'cutoff;stratum;n_positive;n_candidates',figure='FigureS9',panel='a')
add('FigS10a',globalrows(top[top.subset.eq('NEW')&top.model.isin([lig,deg])]),[gsrc],'cutoff;model;subset;full_boundary_tie_fde;n_candidates;n_positives;actual_evaluated_n',figure='FigureS10',panel='a')
add('FigS10b',compoundrows(tie[tie.subset.eq('NEW')&tie.model.isin([lig,deg,tar])],'expected_recall_at_1'),[tsrc],'cutoff;model;subset;expected_recall_at_1;n_positive_compounds',figure='FigureS10',panel='b')
# Figure-level sheets additionally satisfy the journal convention, retaining requested panel sheets.
for i in range(1,6):
 selected=[p for p in payload if p['name'].startswith('Fig'+str(i))]
 frames=[]
 for p in selected:
  d=pd.DataFrame(p['matrix'][1:],columns=p['matrix'][0]);d.insert(0,'panel',p['name']);frames.append(d)
 d=pd.concat(frames,ignore_index=True).fillna('NOT_APPLICABLE');payload.append({'name':'Figure'+str(i),'matrix':[d.columns.tolist()]+d.values.tolist()})
payload=json.loads(json.dumps(payload,default=lambda x:x.item()))
(OUT/'source_data_payload.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2))
pd.DataFrame(manifest).to_csv(OUT/'SOURCE_DATA_MANIFEST.tsv',sep='\t',index=False)
for p in payload:pd.DataFrame(p['matrix'][1:],columns=p['matrix'][0]).to_csv(OUT/(p['name']+'.tsv'),sep='\t',index=False)
print('Source Data:',len(payload),'sheets;',len(manifest),'provenance links')
