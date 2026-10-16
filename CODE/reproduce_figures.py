"""Content-first V2 figures. Python/matplotlib only; SVG and QA PNG, no PDF.

Contract: database counts unchanged; panels distinguish metric definitions,
UNKNOWN chemistry, and addition/removal graph interventions. Point estimates
are descriptive on one frozen seed, not repeated-run means or significance.
Quantitative grids use actual repaired TSV rows and explicit sample sizes.
"""
from repair_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
import shutil
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'14_repair'
RESULT=Path(__file__).resolve().parents[1]/'REPRODUCED'
RESULT.mkdir(exist_ok=True)
FIG=RESULT/'figures';FIG.mkdir(exist_ok=True)
def save(rows,name):
 df=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
 df.to_csv(RESULT/name,sep='\t',index=False);return df
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'svg.fonttype':'path','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
COL=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#777777']
YEARS=[2018,2020,2022];YLAB=['2018','2020','2022\ncensored']
LABEL={'M5_ligand_max':'Ligand similarity','M2_compound_degree':'Compound degree','M1_target_degree':'Target degree','M9_lightgcn':'LightGCN','M7_svd':'SVD','M3_product':'Preferential attachment'}
MODEL_COLOR={'M5_ligand_max':COL[0],'M2_compound_degree':COL[1],'M1_target_degree':COL[2],'M7_svd':COL[3],'M9_lightgcn':COL[4],'M3_product':COL[5]}
manifest=[]
def read(name):return pd.read_csv(OUT/name,sep='\t',dtype={'chembl_version':str})
def panel(fig,p,src,models,metric,subset='ALL_U2',years=YEARS):
 src=';'.join(x if x.startswith(('10_phase1/','11_phase2/','12_phase3/','14_repair/')) else '14_repair/'+x for x in src.split(';'))
 manifest.append(dict(figure=fig,panel=p,source_table=src,models_shown='|'.join(models),metric=metric,years='|'.join(map(str,years)),subset=subset,script='scripts/phase4_figures.py',svg=f'14_repair/figures/{fig}.svg',uncertainty='NONE: descriptive fixed-data point estimate, one frozen training seed; no significance claim'))
def finish(fig,name):
 fig.tight_layout(pad=1.4);fig.savefig(FIG/f'{name}.svg');fig.savefig(FIG/f'{name}.png',dpi=150);plt.close(fig)
def grouped(ax,df,models,metric,title,subset=None):
 if subset is not None:df=df[df.subset.eq(subset)]
 x=np.arange(3);w=.8/len(models)
 for i,m in enumerate(models):
  vals=[float(df[(df.cutoff.eq(y))&df.model.eq(m)][metric].iloc[0]) for y in YEARS]
  ax.bar(x+(i-(len(models)-1)/2)*w,vals,w,label=LABEL.get(m,m),color=MODEL_COLOR.get(m,COL[i]))
 ax.set_xticks(x,YLAB);ax.set_title(title,loc='left');ax.legend(frameon=False)
 if metric=='full_boundary_tie_fde':ax.set_ylim(0,float(df[df.model.isin(models)][metric].max())*1.35)
def note(ax,title,text):
 ax.axis('off');ax.set_title(title,loc='left');ax.text(.03,.85,text,va='top',transform=ax.transAxes,fontsize=9,linespacing=1.5)
def main():
 top=read('TOP_PERCENTILE_TIE_AUDIT.tsv');tie=read('COMPOUND_CENTRIC_TIE_REPAIR.tsv');chem=read('CHEMICAL_CLASSIFICATION_SUMMARY.tsv');effect=read('NESTED_GRAPH_EFFECTS.tsv')
 back=pd.read_csv(ROOT/'10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv',sep='\t');arms=pd.read_csv(ROOT/'11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv',sep='\t');labs=read('AUTHORITATIVE_LABEL_SUMMARY.tsv')
 fig,ax=plt.subplots(2,2,figsize=(7.2,5.2))
 for a,p,title,text in zip(ax.flat,'abcd',['a  Two clocks','b  Historical objects','c  Release anchors','d  Label rule'],['Document year: evidence event\nOperational availability: release object time\nFTP mtime is a proxy, not proof of first access','Strict: recovered historical release\nRollback: ChEMBL 37 filtered by document year','2018: ChEMBL 24.1\n2020: ChEMBL 27\n2022: ChEMBL 31\nChEMBL 25 is after the 2018 cutoff','First observed in recovered release sequence\nNo prior recovered appearance\nGap reappearance excluded\n2022 follow-up ends 2025-09-19']):
  note(a,title,text);panel('Figure1',p,'AUTHORITATIVE_LABEL_SUMMARY.tsv;CHEMBL_RELEASE_DATE_PROVENANCE.tsv',[],'DESIGN',years=YEARS)
 finish(fig,'Figure1')
 fig,ax=plt.subplots(2,2,figsize=(7.2,5.4));x=np.arange(3)
 bottom=np.zeros(3)
 for col,label,color in [('N_overlap','Overlap',COL[0]),('N_strict_only','Strict only',COL[4]),('N_retro_only','Rollback only',COL[1])]:ax[0,0].bar(x,back[col],bottom=bottom,label=label,color=color);bottom+=back[col].to_numpy()
 ax[0,0].set_ylim(0,float(bottom.max())*1.25);ax[0,0].set_xticks(x,YLAB);ax[0,0].legend(frameon=False);ax[0,0].set_title('a  Edge-set composition',loc='left')
 ax[0,1].bar(x,100*back.backfill_fraction,color=COL[1]);ax[0,1].set_xticks(x,YLAB);ax[0,1].set_ylabel('Rollback-only edges (%)');ax[0,1].set_title('b  Reconstruction mismatch',loc='left')
 ax[1,0].bar(x-.18,100*arms.label_erasure_fraction_of_u2_positives,.36,label='U2',color=COL[0]);ax[1,0].bar(x+.18,100*labs.historical_active_erased/labs.historical_active_positive,.36,label='Historical-active',color='0.65');ax[1,0].set_xticks(x,YLAB);ax[1,0].legend(frameon=False);ax[1,0].set_title('c  Two erasure denominators',loc='left');ax[1,0].set_ylabel('Future labels erased (%)');ax[1,0].set_ylim(0,105)
 ax[1,1].bar(x,arms.u2_n/1e6,color=COL[0]);ax[1,1].set_xticks(x,YLAB);ax[1,1].set_ylabel('U2 candidates (millions)');ax[1,1].set_title('d  Strict ranking universe',loc='left')
 for p,src,met in [('a','10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv','edge counts'),('b','10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv','backfill_fraction'),('c','AUTHORITATIVE_LABEL_SUMMARY.tsv;11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv','label erasure'),('d','11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv','u2_n')]:panel('Figure2',p,src,[],met)
 finish(fig,'Figure2')
 fig,ax=plt.subplots(2,2,figsize=(7.2,5.6));mods=['M5_ligand_max','M2_compound_degree']
 grouped(ax[0,0],top,mods,'full_boundary_tie_fde','a  Nominal top 0.1%, full boundary tie','ALL_U2');ax[0,0].set_ylabel('Enrichment')
 grouped(ax[0,1],tie,['M5_ligand_max','M1_target_degree','M2_compound_degree'],'expected_recall_at_1','b  Expected Recall@1','ALL_U2');ax[0,1].set_ylim(0,1)
 grouped(ax[1,0],tie,mods,'top_tie_any_hit_fraction','c  Top-tie any-hit fraction','ALL_U2');ax[1,0].set_ylim(0,1.07)
 four=['M5_ligand_max','M2_compound_degree','M7_svd','M9_lightgcn'];vals=[top[top.cutoff.eq(2022)&top.subset.eq('ALL_U2')&top.model.eq(m)].full_boundary_tie_fde.iloc[0] for m in four]
 ax[1,1].bar(range(4),vals,color=[MODEL_COLOR[m] for m in four]);ax[1,1].set_xticks(range(4),['Ligand','Compound\ndegree','SVD','LightGCN']);ax[1,1].set_ylabel('Enrichment');ax[1,1].set_title('d  2022 censored replication',loc='left')
 for p,src,mm,metric,yrs in [('a','TOP_PERCENTILE_TIE_AUDIT.tsv',mods,'full_boundary_tie_fde',YEARS),('b','COMPOUND_CENTRIC_TIE_REPAIR.tsv',['M5_ligand_max','M1_target_degree','M2_compound_degree'],'expected_recall_at_1',YEARS),('c','COMPOUND_CENTRIC_TIE_REPAIR.tsv',mods,'top_tie_any_hit_fraction',YEARS),('d','TOP_PERCENTILE_TIE_AUDIT.tsv',four,'full_boundary_tie_fde',[2022])]:panel('Figure3',p,src,mm,metric,years=yrs)
 finish(fig,'Figure3')
 bins=read('FINGERPRINT_BIN_COUNTS.tsv')
 fig,ax=plt.subplots(2,2,figsize=(7.2,6));bottom=np.zeros(3);classes=chem.chemical_class.drop_duplicates().tolist()
 for i,cls in enumerate(classes):
  v=[chem[chem.cutoff.eq(y)&chem.chemical_class.eq(cls)].fraction.iloc[0]*100 for y in YEARS]
  ax[0,0].bar(range(3),v,bottom=bottom,color=COL[i],label=cls[0]);bottom+=v
 ax[0,0].set_xticks(range(3),YLAB);ax[0,0].legend(ncol=7,frameon=False,loc='upper left',bbox_to_anchor=(0,-.22),columnspacing=.7,handlelength=1);ax[0,0].set_title('a  Molecular relationship classes A–G',loc='left');ax[0,0].set_ylabel('Tanimoto = 1 pairs (%)')
 p=bins[bins.cutoff.eq(2018)];ax[0,1].bar(range(5),p.n_positive,color=COL[0]);ax[0,1].set_xticks(range(5),['=1','.95–1','.8–.95','.5–.8','<.5']);ax[0,1].set_title('b  2018 future-positive counts',loc='left')
 subs=['ALL_U2','TANIMOTO_LT_0.95','TANIMOTO_LT_0.80','NEW']
 for i,y in enumerate((2018,2020)):
  vals=[top[top.cutoff.eq(y)&top.model.eq('M5_ligand_max')&top.subset.eq(sub)].full_boundary_tie_fde.iloc[0] for sub in subs]
  ax[1,0].plot(range(4),vals,'o-',color=COL[i],label=str(y))
 ax[1,0].set_xticks(range(4),['All','<.95','<.80','New\nscaffold']);ax[1,0].set_ylabel('Enrichment');ax[1,0].set_title('c  Ligand similarity after exclusions',loc='left');ax[1,0].legend(frameon=False)
 grouped(ax[1,1],top,mods,'full_boundary_tie_fde','d  Valid NEW scaffold only','NEW');ax[1,1].set_ylabel('Enrichment')
 for p,src,mm,metric,sub,yrs in [('a','CHEMICAL_CLASSIFICATION_SUMMARY.tsv',[],'fraction','TANIMOTO_1',YEARS),('b','FINGERPRINT_BIN_COUNTS.tsv',[],'n_positive','ALL_U2',[2018]),('c','TOP_PERCENTILE_TIE_AUDIT.tsv',['M5_ligand_max'],'full_boundary_tie_fde','ALL_U2|<.95|<.80|NEW',[2018,2020]),('d','TOP_PERCENTILE_TIE_AUDIT.tsv',mods,'full_boundary_tie_fde','NEW',YEARS)]:panel('Figure4',p,src,mm,metric,sub,yrs)
 finish(fig,'Figure4')
 fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));ax[0,0].bar(range(3),100*arms.label_erasure_fraction_of_u2_positives,color=COL[1]);ax[0,0].set_xticks(range(3),YLAB);ax[0,0].set_title('a  U2 label erasure',loc='left');ax[0,0].set_ylabel('Percent')
 models=['M9_lightgcn','M7_svd','M3_product']
 for j,y in enumerate((2018,2020)):
  for i,effectname in enumerate(('ADD','REMOVE','FULL')):
   a=ax.flat[1+i];g=effect[effect.cutoff.eq(y)&effect.effect.eq(effectname)]
   v=[g[g.model.eq(m)].delta_full_boundary_tie_fde.iloc[0] for m in models]
   a.bar(np.arange(3)+(j-.5)*.34,v,.34,label=str(y),color=COL[j]);a.axhline(0,color='0.5',lw=.6);a.set_xticks(range(3),['LightGCN','SVD','PA']);a.set_title(['b  Pure addition Δ','c  Removal Δ','d  Full reconstruction Δ'][i],loc='left');a.set_ylabel('Change in enrichment');a.legend(frameon=False)
 panel('Figure5','a','11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv',[],'label_erasure_fraction')
 for p,e in zip('bcd',['ADD','REMOVE','FULL']):panel('Figure5',p,'NESTED_GRAPH_EFFECTS.tsv',models,'delta_full_boundary_tie_fde',e,[2018,2020])
 finish(fig,'Figure5')
 # Supplementary panels preserve their original question while replacing invalid statistics.
 for number in range(1,11):
  fig,ax=plt.subplots(figsize=(7.2,3.8));name=f'FigureS{number}';src='';metric='';models=[];subset='ALL_U2'
  if number==1:
   d=read('CHEMBL_RELEASE_DATE_PROVENANCE.tsv');labels=d.chembl_version.astype(str).tolist();xx=np.arange(len(d));ax.plot(xx,np.zeros(len(d)),'o',color=COL[0]);ax.set_xticks(xx,labels);ax.set_yticks([]);ax.set_xlabel('Release identifier (date sources in table)');ax.set_title('Operational release sequence; first public access remains conditional');src='CHEMBL_RELEASE_DATE_PROVENANCE.tsv';metric='release sequence'
  elif number==2:
   ax.bar(range(3),arms.u2_n/1e6,color=COL[0]);ax.set_xticks(range(3),YLAB);ax.set_ylabel('U2 candidates (millions)');ax.set_title('Class cap 300 candidate universe');src='11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv';metric='u2_n'
  elif number==3:
   note(ax,'Confidence-score rule','Primary frozen graph uses score 9 and direct assays.\nNo confidence histogram was frozen.\nFallback rule inspected future-event counts.\nThis choice is not claimed pre-outcome.');src='ANALYSIS_DECISION_TIMELINE.tsv';metric='rule and provenance'
  elif number==4:
   note(ax,'Class-cap sensitivity scope','Cap 500 reproduced cap 300 in the frozen analysis.\nCap 100 was a sensitivity analysis.\nLegacy top-1,000 comparisons are not carried forward\nas tie-corrected model evidence.');src='11_phase2/U2_SENSITIVITY.tsv';metric='candidate cap definition'
  elif number==5:
   note(ax,'Observation windows','Authoritative builder uses recovered first appearance.\nGap reappearances are not discoveries.\nUncorrected legacy 24/48-month presence counts\nare not current endpoints. 2022 ends 2025-09-19.');src='AUTHORITATIVE_LABEL_SUMMARY.tsv';metric='followup_end'
  elif number==6:
   d=read('REWIRE_CONTROL_METRICS.tsv');x=np.arange(2)
   for j,c in enumerate(('STRICT','EXACT_REWIRE')):ax.bar(x+(j-.5)*.34,[d[d.cutoff.eq(y)&d.condition.eq(c)].full_boundary_tie_fde.iloc[0] for y in (2018,2020)],.34,label=c,color=COL[j])
   ax.set_xticks(x,['2018','2020']);ax.set_ylabel('Enrichment');ax.legend(frameon=False);ax.set_title('LightGCN under exact degree-preserving edge swaps');src='REWIRE_CONTROL_METRICS.tsv;REWIRE_DEGREE_AUDIT.tsv';metric='full_boundary_tie_fde';models=['M9_lightgcn']
  elif number==7:
   note(ax,'Uncertainty boundary','New metrics describe fixed data and one frozen model seed.\nNo new cluster interval is inferred from absent bounds.\nLegacy top-1,000 intervals are not used\nto rank repaired methods or claim significance.');src='NESTED_GRAPH_ARM_MANIFEST.tsv';metric='uncertainty statement'
  elif number==8:
   d=top[top.cutoff.eq(2018)&top.subset.eq('ALL_U2')];ax.barh(range(len(d)),d.full_boundary_tie_fde,color=COL[0]);ax.set_yticks(range(len(d)),[LABEL.get(m,m) for m in d.model]);ax.set_xlabel('Enrichment');ax.set_title('All frozen models, 2018, full boundary tie');src='TOP_PERCENTILE_TIE_AUDIT.tsv';metric='full_boundary_tie_fde';models=d.model.tolist()
  elif number==9:
   for i,y in enumerate(YEARS):ax.plot(range(5),bins[bins.cutoff.eq(y)].n_positive,'o-',label=(str(y)+' censored' if y==2022 else str(y)),color=COL[i])
   ax.set_xticks(range(5),['=1','.95–1','.8–.95','.5–.8','<.5']);ax.set_ylabel('Future positives');ax.legend(frameon=False);src='FINGERPRINT_BIN_COUNTS.tsv';metric='n_positive'
  elif number==10:
   fig.clear();ax=fig.subplots(1,2)
   grouped(ax[0],top,mods,'full_boundary_tie_fde','a  NEW scaffold global enrichment','NEW')
   grouped(ax[1],tie,['M5_ligand_max','M1_target_degree','M2_compound_degree'],'expected_recall_at_1','b  NEW scaffold Expected Recall@1','NEW');ax[1].set_ylim(0,1)
   for p,src0,ms,met in [('a','TOP_PERCENTILE_TIE_AUDIT.tsv',mods,'full_boundary_tie_fde'),('b','COMPOUND_CENTRIC_TIE_REPAIR.tsv',['M5_ligand_max','M1_target_degree','M2_compound_degree'],'expected_recall_at_1')]:panel(name,p,src0,ms,met,'NEW')
  if number!=10:panel(name,'a',src,models,metric,subset,years=([2018,2020] if number in (6,7) else [2018] if number==8 else YEARS))
  finish(fig,name)
 save(manifest,'FIGURE_CONTENT_MANIFEST.tsv')
 # Every source table is linked in the machine-readable panel manifest; actual N in top table.
 fig,axes=plt.subplots(5,3,figsize=(18,19))
 for ax,name in zip(axes.flat,[f'Figure{i}' for i in range(1,6)]+[f'FigureS{i}' for i in range(1,11)]):ax.imshow(mpimg.imread(FIG/f'{name}.png'));ax.axis('off');ax.set_title(name)
 fig.tight_layout();fig.savefig(FIG/'QA_contact_sheet.png',dpi=100);plt.close(fig)
if __name__=='__main__':main()
