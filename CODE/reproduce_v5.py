"""Editorial exports from frozen summary tables only; no scientific computation."""
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle,FancyArrowPatch
R=Path(__file__).resolve().parents[1];FIG=R/'REPRODUCED/figures';FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
C=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#777777']
def finish(fig,name):
 fig.tight_layout(pad=1.4)
 for ext in ['svg','png']:fig.savefig(FIG/(name+'_NCOMMS_V2.'+ext),dpi=300)
 plt.close(fig)
def overview():
 fig=plt.figure(figsize=(7.2,4.6));gs=fig.add_gridspec(2,3,height_ratios=[3,1.05],left=.035,right=.99,top=.94,bottom=.08,wspace=.16,hspace=.35)
 axs=[fig.add_subplot(gs[0,i]) for i in range(3)]
 for a in axs:a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
 def text(a,x,y,t,**kw):a.text(x,y,t,ha='center',va='center',fontsize=8,**kw)
 def arrow(a,x,y,xx,yy,color='0.3'):a.add_patch(FancyArrowPatch((x,y),(xx,yy),arrowstyle='-|>',mutation_scale=9,color=color,lw=1))
 a=axs[0];a.set_title('a  Two clocks',loc='left',fontweight='bold')
 # Two parallel time tracks with a shared cutoff, visibly different event positions.
 for y,color in [(.73,C[0]),(.36,C[1])]:arrow(a,.08,y,.94,y,color)
 a.plot([.55,.55],[.2,.88],ls='--',color='0.55',lw=.8);text(a,.55,.94,'Cutoff')
 a.plot(.25,.73,'o',color=C[0],ms=7);a.plot(.79,.36,'s',color=C[1],ms=7)
 text(a,.25,.86,'Evidence event');text(a,.28,.61,'Document /\npublication date')
 text(a,.77,.51,'Database curation');text(a,.72,.2,'Release object →\navailability proxy')
 text(a,.5,.03,'Event time ≠ availability',fontweight='bold')
 a=axs[1];a.set_title('b  Database membership',loc='left',fontweight='bold')
 a.add_patch(Circle((.35,.58),.27,facecolor=C[0],alpha=.12,edgecolor=C[0]));a.add_patch(Circle((.65,.58),.27,facecolor=C[1],alpha=.15,edgecolor=C[1]))
 # Edges, not database-size bubbles; conceptual membership only.
 for x,y,col in [(.2,.61,C[0]),(.28,.45,C[0]),(.49,.65,'0.35'),(.5,.52,'0.35'),(.72,.63,C[1]),(.79,.48,C[1])]:
  a.plot([x-.04,x+.04],[y-.03,y+.03],color=col,lw=2);a.scatter([x-.04,x+.04],[y-.03,y+.03],s=12,color=col)
 text(a,.23,.89,'Strict historical\nChEMBL object');text(a,.78,.89,'Modern ChEMBL 37\nDocument-year rollback')
 text(a,.18,.22,'Strict-only',color=C[0]);text(a,.5,.32,'Overlap');text(a,.82,.22,'Rollback-only',color='#A86000')
 text(a,.5,.03,'Membership differs',fontweight='bold')
 a=axs[2];a.set_title('c  Benchmark consequence',loc='left',fontweight='bold')
 a.plot([.41,.61],[.86,.86],color=C[0],lw=1.5);a.scatter([.41,.61],[.86,.86],s=35,color=C[0]);text(a,.5,.96,'Same candidate–target pair')
 for y,t,col in [(.71,'STRICT: absent at cutoff',C[0]),(.55,'Later appearance → future positive',C[0]),(.34,'ROLLBACK: already included', '#A86000'),(.19,'Treated as known → excluded','#A86000')]:text(a,.5,y,t,color=col)
 arrow(a,.5,.66,.5,.6,C[0]);arrow(a,.5,.29,.5,.24,C[1]);text(a,.5,.03,'LABEL ERASURE',fontweight='bold',color='#A86000')
 a=fig.add_subplot(gs[1,:]);a.axis('off');a.set(xlim=(0,1),ylim=(0,1));a.set_title('d  Historical anchors and censored follow-up',loc='left',fontweight='bold',fontsize=10)
 arrow(a,.05,.5,.96,.5)
 for x,lab,sub in [(.1,'2018 cutoff','ChEMBL 24.1'),(.3,'2019','ChEMBL 25'),(.5,'2020 cutoff','ChEMBL 27'),(.7,'2022 cutoff','ChEMBL 31'),(.93,'19 Sep 2025','ChEMBL 36 proxy\n2022 right-censored')]:
  a.plot(x,.5,'o',ms=5,color=C[0] if 'cutoff' in lab else C[1]);text(a,x,.82,lab);text(a,x,.15,sub)
 # Fixed geometry; no tight_layout needed for the deliberately assembled schematic.
 for ext in ['svg','png']:fig.savefig(FIG/('Figure1_NCOMMS_V2.'+ext),dpi=300)
 plt.close(fig)
overview()
# Reuse quantitative panel structure, read only existing frozen TSVs.
source=(R/'CODE/reproduce_figures.py').read_text()
pre=source[source.index("COL=['"):source.index('def main():')]
pre=pre[:pre.index('def finish(')]+pre[pre.index('def grouped('):]
pre=pre.replace("fontsize=9","fontsize=9")
ns=globals();ns['ROOT']=R;ns['OUT']=R/'14_repair'
exec(pre,ns)
body=source[source.index(" fig,ax=plt.subplots(2,2,figsize=(7.2,5.4))"):source.index(' # Supplementary panels')]
body=body.replace("'c  Top-tie any-hit fraction'","'c  Tie-exact MRR'").replace("'top_tie_any_hit_fraction'","'mrr_tie_expectation'")
body=body.replace("grouped(ax[1,0],tie,mods,'mrr_tie_expectation'","grouped(ax[1,0],tie,['M5_ligand_max','M1_target_degree','M2_compound_degree'],'mrr_tie_expectation'")
body=body.replace("'a  Nominal top 0.1%, full boundary tie'","'a  Global enrichment'")
body=body.replace("'U2'","'Protein-class'").replace('U2 candidates (millions)','Protein-class candidates (millions)').replace('a  U2 label erasure','a  Future-label erasure')
body=body.replace("'a  Molecular relationship classes A–G'","'a  Fingerprint-identical\\ncandidate–target relationships'").replace("'Tanimoto = 1 pairs (%)'","'Fingerprint-identical\\ncandidate–target pairs (%)'")
body=body.replace("figsize=(7.2,6)","figsize=(7.2,6.2)")
body=body.replace("finish(fig,'Figure4')","fig.text(.5,.01,'Fingerprint-identical candidates are a minority of future-positive pairs.',ha='center',fontsize=8);finish(fig,'Figure4')")
# The original panel manifests are not rewritten; new export provenance is in the V2 manifest.
ns['panel']=lambda *args,**kw:None
ns['read']=lambda name:pd.read_csv(R/'14_repair'/name,sep='\t')
ns['top']=ns['read']('TOP_PERCENTILE_TIE_AUDIT.tsv');ns['tie']=ns['read']('COMPOUND_CENTRIC_TIE_REPAIR.tsv');ns['chem']=ns['read']('CHEMICAL_CLASSIFICATION_SUMMARY.tsv');ns['effect']=ns['read']('NESTED_GRAPH_EFFECTS.tsv')
ns['back']=pd.read_csv(R/'10_phase1/CHEMBL_BACKFILL_BY_CUTOFF.tsv',sep='\t');ns['arms']=pd.read_csv(R/'11_phase2/BENCHMARK_ARM_DEFINITIONS.tsv',sep='\t');ns['labs']=ns['read']('AUTHORITATIVE_LABEL_SUMMARY.tsv')
# Source contains one-space function indentation; de-indent without executing its pipeline entrypoint.
exec('\n'.join(line[1:] if line.startswith(' ') else line for line in body.splitlines()),ns)
fig,ax=plt.subplots(figsize=(7.2,3.2));grouped(ax,ns['tie'],['M5_ligand_max','M2_compound_degree'],'top_tie_any_hit_fraction','Top-tie set coverage under tied rankings','ALL_U2');ax.set_ylim(0,1.32);ax.set_ylabel('Top-tie any-hit fraction');fig.text(.5,.01,'Coverage metric, not rank-1 discrimination; degree = 1 by construction.',ha='center',fontsize=8)
fig.tight_layout(pad=1.8)
for ext in ['svg','png']:fig.savefig(FIG/('Supplementary_Figure_TopTieAnyHit.'+ext),dpi=300)
plt.close(fig)
print('Five main figures and supplementary set-coverage figure exported from frozen summaries.')
