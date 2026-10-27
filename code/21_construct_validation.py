"""
21_construct_validation.py — baseline construct validation of the stability-flexibility classification
(Methods "Baseline construct validation"; Results "Construct validation"; Figure S1; Table S6).

Uses the baseline session only (before randomization). One index per construct and task, each aligned so that
higher values mean greater stability (flexibility-pole indices sign-reversed):
  PRL             perseveration (S), mean inverse temperature (S), mean learning rate (F)
  rule switching  pre-switch accuracy (S) and post-switch accuracy (F), each averaged over the two phases
  demand selection task-switch cost in RT (S) and in accuracy (S), averaged over the two phases
                  (sessions whose learning phase was at chance are excluded)
  MSIT            RT interference effect (F)
Test: mean Spearman rho over the 23 cross-task pairs, against 10,000 permutations in which participants are
shuffled independently within each index (missing values stay in place). Within-task pairs are shown but not
tested. Test-retest reliability: Spearman rho between baseline and retest after removing arm means.
Also printed: correlation of each index with a general-accuracy composite (mean z of overall accuracy in the
four tasks).

Output: results/tableS6_construct_pairs.csv, results/tableS6_reliability.csv, figures/FigureS1_construct_validation.png,
        results/21_construct_validation.txt
"""
import os, sys, pandas as pd, numpy as np, warnings
from scipy import stats
from config import DERIVED, RESULTS, FIGURES
warnings.filterwarnings('ignore')
_log = open(os.path.join(RESULTS, "21_construct_validation.txt"), "w")
_print = print
def print(*a, **k):
    _print(*a, **k); _print(*a, file=_log, **k)
rd = lambda f: pd.read_csv(os.path.join(DERIVED, f))
d = rd('allocation.csv'); d['SubjectID'] = d.StudyID.str.replace(r'\D', '', regex=True).astype(int)
d = d[d.alloc.isin(['Placebo', 'Psilocybin'])]
grp = d.set_index('SubjectID').alloc
L = rd('all_measures_long.csv'); L = L[L.SubjectID.isin(set(d.SubjectID))]
# all_measures_long.csv already has demand selection measures of chance-level baseline sessions set to missing
def idx_frame(S):
    S=S.copy(); m=lambda cols:S[cols].mean(axis=1)
    S['PRL_beta']=m(['PRL_betaGain','PRL_betaLoss']); S['PRL_alpha']=m(['PRL_alphaGain','PRL_alphaLoss'])
    S['CF_Pre']=m(['CF_PreSwitchAcc_learning','CF_PreSwitchAcc_transfer']); S['CF_Post']=m(['CF_PostSwitchAcc_learning','CF_PostSwitchAcc_transfer'])
    S['ES_SC_RT']=m(['ES_Learn_SwitchCost_RT','ES_Dec_SwitchCost_RT']); S['ES_SC_Acc']=m(['ES_Learn_SwitchCost_Acc','ES_Dec_SwitchCost_Acc'])
    return S
B=idx_frame(L[L.Session=='pre'].set_index('SubjectID')); P=idx_frame(L[L.Session=='post'].set_index('SubjectID'))
idx={'PRL_Perseveration':1,'PRL_beta':1,'PRL_alpha':-1,'CF_Pre':1,'CF_Post':-1,'ES_SC_RT':1,'ES_SC_Acc':1,'MSIT_InterferenceEffect':-1}
task={'PRL_Perseveration':'Reward learning','PRL_beta':'Reward learning','PRL_alpha':'Reward learning','CF_Pre':'Rule switching','CF_Post':'Rule switching','ES_SC_RT':'Cued switching','ES_SC_Acc':'Cued switching','MSIT_InterferenceEffect':'MSIT'}
lab={'PRL_Perseveration':'Perseveration (S)','PRL_beta':'Inverse temperature (S)','PRL_alpha':'Learning rate (F, reversed)','CF_Pre':'Pre-switch accuracy (S)','CF_Post':'Post-switch accuracy (F, reversed)','ES_SC_RT':'Cued switch cost, RT (S)','ES_SC_Acc':'Cued switch cost, accuracy (S)','MSIT_InterferenceEffect':'MSIT interference, RT (F, reversed)'}
A=pd.DataFrame({k:v*B[k] for k,v in idx.items()})
keys=list(idx); pairs=[(a,b) for i,a in enumerate(keys) for b in keys[i+1:] if task[a]!=task[b]]
def meanr(X):
    R=X.corr(method='spearman'); return np.mean([R.loc[a,b] for a,b in pairs])
obs=meanr(A); rng=np.random.default_rng(2026); perm=[]
for _ in range(10000):
    Q=A.copy()
    for c in Q:
        v=Q[c].values.copy(); msk=~np.isnan(v); v[msk]=rng.permutation(v[msk]); Q[c]=v
    perm.append(meanr(Q))
perm=np.array(perm); p_perm=np.mean(np.abs(perm)>=abs(obs))
rows=[]
for a in keys:
    for b in keys:
        if keys.index(b)<=keys.index(a): continue
        x=A[[a,b]].dropna(); r,p=stats.spearmanr(x[a],x[b]); rows.append(dict(a=a,b=b,cross=task[a]!=task[b],n=len(x),rho=r,p=p))
T=pd.DataFrame(rows)
rel=[]
for c in keys:
    x=pd.concat([B[c].rename('a'),P[c].rename('b'),grp.rename('g')],axis=1).dropna()
    r,p=stats.spearmanr(x.a-x.groupby('g').a.transform('mean'), x.b-x.groupby('g').b.transform('mean')); rel.append(dict(index=c,n=len(x),rho=r,p=p))
REL=pd.DataFrame(rel)
gacc=pd.concat([(B[c]-B[c].mean())/B[c].std() for c in ['PRL_RewardAcc','CF_Accuracy','ES_Learn_Acc','MSIT_InterferenceAcc']],axis=1).mean(axis=1)
G=A.corrwith(gacc,method='spearman')
print(f'mean cross-task aligned rho={obs:.3f}; permutation p={p_perm:.4f} (10000); positive {int((T[T.cross].rho>0).sum())}/{int(T.cross.sum())}; complete cases {A.dropna().shape[0]}')
print(T.round(3).to_string()); print(REL.round(2).to_string()); print('correlation with general accuracy composite:',G.round(2).to_dict())
T.to_csv(os.path.join(RESULTS,'tableS6_construct_pairs.csv'),index=False); REL.to_csv(os.path.join(RESULTS,'tableS6_reliability.csv'),index=False)
# ---- figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'text.color':'#0b0b0b','axes.labelcolor':'#52514e','xtick.color':'#52514e','ytick.color':'#52514e'})
cmap=LinearSegmentedColormap.from_list('div',['#2a78d6','#f0efec','#e34948'])
fig=plt.figure(figsize=(14,6.6))
ax=fig.add_axes([0.15,0.25,0.38,0.62]); cax=fig.add_axes([0.548,0.52,0.009,0.33]); ax2=fig.add_axes([0.77,0.25,0.14,0.62])
R=A.corr(method='spearman'); k=len(keys)
for a_ in (ax,ax2): a_.set_facecolor('#fcfcfb')
for i in range(k):
    for j in range(k):
        if j>=i: continue
        r=R.iloc[i,j]; a,b=keys[j],keys[i]; row=T[((T.a==a)&(T.b==b))].iloc[0]
        cross=task[a]!=task[b]
        ax.add_patch(Rectangle((j,i),1,1,facecolor=cmap((r+1)/2) if cross else '#fcfcfb',edgecolor='#fcfcfb',lw=2))
        if not cross: ax.add_patch(Rectangle((j+0.06,i+0.06),0.88,0.88,facecolor='none',edgecolor='#9a998f',lw=0.9,ls=(0,(2,2))))
        txt=f'{r:+.2f}'.replace('0.','.')
        ax.text(j+.5,i+.42,txt,ha='center',va='center',fontsize=9.5,fontweight='bold' if row.p<.05 else 'normal',color='#0b0b0b' if cross else '#52514e')
        ax.text(j+.5,i+.76,f'n={row.n}',ha='center',va='center',fontsize=6.5,color='#52514e')
ax.set_xlim(0,k-1); ax.set_ylim(k,1)
ax.set_xticks(np.arange(k-1)+.5); ax.set_xticklabels([lab[c] for c in keys[:-1]],rotation=35,ha='right',fontsize=8)
ax.set_yticks(np.arange(1,k)+.5); ax.set_yticklabels([lab[c] for c in keys[1:]],fontsize=8)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.tick_params(length=0)
fig.text(0.01,0.955,'A  Baseline correlations between indices, each aligned so that higher = more stable',fontsize=10.5,fontweight='bold')
fig.text(0.028,0.92,f'Mean cross-task $r_s$ = {obs:+.2f} (permutation p = {p_perm:.3f}); {int((T[T.cross].rho>0).sum())} of {int(T.cross.sum())} cross-task pairs positive. S, stability pole; F, flexibility pole (sign reversed).',fontsize=8.5,color='#52514e')
sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(-1,1)); cb=fig.colorbar(sm,cax=cax,ticks=[-1,-.5,0,.5,1]); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7,length=0); cb.set_label('Spearman $r_s$',fontsize=8)
fig.text(0.545,0.47,'Bold: p < .05\nDashed outline:\nwithin-task pair,\nnot in the test',fontsize=7.5,color='#52514e',va='top')
y=np.arange(k)
ax2.axvline(0,color='#b8b7b1',lw=1)
ax2.hlines(y,0,REL.rho,color='#2a78d6',lw=2)
ax2.scatter(REL.rho,y,s=45,color='#2a78d6',edgecolor='#fcfcfb',linewidth=2,zorder=3)
for yi,(r,n) in enumerate(zip(REL.rho,REL.n)): ax2.text(1.12,yi,f'{r:+.2f} (n={n})'.replace('0.','.'),ha='left',va='center',fontsize=7.5,color='#52514e')
ax2.set_yticks(y); ax2.set_yticklabels([lab[c] for c in keys],fontsize=8); ax2.set_ylim(k-0.5,-0.5); ax2.set_xlim(-1,1); ax2.set_xticks([-1,-.5,0,.5,1])
ax2.set_xlabel('Baseline vs post-dose $r_s$, within arm',fontsize=8)
for s_ in ['top','right']: ax2.spines[s_].set_visible(False)
fig.text(0.605,0.955,'B  Test–retest reliability',fontsize=10.5,fontweight='bold')
fig.text(0.623,0.92,'Low reliability limits cross-task correlations',fontsize=8.5,color='#52514e')
fig.savefig(os.path.join(FIGURES,'FigureS1_construct_validation.png'),dpi=300,facecolor='#fcfcfb')
print('saved')
