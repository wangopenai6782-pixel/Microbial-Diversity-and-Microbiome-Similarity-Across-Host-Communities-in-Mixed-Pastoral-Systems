# Author: Wang Linglong
# Public release script derived from: PCoA_host_species_composition_output\Mantel.txt
# Purpose: reproduce the corresponding analysis/figure in the manuscript.
#
# Inputs should be placed under data/public_figures/ unless paths are supplied
# explicitly in the script. Outputs are written to results/public_figures/ when
# the script uses output_dir/public_output helpers.
import pandas as pd, numpy as np, os, math
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from pathlib import Path

input_dir = Path(os.environ.get('JEV_FIGURE_INPUT_DIR', Path(__file__).resolve().parents[2] / 'data' / 'public_figures'))
jfile = input_dir / 'Host_Jaccard_similarity_matrix.csv'
gfile = input_dir / '最佳拟合模型_修改后gam.csv'
outdir = Path(os.environ.get('JEV_FIGURE_OUTPUT_DIR', Path(__file__).resolve().parents[2] / 'results' / 'public_figures')) / 'Mantel_test_results'
outdir.mkdir(parents=True, exist_ok=True)

jac=pd.read_csv(jfile, index_col=0)
# ensure numeric
jac=jac.apply(pd.to_numeric)
hosts=list(jac.index)

# read gam
gam=pd.read_csv(gfile)
name_map={
 'Bos':'Cattle','Gallus':'Chicken','Hen Harrier':'Circus cyaneus','Circus':'Circus cyaneus','Dauricus':'Daurian ground squirrel',
 'Egret':'Egret','Eagle owl':'Eurasian eagle owl','Sparrow':'Eurasian tree sparrow','Great Bustard':'Great bustard','Homo':'Homo','Mus':'House mouse','mosquito':'Mosquito','Mosquito':'Mosquito','Sus':'Pig','Ory':'Rabbit','Pigeon':'Rock pigeon','Ovis':'Sheep','Ixo':'Tick','Tick':'Tick'
}
gam['h1']=gam['species_1'].map(lambda x: name_map.get(str(x).strip(), str(x).strip()))
gam['h2']=gam['species_2'].map(lambda x: name_map.get(str(x).strip(), str(x).strip()))
print('mapped unique', sorted(set(gam.h1).union(gam.h2)))
missing=set(gam.h1).union(gam.h2)-set(hosts)
print('missing vs hosts', missing)

def overlap_to_score(x):
    x=str(x).strip().lower()
    if x=='complete': return 1.0
    if x=='partial': return 0.5
    if x=='distinct': return 0.0
    return np.nan

gam['habitat_score']=gam['habitat_overlap'].map(overlap_to_score)
gam['diet_score']=gam['diet_overlap'].map(overlap_to_score)
gam['ecological_similarity']=gam[['habitat_score','diet_score']].mean(axis=1)
gam['ecological_distance']=1-gam['ecological_similarity']

# build matrices
n=len(hosts)
microbial=1-jac.values
np.fill_diagonal(microbial,0)
eco=np.full((n,n), np.nan); phylo=np.full((n,n), np.nan); np.fill_diagonal(eco,0); np.fill_diagonal(phylo,0)
idx={h:i for i,h in enumerate(hosts)}
for _,r in gam.iterrows():
    h1,h2=r['h1'],r['h2']
    if h1 in idx and h2 in idx:
        i,j=idx[h1],idx[h2]
        eco[i,j]=eco[j,i]=r['ecological_distance']
        phylo[i,j]=phylo[j,i]=r['phylogenetic_distance']

print('na eco', np.isnan(eco).sum(), 'na phylo', np.isnan(phylo).sum())
if np.isnan(eco).any() or np.isnan(phylo).any():
    for matname,mat in [('eco',eco),('phylo',phylo)]:
        wh=np.argwhere(np.isnan(mat))
        print(matname, [(hosts[i],hosts[j]) for i,j in wh[:20]])

# vectors upper triangle
tri=np.triu_indices(n,1)
vec_microbial=microbial[tri]
vec_eco=eco[tri]
vec_phylo=phylo[tri]

def spearman(x,y):
    return spearmanr(x,y).correlation

def partial_spearman(x,y,z):
    rxy=spearman(x,y); rxz=spearman(x,z); ryz=spearman(y,z)
    den=math.sqrt((1-rxz*rxz)*(1-ryz*ryz))
    return (rxy-rxz*ryz)/den if den>0 else np.nan

def mantel_perm(xmat, yvec, perms=9999, seed=123):
    rng=np.random.default_rng(seed)
    obs=spearman(xmat[tri], yvec)
    count=1
    perm_stats=[]
    for k in range(perms):
        p=rng.permutation(n)
        xv=xmat[p][:,p][tri]
        r=spearman(xv,yvec)
        perm_stats.append(r)
        if abs(r) >= abs(obs): count+=1
    pval=count/(perms+1)
    return obs,pval,np.array(perm_stats)

def partial_mantel_perm(xmat, yvec, zvec, perms=9999, seed=124):
    rng=np.random.default_rng(seed)
    obs=partial_spearman(xmat[tri], yvec, zvec)
    count=1
    perm_stats=[]
    for k in range(perms):
        p=rng.permutation(n)
        xv=xmat[p][:,p][tri]
        r=partial_spearman(xv,yvec,zvec)
        perm_stats.append(r)
        if abs(r) >= abs(obs): count+=1
    pval=count/(perms+1)
    return obs,pval,np.array(perm_stats)

res=[]
obs,pval,_=mantel_perm(microbial, vec_eco, seed=1)
res.append(['Microbial distance vs ecological distance','Mantel',obs,pval,9999])
obs,pval,_=mantel_perm(microbial, vec_phylo, seed=2)
res.append(['Microbial distance vs phylogenetic distance','Mantel',obs,pval,9999])
obs,pval,_=partial_mantel_perm(microbial, vec_eco, vec_phylo, seed=3)
res.append(['Microbial distance vs ecological distance controlling for phylogenetic distance','Partial Mantel',obs,pval,9999])
obs,pval,_=partial_mantel_perm(microbial, vec_phylo, vec_eco, seed=4)
res.append(['Microbial distance vs phylogenetic distance controlling for ecological distance','Partial Mantel',obs,pval,9999])

results=pd.DataFrame(res, columns=['Test','Method','Mantel_r','P_value','Permutations'])
print(results)
results.to_csv(outdir/'Mantel_partial_Mantel_results.csv', index=False)

pd.DataFrame(microbial,index=hosts,columns=hosts).to_csv(outdir/'Microbial_dissimilarity_matrix_1_minus_Jaccard.csv')
pd.DataFrame(eco,index=hosts,columns=hosts).to_csv(outdir/'Ecological_distance_matrix.csv')
pd.DataFrame(phylo,index=hosts,columns=hosts).to_csv(outdir/'Phylogenetic_distance_matrix.csv')
# merged pairs
pairs=[]
for i,j in zip(*tri):
    pairs.append([hosts[i],hosts[j],jac.values[i,j],microbial[i,j],eco[i,j],phylo[i,j]])
pairdf=pd.DataFrame(pairs, columns=['host1','host2','Jaccard_similarity','microbial_distance_1_minus_Jaccard','ecological_distance','phylogenetic_distance'])
pairdf.to_csv(outdir/'Mantel_plot_data.csv', index=False)

# plot
fig, axs=plt.subplots(1,2, figsize=(8,4), dpi=160)
for ax, x, label in zip(axs,[vec_eco,vec_phylo],['Ecological distance','Phylogenetic distance']):
    ax.scatter(x, vec_microbial, s=24, alpha=0.75)
    # trendline
    m,b=np.polyfit(x, vec_microbial, 1)
    xx=np.linspace(np.min(x), np.max(x), 100)
    ax.plot(xx, m*xx+b, lw=1.5)
    ax.set_xlabel(label)
    ax.set_ylabel('Microbial community dissimilarity\n(1 - Jaccard similarity)')
    ax.spines[['top','right']].set_visible(False)
plt.tight_layout()
fig.savefig(outdir/'Mantel_distance_relationships.png', dpi=600)
fig.savefig(outdir/'Mantel_distance_relationships.pdf')
print('saved to', outdir)

