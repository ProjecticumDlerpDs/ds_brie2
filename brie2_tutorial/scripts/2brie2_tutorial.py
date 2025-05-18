# 2brie_tutorial 

# script voor het volgen van de brie2 tutorial 

#laden van packages 
import umap
import os
import brie
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
print(brie.__version__)

# define the path you store the example data
# dat_dir = "./raw_data"
dat_dir = '/home/anne.brussaard/Project-BRIE2/brie2_tutorial/raw_data/msEAE/'

#BRIE2 optie 1: differntiale splicing events
#voor kwantificatie van PSI in downstream analyse
adata = sc.read_h5ad(dat_dir + "/brie_quant_cell.h5ad")
adata

adata.uns['brie_param']

#Gene index veranderen van Ensemebl id naar gene name 
adata.var['old_index'] = adata.var.index
new_index = [adata.var.GeneName[i] + adata.var.GeneID[i][18:] for i in range(adata.shape[1])]
adata.var.index = new_index

#bij index is bytes, veranderen naar str
adata.var.index = [x.decode('utf-8') if type(x) == bytes else x for x in adata.var.index]
adata.obs.index = [x.decode('utf-8') if type(x) == bytes else x for x in adata.obs.index]

#cell annotaties van input covariates
adata.obs['MS'] = ['EAE' if x else 'Ctrl' for x in adata.obsm['Xc'][:, 0]]
adata.obs['isCD1'] = ['CD1_%d' %x for x in adata.obsm['Xc'][:, 1]]

#volcano plot maken voor differntiale splicing events in de terminal 
import matplotlib.pyplot as plt

brie.pl.volcano(adata, y='ELBO_gain', log_y=False, n_anno=16, score_red=7, adjust=True)
plt.xlabel('cell_coeff: effect size on logit(Psi)')
plt.title("MS differential splicing")
plt.savefig("volcano_plot.png", dpi=300, bbox_inches='tight')

#volcano plot maken voor differntial splicing voor in de Rmarkdown 
brie.pl.volcano(adata, y='ELBO_gain', log_y=False, n_anno=16, score_red=7, adjust=True)
plt.xlabel('cell_coeff: effect size on logit(Psi)')
plt.title("MS differential splicing")

DSEs = adata.var.index[adata.varm['ELBO_gain'][:, 0] >= 7]
len(DSEs), len(np.unique([x.split('.')[0] for x in DSEs])), adata.shape

adata.uns['Xc_ids']

#visualiseren van ruwe counts voor DSEs
#voor in de Rmarkdown 
rank_idx = np.argsort(adata.varm['ELBO_gain'][:, 0])[::-1]

fig = plt.figure(figsize=(15, 8))
brie.pl.counts(adata, genes=list(adata.var.index[rank_idx[:6]]),
               color='MS', add_val='ELBO_gain', ncol=3, alpha=0.7)
# fig.savefig(dat_dir + '../../figures/fig_s8_counts.png', dpi=300, bbox_inches='tight')
plt.show()

#visualiseren in de terminal van ruwe counts voor DSEs
# sorteren op ELBO_gain 
rank_idx = np.argsort(adata.varm['ELBO_gain'][:, 0])[::-1]
#figuur maken en de counts plotten 
fig = plt.figure(figsize=(15, 8))
brie.pl.counts(
    adata,
    genes=list(adata.var.index[rank_idx[:6]]),
    color='MS',
    add_val='ELBO_gain',
    ncol=3,
    alpha=0.7
)

#figuur opslaan als afbeelding 
plt.savefig("fig_s8_counts.png", dpi=300, bbox_inches='tight')

#ratio veranderen tussen ambiguous reads
#voor in de Rmarkdown 
fig = plt.figure(figsize=(5, 4))
brie.pl.counts(adata, genes='Ddhd1',
               color='MS', add_val='ELBO_gain',
               layers=['isoform1', 'ambiguous'],
               nrow=1, alpha=0.7)
plt.show()

#voor in de terminal
fig = plt.figure(figsize=(5, 4))
brie.pl.counts(
    adata,
    genes='Ddhd1',
    color='MS',
    add_val='ELBO_gain',
    layers=['isoform1', 'ambiguous'],
    nrow=1,
    alpha=0.7
)
plt.savefig("counts_Ddhd1.png", dpi=300, bbox_inches='tight')

#BRIE2 optie 2: splicing kwantificatie en gebruik
#aggregatie selecteren 
adata_aggr = sc.read_h5ad(dat_dir + "/brie_quant_aggr.h5ad")
adata_aggr

#gene index veranderen 
print(np.mean(adata.var['old_index'] == adata_aggr.var.index))
adata_aggr.var.index = adata.var.index

#van bytes naar str veranderen 
adata_aggr.var.index = [x.decode('utf-8') if type(x) == bytes else x
                        for x in adata_aggr.var.index]
adata_aggr.obs.index = [x.decode('utf-8') if type(x) == bytes else x
                        for x in adata_aggr.obs.index]
                        
# meta data en gen level annotaties toevoegen. 
print(np.mean(adata.obs.index == adata_aggr.obs.index))
adata_aggr.obs['MS'] = adata.obs['MS'].copy()
adata_aggr.obs['isCD1'] = adata.obs['isCD1'].copy()

dat_umap = np.genfromtxt(dat_dir + '/cell_X_umap.tsv', dtype='str', delimiter='\t')

mm = brie.match(adata_aggr.obs.index, dat_umap[:, 0])
idx = mm[mm != None].astype(int)

adata_aggr = adata_aggr[mm != None, :]
adata_aggr.obsm['X_GEX_UMAP'] = dat_umap[idx, 1:3].astype(float)
adata_aggr.obs['cluster'] = dat_umap[idx, 3]
adata_aggr.obs['combine'] = [adata_aggr.obs['cluster'][i] + '-' + adata_aggr.obs['MS'][i]
                            for i in range(adata_aggr.shape[0])]
adata_aggr

#filteren van de cellen 
#Rmarkdown 
plt.hist(np.log10(adata_aggr.X.sum(1)[:, 0] + 1), bins=30)
plt.xlabel("log10(total reads)")
plt.ylabel("Cell frequency")
plt.show()

#terminal 
plt.hist(np.log10(adata_aggr.X.sum(1)[:, 0] + 1), bins=30)
plt.xlabel("log10(total reads)")
plt.ylabel("Cell frequency")
plt.savefig("hist_total_reads.png", dpi=300, bbox_inches='tight')

min_reads = 3000

adata_aggr = adata_aggr[adata_aggr.X.sum(axis=1) > min_reads, :]
adata = adata[adata_aggr.obs.index, :]

#visualisatie splicing phenotypen in gen expressie umap 
#Rmarkdown 
sc.pl.scatter(adata_aggr, basis='GEX_UMAP', color=['cluster', 'MS'], size=30)

#Terminal
sc.pl.scatter(
    adata_aggr,
    basis='GEX_UMAP',
    color=['cluster', 'MS'],
    size=30,
    show=False,        
    save='_scatter.png'  
)

#vanwege de grote hoeveelheid afbeeldingen die voor deze tutorial worden gemaakt is besloten de rest van de tutorial alleen uit te werken via de Rmarkdown. 
