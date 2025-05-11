# 2brie_tutorial 

# script voor het volgen van de brie2 tutorial 

#laden van packages 
import brie
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
scv.logging.print_version()

scv.settings.verbosity = 3               # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True       # set max width size for presenter view
scv.set_figure_params('scvelo')           # for beautified visualization

# define the path you store the example data
# dat_dir = "./raw_data"
dat_dir = '/home/anne.brussaard/Project-BRIE2/brie2_tutorial/raw_data/scNTseq/'

#scVelo dynamical model with total RNAs
adata = scv.read(dat_dir + "/neuron_splicing_totalRNA.h5ad")

scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


#vanaf hier een foutmelding
scv.tl.recover_dynamics(adata, var_names=adata.var_names[:200])





