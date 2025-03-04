import scanpy as sc 
import numpy as np

adata = sc.read_h5ad("/home/sadigungor/Desktop/pathwayScoring/test/test_data/pbmc3k.h5ad")

print(adata)

print(adata.n_obs)
