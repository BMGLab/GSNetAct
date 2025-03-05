import scanpy as sc 
import numpy as np

adata = sc.read_h5ad("/home/sadigungor/Desktop/pathway_scorers/test/test_data/neutrophils_filtered.h5ad")

print("########### PRINT VAR ##########")
print(adata.var)
print("########### PRINT X #############")
print(adata.X)

