from Scoring._annData import createObject

adata_ = "/home/sadigungor/Desktop/testPathway/neutrophils_filtered.h5ad"
jsonFile = "/home/sadigungor/Desktop/testPathway/840.json"

_adata = createObject(adata_,jsonFile,normalized=True)

print(_adata.X)
print(_adata.var)
