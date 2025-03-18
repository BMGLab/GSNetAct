from Scoring._annData import createObject

adata_ = "./test_data/pbmc3k.h5ad"
jsonFile = "./test_data/big_genesets_relations.json"

_adata = createObject(adata_,jsonFile,normalized=True)

print(_adata)
