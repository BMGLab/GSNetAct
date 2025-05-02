from gsnetact import createAdataObject

import pandas as pd

adata_ = "./test_data/pbmc3k.h5ad"
# The path to the h5ad file.

jsonFile = "/home/sadigungor/Desktop/testPathway/840.json"
# The path to the JSON file.

_adata = createAdataObject(adata_, jsonFile, normalized=True)
# Call the createObject function from the package.

df = pd.DataFrame(_adata.X)
# Create a pandas dataframe from the AnnData object's X layer.
df.columns = _adata.var
# Set the column names to the geneset names that are located in the var layer.

df.to_csv("output.csv", sep="\t")
# Create an output file, named  as output.csv.
