from gsnetact import runGSNA, pjson
from anndata import read_h5ad
import pandas as pd


##############################################################
# from gsnetact import makeJson
# makeJson("The Path To Your GeneSets file.")
##############################################################
# ↑↑↑↑↑ Test makeJson module ↑↑↑↑↑

adata_ = read_h5ad("./test_data/pbmc3k.h5ad")
#TODO : adata'nin her birine 10^-6 gibi bir deger eklenecek. adata.X sparse mi degil mi kontrol eden bir
# kondisyon yaz
# Read the anndata object.

jsonFile = pjson("./test_data/big_genesets_relations.json")
# Parse the json file into a pjson() object.

gsnaObject = runGSNA(adata_, jsonFile, normalized=True)
# Call the createObject function from the package.

df = pd.DataFrame(gsnaObject.X)
# Create a pandas dataframe from the AnnData object's X layer.
df.columns = gsnaObject.var
# Set the column names to the geneset names that are located in the var layer.

df.to_csv("output.csv", sep="\t")
# Create an output file, named  as output.csv.
