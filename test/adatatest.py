from gsnetact import runGSNA, pjson
import scanpy as sc
import pandas as pd


##############################################################
# from gsnetact import makeJson
# makeJson("The Path To Your GeneSets file.")
##############################################################
# ↑↑↑↑↑ Test makeJson module ↑↑↑↑↑

adata_ = sc.read_h5ad("/home/sadigungor/Desktop/GSNetAct/test/test_data/pbmc3k.h5ad")
# The path to the h5ad file.

jsonFile = pjson("/home/sadigungor/Desktop/GSNetAct/test/test_data/deneme.json")
# The path to the JSON file.

gsnaObject = runGSNA(adata_, jsonFile, normalized=True)
# Call the createObject function from the package.

df = pd.DataFrame(gsnaObject.X)
# Create a pandas dataframe from the AnnData object's X layer.
df.columns = gsnaObject.var
# Set the column names to the geneset names that are located in the var layer.

df.to_csv("output.csv", sep="\t")
# Create an output file, named  as output.csv.
