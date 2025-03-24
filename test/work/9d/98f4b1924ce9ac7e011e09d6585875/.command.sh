#!/bin/bash -ue
python3 -c '
from Scoring._annData import createObject

import pandas as pd

_adata = createObject("pbmc3k.h5ad","big_genesets_relations.json",normalized=True)

df = pd.DataFrame(_adata.X)
df.columns = _adata.var

df.to_csv("output.csv",sep="	")


'
