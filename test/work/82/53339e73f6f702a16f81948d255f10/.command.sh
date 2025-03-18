#!/bin/bash -ue
python3 -c '
from Scoring._annData import createObject 
print(createObject("pbmc3k.h5ad","big_genesets_relations.json",normalized=True).X)'
