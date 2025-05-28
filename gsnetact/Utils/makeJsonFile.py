from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import requests
import json
import gc


def create_session():
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=0.2, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    return session


def fetch_filtered_string_network(genes,
                                  session,
                                  species="9606",
                                  required_score=100,
                                  limit=1000):

    api_url = "https://string-db.org/api/json/network"
    data = {
        "identifiers": "\n".join(genes),
        "species": species,
        "required_score": required_score,
        "limit": limit,
        "network_type": "functional"
    }
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    
    try:
        response = session.post(api_url, data=data, headers=headers)
        response.raise_for_status()
        network_data = response.json()
        if not network_data:
            print(f"Warning! The Server Returned Null. Genes: {genes}")
            return None
        df = pd.DataFrame(network_data)
        filtered_df = df[(df['preferredName_A'].isin(genes)) & (df['preferredName_B'].isin(genes))]
        return filtered_df.to_dict(orient="records")
    except requests.exceptions.RequestException as e:
        print(f"Data could not be pulled error : {e}")
        return None


def process_gene_set(index, gene_set_name, genes, session):
    print(f"Iteration {index + 1}: Gene Set {gene_set_name}...")
    result = fetch_filtered_string_network(genes, session, limit=len(genes))
    gc.collect()
    return gene_set_name, result


def fetch_all_networks_parallel(gene_sets, max_workers=25):
    results = {}

    with create_session() as session:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_gene_set, index, gene_set_name, genes['geneSymbols'], session): gene_set_name
                for index, (gene_set_name, genes) in enumerate(gene_sets.items())
            }

            for future in futures:
                gene_set_name, result = future.result()
                results[gene_set_name] = result

    return results

# TODO: FARKLI DOSYA BICIMLERI EKLE


def makeJson(msigdbFile, jsonFileName="geneSets.json"):
    with open(msigdbFile, 'r') as f:
        gene_sets = json.load(f)
        all_networks = fetch_all_networks_parallel(gene_sets, max_workers=40)
        relation_dict = {}

    for gene_set_id, interactions in all_networks.items():
        setBuffer = set()
        for interaction in interactions:
            node1 = interaction['preferredName_A']
            node2 = interaction['preferredName_B']
            combined_score = interaction['score']

            setBuffer.add(node1)
            setBuffer.add(node2)

            if gene_set_id not in relation_dict:
                relation_dict[gene_set_id] = {}

            if node1 not in relation_dict[gene_set_id]:
                relation_dict[gene_set_id][node1] = {}
            relation_dict[gene_set_id][node1][node2] = combined_score

            if node2 not in relation_dict[gene_set_id]:
                relation_dict[gene_set_id][node2] = {}
            relation_dict[gene_set_id][node2][node1] = combined_score

        for i in gene_sets[gene_set_id]['geneSymbols']:
            if i not in setBuffer:
                relation_dict[gene_set_id][i] = {}

    with open(jsonFileName, "w") as f:
        json.dump(relation_dict, f, indent=2)
