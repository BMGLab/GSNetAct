import json

# Load the JSON data from the file
with open('/home/sadi/Desktop/GSNetAct/test/test_data/deneme.json', 'r') as f:
    data = json.load(f)

# Create a new dictionary to store the filtered data
filtered_data = {}

# Iterate through the data and filter out values less than 0.4
for gene_set, gene_data in data.items():
    filtered_data[gene_set] = {}
    for gene, relations in gene_data.items():
        filtered_data[gene_set][gene] = {}
        for related_gene, value in relations.items():
            if value >= 0.4:
                filtered_data[gene_set][gene][related_gene] = value

# Write the filtered data to a new JSON file
with open('/home/sadi/Desktop/GSNetAct/test/test_data/deneme_filtered.json', 'w') as f:
    json.dump(filtered_data, f, indent=4)

print("Filtered data has been written to deneme_filtered.json")
