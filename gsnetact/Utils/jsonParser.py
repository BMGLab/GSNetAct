import json


class pjson(dict):

    def __init__(self, jsonFile):

        self.jsonFile = jsonFile
        self.geneSayi = []
        self.geneNamesList = []
        self.iter_count = 0

        with open(jsonFile) as f:
            # Load the json file
            self.js = json.load(f)

        for i in self.js.values():

            for j in i.keys():
                self.geneNamesList.append(j)
                # Get the gene names

            self.geneSayi.append(len(i))
            # Get the gene counts

        self.update(self.js)

    def genesets(self):
        # Return genesets in the JSON file,
        # since they are represented in self.values().
        # Used for easier iteration.

        return self.values()

    @property
    def getAsDict(self):
        return self.js

    @property
    def getGeneNames(self):
        # Acces the gene names
        return self.geneNamesList

    @property
    def getGeneCounts(self):
        # Access the gene counts
        return self.geneSayi

    @property
    def getGeneSetCount(self):
        # Access the gene set count
        return len(self.geneSayi)

    @property
    def getFileInfo(self):
        # Access the file info as a string
        info = ""

        for i, j in enumerate(self.geneSayi):
            info += f"GeneSet : {i+1}, Gene Count : {j} \n"

        return info

    def getFileName(self):
        return self.jsonFile
