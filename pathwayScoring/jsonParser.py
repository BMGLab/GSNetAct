import json 

class pjson(dict):

    def __init__(self,jsonFile):

        self.jsonFile = jsonFile
        self.geneSayi = []
        self.geneNamesList = []

        with open(jsonFile) as f: #Load the json file
            self.js = json.load(f)

        self.update(self.js) # Update the object so that it represents the content of the given JSON file.

        for i in self.js.values():
            
            for j in i.keys():
                self.geneNamesList.append(j) # Get the gene names
            
            self.geneSayi.append(len(i)) # Get the gene counts 

    def __repr__(self):

        return repr(self.js)

    def getGeneNames(self): # Acces the gene names
        return self.geneNamesList

    def getGeneCounts(self): # Access the gene counts
        return self.geneSayi

    def getGeneSetCount(self): # Access the gene set count
        return len(self.geneSayi)

    def getFileInfo(self): # Access the file info as a string
        info = ""

        for i,j in enumerate(self.geneSayi):
            info += f"GeneSet : {i+1}, Gene Count : {j} \n"

        return info

    def getFileName(self):
        return self.jsonFile

if __name__ == "__main__":
    
    newjson = pjson("/home/sadigungor/pathwayScoring/deneme.json")
    
    #print(newjson.getFileInfo())

    #for i in newjson:
    #   print(i,"yeah")

    print(newjson.getGeneNames())
