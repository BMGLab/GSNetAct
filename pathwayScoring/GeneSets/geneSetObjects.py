import numpy as np 

from Items.MatrixItem import MatrixItem
from jsonParser import pjson

class GeneSetMatrix:

    """ Creates and returns the gene set as an incidence matrix. 

        Why an incidence matrix?

            If we think of the gene sets and the relations between genes, they are very similar to graphs. 
            The best way to represent graphs mathematically is by using an incidence matrix.

    """
    
    def __new__(cls,rawGeneSet):

        unique_identifiers = {} # The unique identifiers dictionary to hold the data of the unique position 
                                # of an edge(weight, relation whetever you wanna call it) in the matrix 
                                # for faster calculations.
        n1 = 0
        sayac = 0
        GeneList = []

        for i in rawGeneSet:
            n1 += 1

            for j in rawGeneSet[i]:
                
                unique_hash = hash(i) * hash(j) * hash(rawGeneSet[i][j]) # Creating a unique hash to determine the                                                                         
                                                                         # unique position.

                if unique_hash not in unique_identifiers: # If the edge is unique, give it a new position by giving                                                           
                                                          # it a new column value.
                    unique_identifiers[unique_hash] = sayac
                    sayac += 1
            
                GeneList.append(MatrixItem(float(rawGeneSet[i][j]),n1-1,unique_identifiers[unique_hash]-1))
                # Whether the edge is unique or not, create a MatrixItem to hold the position data and add it to
                # the GeneList.

        n = len(unique_identifiers)

        matrix = np.zeros((n1,n)) # Create the matrix. TODO: Test the sparse matrix approach.

        for item in GeneList: # Build the incidence matrix using the position data in the GeneList.
            w,r,c = item.value,item.row,item.column
            matrix[r,c] = w

        return matrix

class GeneSet:

    """ Holds the properties of a gene set. """
    
    def __init__(self,_id,_rawGeneSet):
        
        self.id = _id
        self.asJson = _rawGeneSet
        self.matrix = GeneSetMatrix(_rawGeneSet)
    
    @property
    def getID(self):

        return self.id
    
    @property
    def getMatrix(self):

        return self.matrix
    
    @property
    def getAsJson(self):

        return self.asJson

    @property
    def getGeneNames(self): # Access gene names as a list
        namesListBuffer = [i for i in self.asJson]

        return namesListBuffer

if __name__ == "__main__":
    newjson = pjson("/home/sadigungor/pathwayScoring/deneme.json")
 
    for y in newjson:
        #a = GeneSetMatrix(newjson[y])

        b = GeneSet("GeneSet1",newjson[y])

        #print(a)
        
        print(b.getAsJson)
        print(b.getMatrix)
        print(b.getID)
        print("## TESTING GENE NAMES LIST ## ")
        print(b.getNamesList)
        break 
