import numpy as np

from ..Utils.MatrixItem import MatrixItem


class GeneSetMatrix:

    """ Creates and returns the gene set as an incidence matrix. 


        If we think of the gene sets and the relations between genes, they are 
            
        very similar to graphs. 

        The best way to represent graphs mathematically is by using an 
            
        incidence matrix. With incidence matrixes,
            
        our mathematical formula can be solved 
            
        for each geneset expeditiously.

        
        The main idea is to create an incidence matrix with specific positions 
        
        for each gene to make calculations
        
        efficiently. In order to do this, 
        
        we give each gene a specific column value.

    """
    
    def __new__(cls, rawGeneSet):

        unique_identifiers = {}
        # The unique identifiers dictionary to hold the data of the unique
        # position of an edge(weight, relation whetever you wanna call it)
        # in the matrix for faster calculations.

        position_in_column = 0
        GeneList = []

        for geneCount, gene in enumerate(rawGeneSet):

            for relatedGene in rawGeneSet[gene]:

                unique_hash = hash(gene) * hash(relatedGene) * hash(rawGeneSet[gene][relatedGene])
                # Creating a unique hash to determine the
                # unique position.

                if unique_hash not in unique_identifiers:
                    # If the edge is unique, give it a new position by giving
                    # it a new column value.

                    unique_identifiers[unique_hash] = position_in_column
                    position_in_column += 1

                GeneList.append(MatrixItem(float(rawGeneSet[gene][relatedGene]),
                                           geneCount-1,
                                           unique_identifiers[unique_hash]-1))
                # Whether the edge is unique or not, create a MatrixItem
                # to hold the position data and add it to
                # the GeneList.

        n = len(unique_identifiers)

        matrix = np.zeros((geneCount, n))
        # Create the matrix. TODO: Test the sparse matrix approach.

        for item in GeneList:
            # Build the incidence matrix using
            # the position data in the GeneList.

            w, r, c = item.weight, item.row, item.column
            matrix[r, c] = w

        return matrix


class GeneSet:

    """ Holds the properties of a gene set. """

    def __init__(self, _id, _rawGeneSet):

        self.id = _id
        self.asJson = _rawGeneSet
        self.matrix = GeneSetMatrix(_rawGeneSet)

    @property
    def getID(self):

        return self.id

    @property
    def getAsJson(self):

        return self.asJson

    @property
    def getGeneNames(self):
        # Access gene names as a list

        namesListBuffer = [i for i in self.asJson]

        return namesListBuffer


def getGSNA(jsonFile):
    # Automatically create geneset objects and return them
    # in a list.

    _list = []

    for i in jsonFile:

        newGeneSet = GeneSet(f"{i}", jsonFile[i])

        _list.append(newGeneSet)

    return _list
