import numpy as np


class GeneSetScore(dict):

    def __init__(self, matrix, geneNamesList):

        self.matrix = matrix

        self.geneNamesList = geneNamesList

        for gene in self.geneNamesList:
            # Create the empty dictionary
            # Which all the scores within are initially zero.
            self[gene] = 0

        num_rows, num_cols = matrix.shape

        row_nz_counts = np.count_nonzero(matrix, axis=1)
        # Get the non-zero element count for each row

        for col in range(num_cols):
            nonzero_rows = np.nonzero(matrix[:, col])[0]
            # Get nonzero rows.

            if nonzero_rows.size == 2:
                # All columns in the GeneSetMatrix object are gonna have only 2
                # elements that are different from 0,
                # control it just in case.

                i, j = nonzero_rows
                # Get the two elements in a column to i,j

                self[geneNamesList[i]] += matrix[i, col] * row_nz_counts[j]
                # Multiply the value in the row i
                # with the nonzero count in the row j and
                # add that to the score of row i.

                self[geneNamesList[j]] += matrix[j, col] * row_nz_counts[i]
                # Vice versa.
