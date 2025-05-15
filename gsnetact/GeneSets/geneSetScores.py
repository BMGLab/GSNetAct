import numpy as np


class GeneSetScore(dict):

    def __init__(self, matrix, geneNamesList, epsilon=10**-6):

        self.epsilon = epsilon

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

        for rowidx, row in enumerate(matrix):
            # Search for rows that contain only zeros,
            # Those are the genes that are not related to
            # any other gene. If you find one, add epsilon to
            # the score of the respected gene, which is 0.

            # This way, genes that are not in relation to others
            # can have an effect to the score.

            if np.count_nonzero(row) == 0:
                self[geneNamesList[rowidx]] += self.epsilon

             
