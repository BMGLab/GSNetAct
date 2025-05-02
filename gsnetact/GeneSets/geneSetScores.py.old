import numpy as np


class GeneSetScore(dict):

    def __init__(self, matrix, geneNamesList):

        self.matrix = matrix

        self.geneNamesList = geneNamesList

        num_rows, num_cols = matrix.shape

        row_nz_counts = np.count_nonzero(matrix, axis=1)
        # Get the non-zero element count for each row

        scores = np.zeros(num_rows)
        # Create an empty matrix to hold the scores.

        for col in range(num_cols):
            nonzero_rows = np.nonzero(matrix[:, col])[0]
            # Get nonzero rows.

            if nonzero_rows.size == 2:
                # All columns in the GeneSetMatrix object are gonna have only 2
                # elements that are different from 0,
                # control it just in case.

                i, j = nonzero_rows
                # Get the two elements in a column to i,j

                scores[i] += matrix[i, col] * row_nz_counts[j]
                # Multiply the value in the row i
                # with the nonzero count in the row j and
                # add that to the score of row i.

                scores[j] += matrix[j, col] * row_nz_counts[i]
                # Vice versa.

        for i, gene in enumerate(geneNamesList):
            # Create the dictionary
            self[gene] = scores[i]
