import numpy as np


class GeneSetScore(dict):

    def __init__(self, geneID, matrix, geneNamesList, epsilon=10**-6):

        self.epsilon = epsilon

        self.matrix = matrix

        self.geneNamesList = geneNamesList

        for gene in self.geneNamesList:
            # Create the empty dictionary
            # Which all the scores within are initially zero.
            self[gene] = 0

        num_rows, num_cols = matrix.shape

        row_nz_counts = matrix.getnnz(axis=1)
        # Get the non-zero element count for each row

        disconnected_genes = np.where(row_nz_counts == 0)[0]
        for idx in disconnected_genes:

            self[geneNamesList[idx]] += self.epsilon

        # Search for rows that contain only zeros,
        # Those are the genes that are not related to
        # any other gene. If you find one, add epsilon to
        # the score of the respected gene, which is 0.
        # This way, genes that are not in relation to others
        # can have an effect to the score.

        # Iterating over columns of a CSC matrix
        for col in range(num_cols):
            start = matrix.indptr[col]
            end = matrix.indptr[col+1]
            rows = matrix.indices[start:end]
            data = matrix.data[start:end]

            if len(rows) == 2:
                i, j = rows
                val_i, val_j = data

                self[geneNamesList[i]] += val_i * row_nz_counts[j]
                # Multiply the value in the row i
                # with the nonzero count in the row j and
                # add that to the score of row i.

                self[geneNamesList[j]] += val_j * row_nz_counts[i]
                # Vice versa.

       
