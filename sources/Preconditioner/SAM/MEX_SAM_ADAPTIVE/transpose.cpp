#include "transpose.h"
#include <cstring>

/**
 * Parallel transposition of a sparse matrix with a symmetric sparsity pattern.
 *
 * @param nrows    Number of rows in the matrix.
 * @param iat      Row pointers of the input matrix (size: nrows + 1).
 * @param ja       Column indices of the input matrix (size: iat[nrows]).
 * @param coef     Matrix coefficients of the input matrix (size: iat[nrows]).
 * @param iat_T    Output row pointers (size: nrows + 1).
 * @param ja_T     Output column indices (size: iat[nrows]).
 * @param coef_T   Output matrix coefficients (size: iat[nrows]).
 * @return         0 on success, -1 on error.
 */

int transpose(const int nrows, const int *const iat, const int *const ja,
              const double *const coef, int *iat_T, int *ja_T, double *coef_T) {

    const int nnz = iat[nrows];

    // Structural symmetry means iat and ja map perfectly to iat_T and ja_T.
    std::memcpy(iat_T, iat, (nrows + 1) * sizeof(int));
    std::memcpy(ja_T,  ja,  nnz * sizeof(int));

    // Loop over the rows
    for (int i = 0; i < nrows; ++i) {
        for (int k = iat[i]; k < iat[i+1]; ++k) {
            const int j = ja[k];

            // Binary search row j for column index i to find the destination pointer
            const int* row_j_start = ja + iat[j];
            const int* row_j_end   = ja + iat[j+1];
            const int* match       = std::lower_bound(row_j_start, row_j_end, i);

            if (match != row_j_end && *match == i) {
                const int target_pos = match - ja;
                coef_T[target_pos] = coef[k];
            }
        }
    }

    return 0;
}
