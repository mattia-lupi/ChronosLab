#include "transpose.h"
#include "merge_row_patt.h"

// Symmetrizes the pattern of a non-symmetric matrix by padding with zeroes missing enries
int SymmetrizePattern(const int nrows, int *& iat, int *& ja, double *& coef,
                      double *&coef_T);
