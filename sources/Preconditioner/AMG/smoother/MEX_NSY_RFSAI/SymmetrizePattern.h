#pragma once

// Symmetrizes the pattern of a non-symmetric matrix by padding with zeroes missing enries
int SymmetrizePattern(const int nrows, int *& iat, int *& ja, double *& coef,
                      double *&coef_T, const int num_threads);
