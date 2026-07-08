#include <iostream>
#include <algorithm> // to use: fill_n

int transpose(const int nrows, const int ncols, const int *const iat, const int *const ja,
              const double *const coef, int *&iat_T, int *&ja_T, double *&coef_T);
