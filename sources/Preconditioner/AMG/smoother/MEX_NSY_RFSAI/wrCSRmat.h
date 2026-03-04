#pragma once
#include <cstdio>

void wrCSRmat(FILE *ofile, const bool patt, const int nn_A, const int *const iat_A,
              const int *const ja_A, const double *const coef_A);
