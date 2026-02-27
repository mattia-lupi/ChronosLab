#pragma once

int EXTI_prolongation(const int level, const int nthreads,
                      const int *vecstart, const int nn_A, const int nt_A,
                      const int *const iat_A, const int *const ja_A,
                      const double *const coef_A, const int *const coef_S,
                      const int *const iat_C, const int *const ja_C,
                      const double *const coef_C, const int *const fcnodes,
                      const int nr_I, const int nc_I, int &nt_I, int *&iat_I,
                      int *&ja_I, double *&coef_I);
