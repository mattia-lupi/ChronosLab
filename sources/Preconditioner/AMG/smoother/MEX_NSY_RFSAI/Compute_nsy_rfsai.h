#pragma once

int Compute_nsy_rfsai(const int nstep, const int step_size, const double eps,
                      const int nn_A, const int *iat_A, const int *ja_A, 
                      const double *coef_A,
                      int *&iat_FL, int *&ja_FL, double *&coef_FL,
                      int *&iat_FU, int *&ja_FU, double *&coef_FU,
                      const int num_threads);

