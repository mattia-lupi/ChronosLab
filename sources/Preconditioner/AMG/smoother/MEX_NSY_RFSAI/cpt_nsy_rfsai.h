#ifndef CPT_NSY_RFSAI_H
#define CPT_NSY_RFSAI_H

int cpt_nsy_rfsai(const int nstep, const int step_size, const double eps,
                  const int nn_A, const int nt_A, const double *diag_A,
                  const int *iat_A, const int *ja_A, const double *coef_A,
                  const double *coef_AT, int *&iat_FL, int *&ja_FL,
                  double *&coef_FL, double *&coef_FUT,
                  const int num_threads = 1);

#endif
