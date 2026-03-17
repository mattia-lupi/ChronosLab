#include <iostream>  // to use: cout,endl
#include <algorithm> // to use: fill_n
#include <lapacke.h>  // to use: dgetrf,dgetrs

#include "DEBUG.h"
#include "KapGrad_NSY.h"
#include "gather_fullsys.h"
#include "inl_blas1.h"

int cpt_nsy_rfsai(const int nstep, const int step_size, const double eps, const int nn_A,
                  const int nt_A, const double *diag_A, const int *iat_A, const int *ja_A,
                  const double *coef_A, const double *coef_AT, int *&iat_FL, int *&ja_FL,
                  double *&coef_FL, double *&coef_FUT);
