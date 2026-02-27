#include "DEBUG.h"

void gather_fullsys(const int irow, const int mrow, const int *const vecinc,
                    const int nn, const int *const iat, const int *const ja,
                    const double *const coef, double *full_A, double *rhs_L,
                    double *rhs_U, bool &null_L, bool &null_U);
