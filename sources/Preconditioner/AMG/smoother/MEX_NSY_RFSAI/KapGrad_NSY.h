#pragma once

#include <math.h>    // to use: fabs
#include <algorithm> // to use: min

#include "DEBUG.h"
#include "ri_sortsplit_nsy.h"
#include "heapsort.h"

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

// Compute the gradient of the Kaporin Number
void KapGrad_NSY(const int istep, const int irow, int& mrow, const int jendbloc,
                 const int lfil, const int* RESTRICT iat, const int* RESTRICT ja,
                 const double* RESTRICT coef, const double* RESTRICT coef_T,
                 const double* RESTRICT sol_L, const double* RESTRICT sol_U,
                 int* RESTRICT IWN, int* RESTRICT JWN,
                 double* RESTRICT WR_L, double* RESTRICT WR_U);
