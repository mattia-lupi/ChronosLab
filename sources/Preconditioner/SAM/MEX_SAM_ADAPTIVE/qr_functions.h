#pragma once

#include "lapack.h"
#include "precision.h"

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void computeFirstQR(double *Ahat, lapack_int sizeI, lapack_int sizeJ, double *R, 
                    double *Rtriang, double *tau, double *work, lapack_int lwork, 
                    lapack_int &info);

void applyFirstQt(double *Ahat, lapack_int sizeI, lapack_int sizeJ, double *tau, 
                  double *a0k, double *work, lapack_int lwork, lapack_int &info);

void applyR(lapack_int sizeJ, double *R, double *a0k, lapack_int &info);

void applyQt(iReg t, const lapack_int* RESTRICT sizeJ,
             const lapack_int* RESTRICT sizeI, lapack_int *qStart,
             double *Ahat, double *tau, double *a0k, lapack_int nRowsRHS,
             lapack_int ncolsRHS, double *work, lapack_int lwork, lapack_int &info);

void computeNewQR(iReg t, lapack_int *sizeI, lapack_int *sizeJ, lapack_int *qStart, 
                  double *Ahat, double *tau, double *R, double *Rtriang, double *work, 
                  lapack_int lwork, lapack_int &info);
