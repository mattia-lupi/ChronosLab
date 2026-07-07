#pragma once

#include "lapacke.h"
#include "precision.h"

void computeFirstQR(double *Ahat, lapack_int sizeI, lapack_int sizeJ, double *R, 
                    double *Rtriang, double *tau, double *work, lapack_int lwork, 
                    lapack_int &info);

void applyFirstQt(double *Ahat, lapack_int sizeI, lapack_int sizeJ, double *tau, 
                  double *a0k, double *work, lapack_int lwork, lapack_int &info);

void applyR(lapack_int sizeJ, double *R, double *a0k, lapack_int &info);

void applyQt(iReg t, lapack_int *sizeJ, lapack_int *sizeI, lapack_int *qStart, 
             double *Ahat, double *tau, double *a0k, lapack_int nRowsRHS, 
             lapack_int ncolsRHS, double *work, lapack_int lwork, lapack_int &info);

void computeNewQR(iReg t, lapack_int *sizeI, lapack_int *sizeJ, lapack_int *qStart, 
                  double *Ahat, double *tau, double *R, double *Rtriang, double *work, 
                  lapack_int lwork, lapack_int &info);
