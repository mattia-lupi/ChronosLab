#pragma once

#include "precision.h"

void computeFirstQR(double *Ahat, iReg sizeI, iReg sizeJ, double *R, double *Rtriang, double *tau, double *work, iExt lwork, iReg &info);

void applyFirstQt(double *Ahat, iReg sizeI, iReg sizeJ, double *tau, double *a0k, double *work, iExt lwork, iReg &info);

void applyR(iReg sizeJ, double *R, double *a0k, iReg &info);

void applyQt(iReg t, iReg *sizeJ, iReg *sizeI, iExt *qStart, double *Ahat, 
             double *tau, double *a0k, iReg nrowsA0k, iReg ncolsA0k, double *work, iExt lwork, iReg &info);

void computeNewQR(iReg t, iReg *sizeI, iReg *sizeJ, iExt *qStart, double *Ahat, double *tau, double *R, 
                  double *Rtriang, double *work, iExt lwork, iReg &info);
