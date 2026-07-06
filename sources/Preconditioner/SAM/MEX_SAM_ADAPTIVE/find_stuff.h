#pragma once

#include "precision.h"

void fullA0k(iExt nn_A, iExt *iat0, iReg *ja0, double *coef0, iExt k, double *A0k);

void findNonZeroInColJ(iReg *J, iExt *iatk, iReg *jak, iReg n2, iReg *I, iReg &sizeI);

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, iExt *iat0, iReg *ja0, double *coef0, iExt k);

void getAhat( iReg *I, iReg sizeI, iReg *J, iReg Jstart, iReg Jend,
             iExt *iatk, iReg *jak, double *coefk, double *Ahat, iExt &Astart);

void getAJ(iReg *J, iReg Jsize, iExt nn_A, iExt *iatk, iReg *jak, double *coefk, double *AJ);

void fillL(iReg *L, double *res, iExt nn_A, iReg &usedL);

void findJtilde(iReg *Jtilde, iReg &JtildeSize, iReg *L, iReg sizeL, iExt *iatk, iReg *jak, iReg *J, iReg sizeJ);

void fullAJtilde(iExt nn_A, iExt *iatk, iReg *jak, double *coefk, iReg *Jtilde, iReg JtildeSize, double *AJtilde);
