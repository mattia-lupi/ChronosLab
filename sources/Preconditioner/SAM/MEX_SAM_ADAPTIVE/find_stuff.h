#pragma once

#ifdef I
#undef I
#endif

#include "precision.h"

void fullA0k(const iExt nn_A, const iExt __restrict *iat0,                                 
             const iReg __restrict *ja0, const double __restrict *coef0,                   
             const iExt k, double *A0k);

void findNonZeroInColJ(const iReg __restrict *J, const iExt __restrict *iatk,              
                       const iReg __restrict *jak, const iReg n2, iReg *I, iReg &sizeI);

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, iExt *iat0, iReg *ja0, 
            double *coef0, iExt k);

void getAhat( iReg *I, iReg sizeI, iReg *J, iReg Jstart, iReg Jend,
             iExt *iatk, iReg *jak, double *coefk, double *Ahat, iExt &Astart);

void getAJ(iReg *J, iReg Jstart, iReg Jend, iExt nn_A, iExt *iatk, iReg *jak, 
           double *coefk, iExt *jatAJ, iReg *iaAJ, double *coefAJ);

void fillL(iReg *__restrict L, const double *__restrict res, iExt nn_A, iReg &usedL);

void findJtilde(iReg *Jtilde, iReg &JtildeSize, const iReg __restrict *L,                  
                const iReg sizeL, const iExt __restrict *iatk, const iReg __restrict *jak, 
                iReg *J, iReg sizeJ);

void getAJtilde(iExt nn_A, iExt *iatk, iReg *jak, double *coefk, iReg *Jtilde,
                iReg JtildeSize, iExt* jatAJtilde, iReg *iaAJtilde, double *coefAJtilde, 
                iExt *workspace);
