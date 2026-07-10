#pragma once

#ifdef I
#undef I
#endif

#include "precision.h"


#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif


void fullA0k(const iExt nn_A, const iExt* RESTRICT iat0,
             const iReg* RESTRICT ja0, const double* RESTRICT coef0,
             const iExt k, double *A0k);

void findNonZeroInColJ(const iReg* RESTRICT J, const iExt* RESTRICT iatk,
                       const iReg* RESTRICT jak, const iReg n2, iReg *I, iReg &sizeI);

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, iExt *iat0, iReg *ja0,
            double *coef0, iExt k);

void getAhat( iReg *I, iReg sizeI, iReg *J, iReg Jstart, iReg Jend,
             iExt *iatk, iReg *jak, double *coefk, double *Ahat, iExt &Astart);

void getAJ(iReg *J, iReg Jstart, iReg Jend, iExt *iatk, iReg *jak,
           double *coefk, iExt *jatAJ, iReg *iaAJ, double *coefAJ);

void fillL(iReg * RESTRICT L, const double * RESTRICT res, iExt nn_A, iReg &usedL);

void findJtilde(iReg *Jtilde, iReg &JtildeSize, const iReg* RESTRICT L,
                const iReg sizeL, const iExt* RESTRICT iatk, const iReg* RESTRICT jak,
                iReg *J, iReg sizeJ);

void getAJtilde(iExt nn_A, iExt *iatk, iReg *jak, double *coefk, iReg *Jtilde,
                iReg JtildeSize, iExt* jatAJtilde, iReg *iaAJtilde, double *coefAJtilde,
                iExt *workspace);
