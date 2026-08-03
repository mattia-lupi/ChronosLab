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

double fullA0k(const iExt * RESTRICT iat0,
             const iReg * RESTRICT ja0, const double * RESTRICT coef0,
             const iExt k, double * RESTRICT A0k,
             iReg * RESTRICT A0k_idx, iReg &A0k_nnz);

void findNonZeroInColJ(const iReg *RESTRICT J, const iExt *RESTRICT iatk,
                       const iReg *RESTRICT jak, const iReg n2, 
                       iReg *RESTRICT I, iReg &sizeI,
                       iReg *RESTRICT visited, const int t);

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, double* fullA0k);

void getAhat(iReg * RESTRICT I, iReg sizeI, iReg * RESTRICT J, iReg Jstart, 
             iReg Jend, iExt * RESTRICT iatk, iReg * RESTRICT jak, 
             double * RESTRICT coefk, double * RESTRICT Ahat, iReg &Astart);

void getAJ(iReg *J, iReg Jstart, iReg Jend, iExt *jatk, iReg *iak, double *coefk,
           const iReg **iaAJ, const double **coefAJ, iExt *jatAJ);

void findJtilde(iReg *Jtilde, iReg &JtildeSize,                                            
                const iReg* RESTRICT L, const iReg sizeL,                                  
                const iExt* RESTRICT iatk, const iReg* RESTRICT jak,                       
                const iReg* RESTRICT J, const iReg sizeJ,                                  
                uint8_t* RESTRICT seen);
