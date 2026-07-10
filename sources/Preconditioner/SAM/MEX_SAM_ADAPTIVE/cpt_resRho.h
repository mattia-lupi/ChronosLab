#pragma once
#include "precision.h"

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void cptRhoJ2(const iReg JtildeSize, double *normColJ, const iExt* RESTRICT jatAJtilde, 
              const iReg* RESTRICT iaAJtilde, const double* RESTRICT coefAJtilde, 
              double *res, double *tmpRes, const double normRes);

iReg minIdx(double *rhoJ2, iReg JtildeSize);

void cptRes(iReg sizeJ, const iExt * RESTRICT jatAJ, const iReg * RESTRICT iaAJ,           
            const double * RESTRICT coefAJ, const double * RESTRICT mHat,                  
            const iReg * RESTRICT A0k_idx, const double * RESTRICT A0k, iReg A0k_nnz,      
            double * RESTRICT res, iReg * RESTRICT L, iReg &usedL,                         
            double &resRelNorm, double &resNorm,                                           
            int* RESTRICT ws_idx, double* RESTRICT ws_val);
