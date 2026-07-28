#pragma once
#include "precision.h"

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void cptRhoJ2(const iReg JtildeSize, 
              double * RESTRICT normColJ, 
              const iExt* RESTRICT jatAJtilde,
              const iReg * const * RESTRICT iaAJtilde, 
              const double * const * RESTRICT coefAJtilde,
              const double * RESTRICT res,
              double * RESTRICT tmpRes,
              const double normRes, const double * RESTRICT colANorm,                       
              const iReg * RESTRICT Jtilde);

iReg minIdx(double *rhoJ2, iReg JtildeSize);

void cptRes(iReg sizeJ, const iExt * RESTRICT jatAJ, 
            const iReg * const * RESTRICT iaAJ, 
            const double * const * RESTRICT coefAJ, 
            const double * RESTRICT mHat,
            const iReg * RESTRICT A0k_idx, const double * RESTRICT A0k, iReg A0k_nnz,
            double * RESTRICT res, iReg * RESTRICT L, iReg &usedL,
            double &resRelNorm, double &resNorm,
            int* RESTRICT ws_idx, double* RESTRICT ws_val);
