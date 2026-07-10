#pragma once

#include "precision.h"

void cptRhoJ2(const iReg JtildeSize, double *normColJ, const iExt __restrict *jatAJtilde, 
              const iReg __restrict *iaAJtilde, const double __restrict *coefAJtilde, 
              double *res, double *tmpRes, const double normRes);

iReg minIdx(double *rhoJ2, iReg JtildeSize);

void cptRes(iReg nn_A, iReg sizeJ, const double * __restrict A0k,
            const iExt * __restrict jatAJ, const iReg * __restrict iaAJ,
            const double * __restrict coefAJ, const double * __restrict mHat,
            double * __restrict res, double &resRelNorm, double &resNorm,
            int* __restrict ws_idx, double* __restrict ws_val);
