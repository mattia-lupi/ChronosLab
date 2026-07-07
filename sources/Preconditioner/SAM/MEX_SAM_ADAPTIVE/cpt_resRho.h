#pragma once

#include "precision.h"

void cptRes(iExt nn_A, iReg sizeJ, double *A0k, double *AJ, double *mHat, double *res, double &resRelNorm, double &resNorm);

void cptRhoJ2(iReg JtildeSize, double *normColJ, double *AJtilde, iExt nn_A, double *res, double *tmpRes, double normRes);

iReg minIdx(double *rhoJ2, iReg JtildeSize);
