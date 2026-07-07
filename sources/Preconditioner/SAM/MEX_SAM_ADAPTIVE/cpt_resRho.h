#pragma once

#include <iostream>
#include <cstring>
#include "precision.h"

void cptRes(ptrdiff_t nn_A, ptrdiff_t sizeJ, double *A0k, double *AJ, double *mHat, double *res, double &resRelNorm, double &resNorm);

void cptRhoJ2(ptrdiff_t JtildeSize, double *normColJ, double *AJtilde, ptrdiff_t nn_A, double *res, double *tmpRes, double normRes);

ptrdiff_t minIdx(double *rhoJ2, ptrdiff_t JtildeSize);
