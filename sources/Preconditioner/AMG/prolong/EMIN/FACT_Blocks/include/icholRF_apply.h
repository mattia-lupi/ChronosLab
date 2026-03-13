#pragma once
#include "precision.h"

void icholRF_apply(const iReg nn, const iExt* __restrict__ iU, const iReg* __restrict__ jU,
                   const rExt* __restrict__ coef_U, const rExt* __restrict__ D_inv,
                   const rExt* __restrict__ vec, rExt* __restrict__ pvec);
