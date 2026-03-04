#include "icholRF_apply.h"

void icholRF_apply(const iReg nn, const iExt* __restrict__ iU, const iReg* __restrict__ jU,
                   const rExt* __restrict__ coef_U, const rExt* __restrict__ D_inv,
                   const rExt* __restrict__ vec, rExt* __restrict__ pvec){

   // Initialize pvec
   for (iReg i = 0; i < nn; i++) pvec[i] = 0.0;

   // Forward substitution
   iExt iend = iU[0];
   for (iReg k = 0; k < nn; k++){
      iExt istart = iend;
      iend = iU[k+1];
      pvec[k] = vec[k] - pvec[k];
      for (iExt m = istart; m < iend; m++) pvec[jU[m]] += coef_U[m]*pvec[k];
   }

   // Scale by D_inv
   for (iReg i = 0; i < nn; i++) pvec[i] *= D_inv[i];

   // Backward substitution
   iend = iU[nn];
   for (iReg k = nn-1; k >= 0; k--){
      iExt istart = iend-1;
      iend = iU[k];
      rExt a = 0.0;
      for (iExt m = istart; m >= iend; m--) a += coef_U[m]*pvec[jU[m]];
      pvec[k] -= a;
   }

}
