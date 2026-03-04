#include "omp.h"

void icholRF_apply(const int np, const int nblk, const int* __restrict__ pt_blk,
                   const int* __restrict__ iU, const int* __restrict__ jU,
                   const double* __restrict__ coef_U, const double* __restrict__ D_inv,
                   const double* __restrict__ vec, double* __restrict__ pvec){

   #pragma omp parallel for num_threads(np)
   for (int iblk = 0; iblk < nblk; iblk++){
      int firstrow = pt_blk[iblk];
      int lastrow = pt_blk[iblk+1];

      // Initialize pvec
      for (int i = firstrow; i < lastrow; i++) pvec[i] = 0.0;

      // Forward substitution
      int iend = iU[firstrow];
      for (int k = firstrow; k < lastrow; k++){
         int istart = iend;
         iend = iU[k+1];
         pvec[k] = vec[k] - pvec[k];
         for (int m = istart; m < iend; m++) pvec[jU[m]] += coef_U[m]*pvec[k];
      }

      // Scale by D_inv
      for (int i = firstrow; i < lastrow; i++) pvec[i] *= D_inv[i];

      // Backward substitution
      iend = iU[lastrow];
      for (int k = lastrow-1; k >= firstrow; k--){
         int istart = iend-1;
         iend = iU[k];
         double a = 0.0;
         for (int m = istart; m >= iend; m--) a += coef_U[m]*pvec[jU[m]];
         pvec[k] -= a;
      }

   }

}
