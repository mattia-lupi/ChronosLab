#include "spmatv_blk.h"

void spmatv_blk(const int np, const int nblk, const int* RESTRICT pt_blk,
                const int* RESTRICT iat, const int* RESTRICT ja,
                const double* RESTRICT coef, const double* RESTRICT v_in,
                double* RESTRICT v_out){

   #pragma omp parallel for num_threads(np) 
   for (int i = 0; i < pt_blk[nblk]; i++) v_out[i] = 0.0;

   #pragma omp parallel for num_threads(np)
   for (int iblk = 0; iblk < nblk; iblk++){
      int istart = pt_blk[iblk];
      int iend = pt_blk[iblk+1];
      int jend = iat[istart];
      for (int i = istart; i < iend; i++){
         int jstart = jend;
         jend = iat[i+1];
         v_out[i] += coef[jstart]*v_in[ja[jstart]];
         for(int j = jstart+1; j < jend; j++){
            v_out[i]     += coef[j]*v_in[ja[j]];
            v_out[ja[j]] += coef[j]*v_in[i];
         }
      }
   }

}
