#include "parm_EMIN.h"
#include <omp.h>

void load_Jacobi(const int np, const int sol_type, const int nn, const int nt,
                 const int *iat_patt, const int *ja_Tpatt,
                 const int *iat, const int *ja, const double *coef, double *scr,
                 double *D_inv){

   // Load inverse diagonal entries in scr
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn; i++){
      int ind = iat[i];
      while (ja[ind] < i) ind++;
      scr[i] = 1.0 / coef[ind];
   }

   // Load diagonal entries into D_inv:
   // If sol_type == SPMAT, D_inv is arranged in row-major order conforming to iat_patt.
   // Entry (i, k) corresponds to fine row i, so its diagonal is scr[i].
   if (sol_type == SPMAT) {
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn; i++){
         const int istart = iat_patt[i];
         const int iend = iat_patt[i+1];
         const double di = scr[i];
         for (int j = istart; j < iend; j++){
            D_inv[j] = di;
         }
      }
   } else {
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nt; i++){
         D_inv[i] = scr[ja_Tpatt[i]];
      }
   }

}
