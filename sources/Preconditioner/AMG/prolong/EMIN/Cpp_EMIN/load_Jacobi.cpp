#include "get_diagpos.h"
#include <cmath>

void load_Jacobi(const int np, const int nn, const int *iat_patt,
                 const int *iat, const int *ja, const double *coef, double *scr,
                 double *D_inv){

   // Load diagonal entries in scr and broadcast into row-major D_inv
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn; i++){
      int ind = get_diagpos(i, iat, ja);
      double const inv_d = (ind >= 0 && std::abs(coef[ind]) > 1e-15) ? (1.0 / coef[ind]) : 1.0;
      scr[i] = inv_d;
      int const istart = iat_patt[i];
      int const iend   = iat_patt[i + 1];
      for (int k = istart; k < iend; k++) {
         D_inv[k] = inv_d;
      }
   }

}
