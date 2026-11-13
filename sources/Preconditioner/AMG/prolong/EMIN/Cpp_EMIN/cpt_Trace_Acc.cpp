#include "omp.h"
#include "get_diagpos.h"

double cpt_Trace_Acc(const int np, const int nn, const int *fcnode, const int *iat,
                     const int *ja, const double *coef){

   double tr = 0.0;
   #pragma omp parallel for num_threads(np) reduction(+:tr)
   for (int irow = 0; irow < nn; irow++){
      if (fcnode[irow] >= 0){
         int ind = get_diagpos(irow, iat, ja);
         if (ind >= 0) tr += coef[ind];
      }
   }

   return tr;

}
