#include "omp.h"
#include "ddot_par.h"

double ddot_par(const int np, const int nn, const double* RESTRICT x,
                const double* RESTRICT y, double* RESTRICT reduc){

   (void)reduc;
   double ddot = 0.0;
   #pragma omp parallel for num_threads(np) reduction(+:ddot)
   for ( int i = 0; i < nn; i++ ){
      ddot += x[i] * y[i];
   }

   return ddot;

}
