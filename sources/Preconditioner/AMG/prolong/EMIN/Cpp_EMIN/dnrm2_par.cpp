#include <cmath>
#include "omp.h"
#include "dnrm2_par.h"

double dnrm2_par(const int np, const int nn, const double* RESTRICT x,
                 double* RESTRICT reduc){

   (void)reduc;
   double dnrm2 = 0.0;
   #pragma omp parallel for num_threads(np) reduction(+:dnrm2)
   for ( int i = 0; i < nn; i++ ){
      dnrm2 += x[i] * x[i];
   }

   return std::sqrt(dnrm2);

}
