#include <cmath>
#include "omp.h"

double dnrm2_par(const int np, const int nn, const double* __restrict__ x,
                 double* __restrict__ reduc){

   #pragma omp parallel num_threads(np)
   {
      double dnrm2 = 0.0;
      #pragma omp for
      for ( int i = 0; i < nn; i++ ) dnrm2 += x[i] * x[i];
      int myid = omp_get_thread_num();
      reduc[myid] = dnrm2;
   }

   double dnrm2 = 0.0;
   for ( int ip = 0; ip < np; ip++ ) dnrm2 += reduc[ip];

   return sqrt(dnrm2);

}
