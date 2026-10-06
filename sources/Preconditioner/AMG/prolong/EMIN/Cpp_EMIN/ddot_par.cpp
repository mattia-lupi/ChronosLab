#include "omp.h"
#include "ddot_par.h"

double ddot_par(const int np, const int nn, const double* RESTRICT x,
                const double* RESTRICT y, double* RESTRICT reduc){

   #pragma omp parallel num_threads(np)
   {
      double ddot = 0.0;
      #pragma omp for
      for ( int i = 0; i < nn; i++ ){
         ddot += x[i] * y[i];
      }
      int myid = omp_get_thread_num();
      reduc[myid] = ddot;
   }

   double ddot = 0.0;
   for ( int ip = 0; ip < np; ip++ ) ddot += reduc[ip];

   return ddot;

}
