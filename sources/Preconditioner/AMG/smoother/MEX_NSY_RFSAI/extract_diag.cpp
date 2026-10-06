#include <iostream>
#include "extract_diag.h"
#include <cstdlib>
#include <omp.h>

void extract_diag(const int nrows, const int *const iat, const int *const ja,
                  const double *const coef, double *&diag, const int num_threads){

   // Allocate output
   diag = (double*) malloc(nrows * sizeof(double));
   if (diag == nullptr){
      // Allocation error
      std::cout << "Allocation Error in extract_diag" << std::endl;
      return;
   }

   // Loop over matrix rows
   #pragma omp parallel for schedule(static) num_threads(num_threads)
   for (int i = 0; i < nrows; i++){
      double d = 0.0;
      int istart = iat[i];
      int iend = iat[i+1];
      for (int j = istart; j < iend; j++){
         if (ja[j] == i){
            d = coef[j];
            break;
         }
      }
      diag[i] = d;
   }
}
