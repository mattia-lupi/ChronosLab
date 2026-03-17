#include <iostream>
#include "extract_diag.h"

void extract_diag(const int nrows, const int *const iat, const int *const ja,
                  const double *const coef, double *&diag){

   // Allocate output
   diag = (double*) malloc(nrows * sizeof(double));
   if (diag == nullptr){
      // Allocation error
      std::cout << "Allocation Error in extract_diag" << std::endl;
      return;
   }
   
   // Loop over matrix rows
   int iend = iat[0];
   for ( int i = 0; i < nrows; i++ ){
      int istart = iend;
      iend = iat[i+1];
      int jcol = ja[istart];
      if (jcol == i) {
         diag[i] = coef[istart];
         continue;
      }
      diag[i] = 0.0;
      while (istart < iend-1){
         jcol = ja[++istart];
         if (jcol == i) {
            diag[i] = coef[istart];
            break;
         }
      }
   }

}
