#include "transpose.h"

int transpose(const int nrows, const int ncols, const int *const iat, const int *const ja,
              const double *const coef, int *&iat_T, int *&ja_T, double *&coef_T){

   // Allocate output and scracth
   iat_T = (int*) malloc((ncols+1) * sizeof(int));
   int nterm = iat[nrows];
   ja_T = (int*) malloc((nterm) * sizeof(int));
   coef_T = (double*) malloc((nterm) * sizeof(double));
   int *ISCR = (int*) malloc((ncols+1) * sizeof(int));
   if (iat_T == nullptr || ja_T == nullptr || coef_T == nullptr || ISCR == nullptr){
      // Allocation error
      cout << "Allocation Error in transpose" << endl;
      return 1;
   }
   
   // Initialize pointers
   fill_n(iat_T,ncols+1,0);

   // Count non-zeroes for each column of the input matrix
   for ( int i = 0; i < nrows; i++ ){
      for ( int j = iat[i]; j < iat[i+1]; j++ ) iat_T[ja[j]]++;
   }

   // Set pointers
   ISCR[0] = 0;
   for ( int i = 1; i < ncols+1; i++ ) ISCR[i] = ISCR[i-1] + iat_T[i-1];
   for ( int i = 0; i < ncols+1; i++ ) iat_T[i] = ISCR[i];

   // Transpose column indices and coefficients
   for ( int i = 0; i < nrows; i++ ){
      for ( int j = iat[i]; j < iat[i+1]; j++ ){
         int ind  = ISCR[ja[j]];
         ja_T[ind] = i;
         coef_T[ind] = coef[j];
         ISCR[ja[j]] = ind+1;
      }
   }

   // Deallocate scratch
   free(ISCR);

   return 0;

}
