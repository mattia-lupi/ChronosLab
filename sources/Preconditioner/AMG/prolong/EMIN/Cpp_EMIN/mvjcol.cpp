#include "mvjcol.h"
// Transpose pattern indices.
void mvjcol( const int firstrow, const int nrows, const int* RESTRICT iat,
             const int* RESTRICT ja, int* RESTRICT ja_T, int* RESTRICT punt,
             int* __restrict perm){

   // Transpose local stripe
   int shift = firstrow;
   for ( int i = 0; i < nrows; i++ ) {
      for ( int j = iat[i]; j < iat[i+1]; j++ ) {
         int irow = ja[j];
         perm[j] = punt[irow];
         ja_T[punt[irow]] = i + shift;
         punt[irow]++;
      }
   }

}
