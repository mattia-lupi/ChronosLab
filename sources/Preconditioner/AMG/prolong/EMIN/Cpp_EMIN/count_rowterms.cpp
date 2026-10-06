
#include "count_rowterms.h"

//----------------------------------------------------------------------------------------

// Counts the number of non-zeroes assigned to a thread in each row
void count_rowterms( const int nequ, const int nterm, const int* RESTRICT ja,
                     int* RESTRICT WI ){

   // Initialize WI
   for ( int i = 0; i < nequ; i++ ) {
      WI[i] = 0;
   }

   // Count non-zeroes
   for ( int i = 0; i < nterm; i++ ) {
      WI[ja[i]]++;
   }

}
