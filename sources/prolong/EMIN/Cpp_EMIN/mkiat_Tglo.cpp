
#include "mkiat_Tglo.h"

//----------------------------------------------------------------------------------------

// Computes the global iat_T and update WI pointers.
void mkiat_Tglo (const int myid, const int nrows, const int nequ, const int nthreads,
                 const int firstrow, int** __restrict__ WI,
                 const int* __restrict__  nnz, int* __restrict__ iat_T  ){

   // Compute the number of non-zeros belonging to previous threads
   int ntprec = 0;
   for ( int i = 0; i < myid; i++ ){
      ntprec += nnz[i];
   }

   // Update iat_T
   for ( int i = 0; i < nrows; i++ ){
      iat_T[i] += ntprec;
   }

   // Account for (nrows+1)-th component for last threads
   if ( myid + 1 == nthreads ) {
      iat_T[nrows] = ntprec + nnz[myid];
   }

   // Update WI1
   for ( int i = firstrow; i < firstrow+nrows; i++ ){
      WI[0][i] = iat_T[i-firstrow];
   }
   for ( int j = 1; j < nthreads; j++ ){
      for ( int i = firstrow; i < firstrow+nrows; i++ ){
         WI[j][i] +=  WI[j-1][i];
      }
   }

}
