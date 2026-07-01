#include "lapack.h"
#include "cblas.h"
#include "qr_functions.h"



void applyFirstQt(double *Ahat, int sizeI, int sizeJ, double *tau, double *a0k, double *work, int lwork, int &info){
   char side = 'L';
   char trans = 'T';

   // Check why I needed to add the fortran string length sizes
   dormqr_(&side, &trans, &sizeI, &sizeJ, &sizeJ, Ahat, &sizeI, tau, a0k, &sizeI, work, &lwork, &info,3,3);
   if (info != 0){
      printf("Exit at first Qt apply due to error %d\n", info);
      return;
   }

   return;
}

void applyR(int sizeJ, double *R, double *a0k, int &info){
   char uplo = 'U';
   char trans = 'N';
   char diag = 'N';
   int nrhs = 1;

   // Check why I needed to add the fortran string length sizes
   dtrtrs_(&uplo, &trans, &diag, &sizeJ, &nrhs, R, &sizeJ, a0k, &sizeJ, &info,3,3,3);
   if (info != 0){
      printf("Exit at R apply due to error %d\n", info);
      return;
   }

   return;
}

void applyOldQ(int t, int oldSizeJ, int sizeJ, int oldSizeI, int sizeI, 
               double *Ahat, double *tau, double *a0k, double *work, int lwork, int &info){

   if(t == 1){
      // If it's the first step then it's easy to reapply the Q'
      applyFirstQt(Ahat, oldSizeI, sizeJ-oldSizeJ, tau, Ahat + oldSizeJ*oldSizeI, work, lwork, info);
   }
   else{
      // Apply iteratively the Qs in the correct way

   }

   return;
}

void computeFirstQR(double *Ahat, int sizeI, int sizeJ, double *R, double *Rtriang, double *tau, double *work, int lwork, int &info){
   // Compute QR Factorization
   dgeqrf_(&sizeI, &sizeJ, Ahat, &sizeI, tau, work, &lwork, &info);
   if (info != 0){
      printf("Exit at first QR due to error %d\n", info);
      return;
   }

   // Copy R into a separate vector to store and use also the updated R later
   R[0] = Ahat[0];
   Rtriang[0] = Ahat[0];

   return;
}

// colSizeB2 = J_add size
// rowSizeB2 = J_add + I_add size
// startB2   = oldSizeI*n2old + n2-n2old
// rowSizeB  = sizeI
void computeNewQR(int rowSizeB2, int colSizeB2, int startB2, int rowSizeB, int oldSizeTau, double *Ahat,
                  double *tau, double *R, int sqrtStartR, double *Rtriang, double *work, int lwork, int &info){
   // Compute QR Factorization on the new part of the matrix
   dgeqrf_(&rowSizeB2, &colSizeB2, Ahat+startB2, &rowSizeB, tau+oldSizeTau, work, &lwork, &info);
   if (info != 0){
      printf("Exit at first QR due to error %d\n", info);
      return;
   }
   // Get the filled size of Rtriang
   int startRtri = 0.5*sqrtStartR*(sqrtStartR+1);

   // Copy the data into the triangular R first (colmajor)
   int sizeAddR = colSizeB2*rowSizeB;
   int startB = startB2 - rowSizeB + rowSizeB2;
   
   // Loop over the columns
   for (int j = 0; j < colSizeB2; ++j){
      // Loop over rows
      for (int i = 0; i < rowSizeB2; ++i){
         // Save the new values carefully adjusting for the leading dimension
         Rtriang[startRtri] = Ahat[startB + i + j*rowSizeB];
         // Update the starting point
         startRtri++;
      }
   }

   // Set counters
   int newSizeR = sqrtStartR + colSizeB2;
   int count = 0;
   
   // Loop over the columns of the new matrix
   for (int i = 0; i < newSizeR; ++i){
      for (int j = 0; j <= i; ++j){
         R[i*newSizeR + j] = Rtriang[count];
         count++;
      }
   }
   
   return;
}
