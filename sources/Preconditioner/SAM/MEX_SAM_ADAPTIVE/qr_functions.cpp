#include "qr_functions.h"
#include <iostream>

void computeFirstQR(double *Ahat, lapack_int sizeI, lapack_int sizeJ, double *R, 
                    double *Rtriang, double *tau, double *work, lapack_int lwork, 
                    lapack_int &info){
   // Compute QR Factorization
   dgeqrf_(&sizeI, &sizeJ, Ahat, &sizeI, tau, work, &lwork, &info);

   if (info != 0){
      printf("Exit at first QR due to error %d\n", static_cast<int>(info));
      return;
   }

   // Copy R into a separate vector to store and use also the updated R later
   R[0] = Ahat[0];
   Rtriang[0] = Ahat[0];

   return;
}

void applyFirstQt(double *Ahat, lapack_int sizeI, lapack_int sizeJ, double *tau, 
                  double *a0k, double *work, lapack_int lwork, lapack_int &info){
   char side = 'L';
   char trans = 'T';

   #if defined(__clang__) && !defined(MATLAB_MEX_FILE)
      dormqr_(&side, &trans, &sizeI, &sizeJ, &sizeJ,
              Ahat, &sizeI, tau, a0k, &sizeI, work, &lwork, &info,1,1);
   #else
      dormqr_(&side, &trans, &sizeI, &sizeJ, &sizeJ,
              Ahat, &sizeI, tau, a0k, &sizeI, work, &lwork, &info);
   #endif
   if (info != 0){
      printf("Exit at first Qt apply due to error %d\n", static_cast<int>(info));
      return;
   }

   return;
}

void applyR(lapack_int sizeJ, double *R, double *a0k, lapack_int &info){
   char uplo = 'U';
   char trans = 'N';
   char diag = 'N';
   lapack_int nrhs = 1;

   #if defined(__clang__) && !defined(MATLAB_MEX_FILE)
      dtrtrs_(&uplo, &trans, &diag, &sizeJ, &nrhs, 
              R, &sizeJ, a0k, &sizeJ, &info,1,1,1); 
   #else
      dtrtrs_(&uplo, &trans, &diag, &sizeJ, &nrhs, 
              R, &sizeJ, a0k, &sizeJ, &info); 
   #endif

   if (info != 0) {
       printf("Exit at R apply due to error %d\n", static_cast<int>(info));
       return;
   }

   return;
}

void applyQt(iReg t, const lapack_int* RESTRICT sizeJ, 
             const lapack_int* RESTRICT sizeI, lapack_int *qStart, 
             double *Ahat, double *tau, double *a0k, lapack_int nRowsRHS, 
             lapack_int ncolsRHS, double *work, lapack_int lwork, lapack_int &info){

   // Apply iteratively the Qs in the correct way
   char side = 'L';
   char trans = 'T';

   // Loop over all the matrices Q I have computed and apply them to the correct point
   lapack_int nrows, ncols, nrefl, LDA, LDC, ofA0k = 0, ofTau = 0; 
   for (iReg i = 0; i < t; ++i){
      // Get the values for the new iteration matrix
      ofA0k = sizeJ[i];
      nrows = sizeI[i+1] - ofA0k;
      ncols = ncolsRHS;
      nrefl = sizeJ[i+1] - ofA0k;
      LDA   = nrows;
      LDC   = nRowsRHS;
      ofTau = ofA0k;

      #if defined(__clang__) && !defined(MATLAB_MEX_FILE)
         dormqr_(&side, &trans, &nrows, &ncols, &nrefl, Ahat + qStart[i] + ofA0k, 
                 &LDA, tau + ofTau, a0k + ofA0k, &LDC, work, &lwork, &info,1,1);
      #else
         dormqr_(&side, &trans, &nrows, &ncols, &nrefl, Ahat + qStart[i] + ofA0k, 
                 &LDA, tau + ofTau, a0k + ofA0k, &LDC, work, &lwork, &info);
      #endif

      if (info != 0){
         printf("Exit at Qt apply due to error %d\n", static_cast<int>(info));
         printf("nrows(3) %d, ncols(4) %d, nrefl(5) %d, lda(7) %d, ldc(10) %d\n", 
                static_cast<int>(nrows), static_cast<int>(ncols), static_cast<int>(nrefl), 
                static_cast<int>(LDA), static_cast<int>(LDC));
         return;
      }
   }
   return;
}


// colSizeB2 = J_add size
// rowSizeB2 = J_add + I_add size
// startB2   = oldSizeI*n2old + n2-n2old
// rowSizeB  = sizeI
void computeNewQR(iReg t, lapack_int *sizeI, lapack_int *sizeJ, lapack_int *qStart, 
                  double *Ahat, double *tau, double *R, double *Rtriang, double *work, 
                  lapack_int lwork, lapack_int &info){
   // Compute QR Factorization on the new part of the matrix
   lapack_int oldSizeTau = sizeJ[t];
   lapack_int colSizeB2  = sizeJ[t+1] - oldSizeTau;
   lapack_int rowSizeB   = sizeI[t+1];
   lapack_int rowSizeB2  = sizeI[t+1] - oldSizeTau;
   lapack_int sqrtStartR = oldSizeTau;
   // printf("sizeIt %d, sizeJt %d\n", sizeI[t+1],oldSizeTau);
   // Set the starting point for the new QR factorization
   qStart[t+1]    = qStart[t] + rowSizeB*colSizeB2;
   lapack_int startB2    = qStart[t] + oldSizeTau;

   // Get the filled size of Rtriang
   lapack_int startRtri = 0.5*sqrtStartR*(sqrtStartR+1);

   // Copy the data into the triangular R first (colmajor)
   lapack_int startB = qStart[t];
   // printf("startB2 %d startB %d at t %d\n", startB2,startB, t);
   
   // Compute the new QR factorization
   dgeqrf_(&rowSizeB2, &colSizeB2, Ahat + startB2, &rowSizeB,
          tau + oldSizeTau, work, &lwork, &info);

   if (info != 0){
      printf("Exit at first QR due to error %d\n", static_cast<int>(info));
      return;
   }

   // Loop over the columns
   for (iReg j = 0; j < colSizeB2; ++j) {
      // Loop over rows, stopping EXACTLY at the diagonal of the current column
      for (iReg i = 0; i <= oldSizeTau + j; ++i) { 
         // This now grabs ONLY valid R entries (both old and new)
         Rtriang[startRtri] = Ahat[startB + i + j * rowSizeB];
         // Update the starting point for packed storage
         startRtri++;
      }
   }

   // Set counters
   iReg newSizeR = sqrtStartR + colSizeB2;
   iReg count = 0;
   
   // Loop over the columns of the new matrix
   for (iReg i = 0; i < newSizeR; ++i){
      for (iReg j = 0; j <= i; ++j){
         R[i*newSizeR + j] = Rtriang[count];
         count++;
      }
   }
   
   return;
}
