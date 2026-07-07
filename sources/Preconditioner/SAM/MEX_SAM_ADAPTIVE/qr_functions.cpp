#include "lapack.h"
#include "precision.h"
#include "qr_functions.h"
#include <iostream>


// void print_vector1(const char* desc, ptrdiff_t n, double* vec) {
//     printf("\n--- %s ---\n", desc);
//     for (ptrdiff_t i = 0; i < n; i++) {
//         printf("%10.4g\n", vec[i]);
//     }
// }

void computeFirstQR(double *Ahat, ptrdiff_t sizeI, ptrdiff_t sizeJ, double *R, double *Rtriang, double *tau, double *work, ptrdiff_t lwork, ptrdiff_t &info){
   // Compute QR Factorization
   dgeqrf_(&sizeI, &sizeJ, Ahat, &sizeI, tau, work, &lwork, &info);
   if (info != 0){
      printf("Exit at first QR due to error %td\n", info);
      return;
   }

   // Copy R into a separate vector to store and use also the updated R later
   R[0] = Ahat[0];
   Rtriang[0] = Ahat[0];

   return;
}

void applyFirstQt(double *Ahat, ptrdiff_t sizeI, ptrdiff_t sizeJ, double *tau, double *a0k, double *work, ptrdiff_t lwork, ptrdiff_t &info){
   char side = 'L';
   char trans = 'T';

   // Check why I needed to add the fortran string length sizes
   dormqr_(&side, &trans, &sizeI, &sizeJ, &sizeJ, Ahat, &sizeI, tau, a0k, &sizeI, work, &lwork, &info);//,3,3
   if (info != 0){
      printf("Exit at first Qt apply due to error %td\n", info);
      return;
   }

   return;
}

void applyR(ptrdiff_t sizeJ, double *R, double *a0k, ptrdiff_t &info){
   char uplo = 'U';
   char trans = 'N';
   char diag = 'N';
   ptrdiff_t nrhs = 1;

   // Check why I needed to add the fortran string length sizes
   dtrtrs_(&uplo, &trans, &diag, &sizeJ, &nrhs, R, &sizeJ, a0k, &sizeJ, &info);// ,3,3,3
   if (info != 0){
      printf("Exit at R apply due to error %td\n", info);
      return;
   }

   return;
}

void applyQt(ptrdiff_t t, ptrdiff_t *sizeJ, ptrdiff_t *sizeI, ptrdiff_t *qStart, double *Ahat, 
             double *tau, double *a0k, ptrdiff_t nRowsRHS, ptrdiff_t ncolsRHS, double *work, ptrdiff_t lwork, ptrdiff_t &info){

   // Apply iteratively the Qs in the correct way
   char side = 'L';
   char trans = 'T';

   // Loop over all the matrices Q I have computed and apply them to the correct point
   ptrdiff_t nrows, ncols, nrefl, LDA, LDC, ofA0k = 0, ofTau = 0; 
   for (ptrdiff_t i = 0; i < t; ++i){
      // Get the values for the new iteration matrix
      nrows = sizeI[i+1] - sizeJ[i];
      ncols = ncolsRHS;
      nrefl = sizeJ[i+1] - sizeJ[i];
      LDA   = nrows;
      LDC   = nRowsRHS;
      ofA0k = sizeJ[i];
      ofTau = sizeJ[i];
//       printf("nrows %d, ncols %d, nrefl %d, LDA %d, LDC %d, ofA0k %d, ofTau %d\n",nrows,ncols,nrefl,LDA,LDC,ofA0k,ofTau);
//       printf("it %d, ahat(1) %f, tau(1) %f, rhs(1) %f\n",i,Ahat[qStart[i] + sizeJ[i]],tau[ofTau],a0k[ofA0k]);

      // Check why I needed to add the fortran string length sizes
      dormqr_(&side, &trans, &nrows, &ncols, &nrefl, Ahat + qStart[i] + sizeJ[i], 
              &LDA, tau + ofTau, a0k + ofA0k, &LDC, work, &lwork, &info);// ,3,3
      // print_vector1("chat vec iter",nrowsA0k,a0k);// check why +i above
      if (info != 0){
         printf("Exit at Qt apply due to error %td\n", info);
         printf("nrows(3) %td, ncols(4) %td, nrefl(5) %td, lda(7) %td, ldc(10) %td\n", nrows, ncols, nrefl, LDA, LDC);
         return;
      }
   }
   return;
}


// colSizeB2 = J_add size
// rowSizeB2 = J_add + I_add size
// startB2   = oldSizeI*n2old + n2-n2old
// rowSizeB  = sizeI
void computeNewQR(ptrdiff_t t, ptrdiff_t *sizeI, ptrdiff_t *sizeJ, ptrdiff_t *qStart, double *Ahat, double *tau, double *R, 
                  double *Rtriang, double *work, ptrdiff_t lwork, ptrdiff_t &info){
   // Compute QR Factorization on the new part of the matrix
   ptrdiff_t oldSizeTau = sizeJ[t];
   ptrdiff_t colSizeB2  = sizeJ[t+1] - oldSizeTau;
   ptrdiff_t rowSizeB   = sizeI[t+1];
   ptrdiff_t rowSizeB2  = sizeI[t+1] - oldSizeTau;
   ptrdiff_t sqrtStartR = oldSizeTau;
   // printf("sizeIt %d, sizeJt %d\n", sizeI[t+1],oldSizeTau);
   // Set the starting point for the new QR factorization
   qStart[t+1]    = qStart[t] + rowSizeB*colSizeB2;
   ptrdiff_t startB2    = qStart[t] + oldSizeTau;

   // Get the filled size of Rtriang
   ptrdiff_t startRtri = 0.5*sqrtStartR*(sqrtStartR+1);

   // Copy the data into the triangular R first (colmajor)
   ptrdiff_t startB = qStart[t];
   // printf("startB2 %d startB %d at t %d\n", startB2,startB, t);
   
   // Compute the new QR factorization
   dgeqrf_(&rowSizeB2, &colSizeB2, Ahat+startB2, &rowSizeB, tau+oldSizeTau, work, &lwork, &info);
   if (info != 0){
      printf("Exit at first QR due to error %td\n", info);
      return;
   }
   
   // Loop over the columns
   for (ptrdiff_t j = 0; j < colSizeB2; ++j){
      // Loop over rows
      for (ptrdiff_t i = 0; i < rowSizeB; ++i){
         // Save the new values carefully adjusting for the leading dimension
         Rtriang[startRtri] = Ahat[startB + i + j*rowSizeB];
         // Update the starting point
         startRtri++;
      }
   }

   // Set counters
   ptrdiff_t newSizeR = sqrtStartR + colSizeB2;
   ptrdiff_t count = 0;
   
   // Loop over the columns of the new matrix
   for (ptrdiff_t i = 0; i < newSizeR; ++i){
      for (ptrdiff_t j = 0; j <= i; ++j){
         R[i*newSizeR + j] = Rtriang[count];
         count++;
      }
   }
   
   return;
}
