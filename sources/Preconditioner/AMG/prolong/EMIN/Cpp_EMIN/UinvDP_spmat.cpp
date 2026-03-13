#include <omp.h>
#include "UinvDP_spmat.h"

// Second part of symmetric Gauss-Seidel: (D^-1*U)^-1 * P
// Computes local MxM. C = P^T*triu(D^-1*A)^-1 with the pattern of C already defined by
// the pattern of B.
// NOTE: works in place!
// It assumes A to be symmetric (otherwise a transposition is needed before this call)

void UinvDP_spmat(const iReg nthreads, const iReg nrows_A, const iExt* iat_A,
                  const iReg *ja_A, rExt *coef_A, const iReg ncols_B, const iExt *iat_B,
                  const iReg *ja_B, const rExt *coef_B, iReg* WNALL){

   // Open parallel region
   #pragma omp parallel num_threads(nthreads)
   {

      // Get thread ID
      int mythid = omp_get_thread_num();

      // Each thread gets its own part of WNALL
      iReg* WN = WNALL + mythid*ncols_B;

      // Initialize WN
      for (iReg i = 0; i < ncols_B; i++){
         WN[i] = -1;
      }

      // Loop over A rows assigned to the thread
      #pragma omp for
      for (iReg irow = 0; irow < nrows_A; irow++){

         // Poiter to the first entry of current row
         iExt ind_C = iat_A[irow];

         // Zero-out entries
         iReg nn_row = 0;
         iExt Cstop = iat_A[irow+1];
         for (iExt j = ind_C; j < Cstop; j++){
            iReg jcol = ja_A[j];
            WN[jcol] = nn_row;
            nn_row++;
         }

         // Explore the i-th row of A
         iExt ISTRT = iat_A[irow];
         iExt ISTOP = iat_A[irow+1];
         for (iExt j = ISTOP-1; j >= ISTRT; j--){

            iReg jjA = ja_A[j];
            rExt coef = coef_A[j];

            // Explore the jjA-th row of B
            iExt JSTRT = iat_B[jjA];
            iExt JSTOP = iat_B[jjA+1];

            // Get diagonal coefficient
            rExt diag;
            for (iExt k = JSTOP-1; k >= JSTRT; k--){
               if( ja_B[k] == jjA )
               {
                 diag = coef_B[k];
               }
            }
            for (iExt k = JSTOP-1; k >= JSTRT; k--){

               // Column index of B
               iReg pos = WN[ja_B[k]];

               if (pos >= 0) {
                  // Column in the list, sum this contribution
                  if( ja_B[k] > jjA )
                  {
                    coef -= coef_A[ind_C+pos] * coef_B[k] / diag;
                  }
                  else if( ja_B[k] == jjA )
                  {
                    coef_A[ind_C+pos] = coef;
                  }
               }
            }
         }

         // Sparse resetting of the list
         for ( iReg j = 0; j < nn_row; j ++ ) {
            WN[ja_A[ind_C+j]] = -1;
         }

      }

   }

}
