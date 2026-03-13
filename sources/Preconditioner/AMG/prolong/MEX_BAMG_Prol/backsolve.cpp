
#include "backsolve.h"

//----------------------------------------------------------------------------------------

// Performs a back substitution with the upper triangular (non unitary) matrix R
// solving Z_in = Z_out*R
void backsolve(const iReg n, const iReg m, const rExt *const *const R, rExt **Z){

   // Reverse loop on columns of Z
   for (iReg i = n-1; i > -1; i--){
      // Loop over the rows of Z
      for (iReg k = m-1; k > -1; k--){
         // Scale by R[i,i]
         Z[k][i] /= R[i][i];
         // Combine column j with column i
         for (iReg j = i-1; j > -1; j--){
            Z[k][j] -= R[i][j]*Z[k][i];
         }
      }
   }

   return;
}
