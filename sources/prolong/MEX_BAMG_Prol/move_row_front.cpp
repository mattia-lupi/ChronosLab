#include "move_row_front.h"

//----------------------------------------------------------------------------------------

// Moves to the first position the j-th row of A whose norm ( A(j,icol:end) ) is maximal
void move_row_front(const iReg jcol, const iReg n, const iReg m, rExt **A,
                    iReg &irow, rExt &maxnorm){

   // Find the row with max norm
   irow = 0;
   maxnorm = inl_dnrm2(m-jcol,&(A[0][jcol]),1);
   for (iReg i = 1; i < n; i++){
      rExt norm = inl_dnrm2(m-jcol,&(A[i][jcol]),1);
      if (norm > maxnorm){
         maxnorm = norm;
         irow = i;
      }
   }
   if (irow != 0){
      // Swap columns
      for (iReg j = 0; j < m; j++){
         rExt tmp = A[irow][j];
         A[irow][j] = A[0][j];
         A[0][j] = tmp;
      }
   }

   return;
}
