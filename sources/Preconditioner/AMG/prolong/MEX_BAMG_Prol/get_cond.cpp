
#include "get_cond.h"

//----------------------------------------------------------------------------------------

// Retrieves rank and conditioning of an input matrix A of size n x m using a rank
// revealing QR based on Householder rotations.
void get_cond(const rExt condmax, const iReg n, const iReg m, rExt **A,
              iReg &rank, rExt &cond){

   iReg maxrank = std::min(n,m);
   rExt wrot[m];

   // Set the initial guess for the rank
   rank = 1;

   // Swap the maximum norm row to the first position
   // Its norm is an estimate of the first SVD
   iReg irow;
   rExt max_sv_est;
   cond = 1.0;
   move_row_front(0,n,m,A,irow,max_sv_est);

   // Compute the corresponding Householder rotation
   mk_HouHolVec(m,A[0],wrot);

   // Apply the rotation to the other rows of A
   for (iReg i=1; i<n; i++) Apply_HouHol_Rot(m,wrot,A[i]);

   while (rank < maxrank){

      rExt norm;

      // Increase the rank representing the row index of row to process
      rank += 1;

      // Swap the maximum norm row to the first position
      move_row_front(rank-1,n-rank+1,m,&(A[rank-1]),irow,norm);

      // Check conditioning
      rExt cond_new = max_sv_est / fabs(norm);
      if (cond_new > condmax){
         rank -= 1;
         break;
      }
      cond = cond_new;

      // Compute the corresponding Householder rotation
      mk_HouHolVec(m-rank+1,&(A[rank-1][rank-1]),wrot);

      // Apply the rotation to the other rows of A
      for (iReg i=rank; i<n; i++) Apply_HouHol_Rot(m-rank+1,wrot,&(A[i][rank-1]));

   }

}
