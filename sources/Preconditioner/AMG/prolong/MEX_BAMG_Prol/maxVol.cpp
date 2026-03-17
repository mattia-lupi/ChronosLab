
#include "maxVol.h"

//----------------------------------------------------------------------------------------

// Applies the max volume algorithm to retrieve the best possible basis among the
// rows of the input matrix A[n,m].
void maxVol(const iReg mmax, const rExt condmax, const iReg itmax, const rExt delta,
            const iReg n, const iReg m, rExt **A, iReg &rank, iReg *list){

   iReg maxrank = std::min(n,m);
   maxrank = std::min(maxrank,mmax);
   rExt wrot[m];

   // Set the initial guess for the rank
   rank = 1;

   //-------------------------------------------------------------------------------------
   // Stage 1 - Find the rank of the matrix and a set of independent rows using a
   // rank reavealing QR
   //-------------------------------------------------------------------------------------

   // Swap the maximum norm row to the first position
   // Its norm is an estimate of the first SVD
   iReg irow;
   rExt max_sv_est;
   move_row_front(0,n,m,A,irow,max_sv_est);
   SWAP(list[0],list[irow]);

   // Compute the corresponding Householder rotation
   mk_HouHolVec(m,A[0],wrot);

   // Apply the rotation to the other rows of A
   for (iReg i=0; i<n; i++) Apply_HouHol_Rot(m,wrot,A[i]);

   while (rank < maxrank){

      rExt norm;

      // Increase the rank representing the row index of row to process
      rank += 1;

      // Swap the maximum norm row to the first position
      move_row_front(rank-1,n-rank+1,m,&(A[rank-1]),irow,norm);
      SWAP(list[rank-1],list[irow+rank-1]);

      // Check conditioning
      rExt cond_new = max_sv_est / abs(norm);
      if (cond_new > condmax){
         rank -= 1;
         break;
      }

      // Compute the corresponding Householder rotation
      mk_HouHolVec(m-rank+1,&(A[rank-1][rank-1]),wrot);

      // Apply the rotation to the other rows of A
      for (iReg i=rank-1; i<n; i++) Apply_HouHol_Rot(m-rank+1,wrot,&(A[i][rank-1]));

   }

   //-------------------------------------------------------------------------------------
   // Stage 2 - Permute rows of A in order to find the best basis
   //-------------------------------------------------------------------------------------

   // The matrix A is now partitioned as [R,B]' with R upper triangular of size rank
   // Now compute Z = [inv(R)*B]' = B'*inv(R')
   backsolve(rank,n-rank,&(A[0]),&(A[rank]));

   // Permute row in order to find the best basis
   try {
      maxVol_inner(itmax,delta,n-rank,rank,list,&(A[rank]));
   } catch (linsol_error) {
      throw linsol_error ("maxVol","permuting row in order to find the best basis");
   }

}
