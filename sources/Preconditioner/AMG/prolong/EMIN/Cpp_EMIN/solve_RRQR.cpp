#include "blas.h"     // to use: DGEMV
#include "lapacke.h"  // to use: DGEQRF and DORGQR

/*****************************************************************************************
 *
 * This function gathers from the prolongation pattern (given by rows) and the test vector
 * array TV, the block diagonal matrix B that is used to enforce the constraint. B is
 * immediately factorized with QR and only Q is returned.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> allocation error for private scratches
 * 3 ---> allocation error for the final output
 * 4 ---> lapack error
 *
*****************************************************************************************/

int solve_RRQR( int const istart_patt, int const iend_patt, int const * const ja_patt,
                int const * const c2glo, int const ntv, double const * const * const TV,
                int const nr_BB_loc, int const l_mm, int const l_nn, double * const BBT,
                lapack_int * const JPVT, int const ind_tau, double * const tau, int const lwork,
                double * const work, double const condmax, int const ind_g,
                double * const g_scr, double * const Rmat, double * const RRT,
                double * const delta, int const ind_BB, double * const BB_scr,
                double * const coef_P0 )
{
   int ierr = 0;

   //-------------------------------------------------------------+
   // Local BB is rank deficient use rank-revealing factorization |
   //-------------------------------------------------------------+

   int i_lpk_1, i_lpk_2, i_lpk_3, i_lpk_4, i_lpk_5;

   // Copy TV into BBT (BB_scr transposed)
   int k = 0;
   for (int i = istart_patt; i < iend_patt; i++){
      int i_F = c2glo[ja_patt[i]];
      for (int j = 0; j < ntv; j++){
         BBT[k] = TV[i_F][j];
         k++;
      }
   }

   // Perform QR with column pivoting on BBT
   for (int i = 0; i < nr_BB_loc; i++) JPVT[i] = 0;
   i_lpk_2 = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR,l_mm,l_nn,BBT,l_mm,JPVT,
                                 &(tau[ind_tau]),work,lwork);

   // Find rank of BBT
   int rank = 0;
   while (fabs(BBT[0]/BBT[rank*ntv+rank]) < condmax){
      rank++;
      if (rank == ntv) break;
   }

   // Multiply g by QT: g = QT*g
   lapack_int l_rank = static_cast<lapack_int>(rank);
   i_lpk_3 = LAPACKE_dormqr_work(LAPACK_COL_MAJOR,'L','T',l_mm,1,l_mm,BBT,
                                 l_mm,&(tau[ind_tau]),&(g_scr[ind_g]),
                                 l_mm,work,lwork);

   /* form of op(A) & op(B) to use in matrix vector multiplication */
   char const *chn = "N", *cht = "T", *chl = "L", *chu = "U";
   /* scalar values to use in dgemv */
   double const one = 1.0, mone = -1.0, zero = 0.0;
   lapack_int const oneint = 1;

   // Compute R*RT
   for (int j = 0; j < nr_BB_loc; j++){
      int col = ntv*j;
      for (int i = 0; i < rank; i++) Rmat[col+i] = (j>=i) ? BBT[col+i]:0.0;
   }
   lapack_int b_m = static_cast<lapack_int>( rank );
   lapack_int b_n = static_cast<lapack_int>( rank );
   lapack_int b_k = static_cast<lapack_int>( nr_BB_loc );
   lapack_int lda = static_cast<lapack_int>( l_mm );
   lapack_int ldb = static_cast<lapack_int>( l_mm );
   lapack_int ldc = static_cast<lapack_int>( l_mm );
   dgemm(chn,cht,&b_m,&b_n,&b_k,&one,Rmat,&lda,Rmat,&ldb,&zero,RRT,&ldc);

   // Compute g = inv(mat_RRT)*g;
   i_lpk_4 = LAPACKE_dpotrf_work(LAPACK_COL_MAJOR,'U',l_rank,RRT,l_mm);
   i_lpk_5 = LAPACKE_dpotrs_work(LAPACK_COL_MAJOR,'U',l_rank,1,RRT,l_mm,
                                 &(g_scr[ind_g]),l_mm);
   if (i_lpk_1 || i_lpk_2 || i_lpk_3 || i_lpk_4 || i_lpk_5){
      #pragma omp atomic write
      ierr = 4;
      return ierr;
   }

   // Compute delta = mat_R^T * res
   b_m = static_cast<lapack_int>( rank );
   b_n = static_cast<lapack_int>( nr_BB_loc );
   lda = static_cast<lapack_int>( ntv );
   dgemv(cht,&b_m,&b_n,&one,Rmat,&lda,&(g_scr[ind_g]),&oneint,&zero,
         delta,&oneint);

   // Update p0
   for (int i = 0; i < nr_BB_loc; i++){
      coef_P0[istart_patt+JPVT[i]-1] +=  delta[i];
   }

   // Create the orthogonal projector Q (temporarily storing it permuted
   // in Rmat
   b_m = static_cast<lapack_int>( rank );
   b_n = static_cast<lapack_int>( nr_BB_loc );
   lda = static_cast<lapack_int>( ntv );
   ldb = static_cast<lapack_int>( ntv );
   dtrsm(chl,chu,cht,chn,&b_m,&b_n,&one,RRT,&lda,Rmat,&ldb);

   // Permute columns of Q with PVT and store it in BB_scr
   for (int j = 0; j < nr_BB_loc; j++){
      int col = j*ntv;
      int pcol = ind_BB + (JPVT[j]-1)*ntv;
      for (int i = 0; i < rank; i++) BB_scr[pcol+i] = Rmat[col+i];
      for (int i = rank; i < ntv; i++) BB_scr[pcol+i] = 0.0;
   }

   return ierr;
}
