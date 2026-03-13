#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include "emin_blas.h"     // to use: DGEMV
//////////////////////////////////
#include <iostream>
#include <stdio.h>
using namespace std;
//////////////////////////////////
#include "gather_f.h"

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
int gather_B_Z(const int np, const int nn, const int nn_C, const int ntv,
               const int *fcnode, const int *iat_patt, const int *ja_patt,
               const int *iat_Tpatt, const int *ja_Tpatt,
               const int *iat_A, const int *ja_A, const double *coef_A,
               const double *const *TV, int *&pt_Z, int *&pt_col_Z, double *&mat_Z,
               double *&vec_f, double *coef_P0){
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// g E tau POSSONO ESSERE DEI VETTORI LOCALI DI DIMENSIONE NTV CHE POI VENGONO CANCELLATI
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Init error code
   int ierr = 0;

   // Pointer to the part of Z entries belonging to each thread
   int pt_THpart[np+1];
   // Pointer to the part of Z columns belonging to each thread
   int pt_THcol[np+1];

   // Compute the number of entries in Z and max block size
   int nrmax_blk = 0;
   #pragma omp parallel num_threads(np) reduction(max:nrmax_blk)
   {
      // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstcol, ncolth, lastcol;
      if (mythid <= resto) {
         ncolth = bsize+1;
         firstcol = mythid*ncolth;
         if (mythid == resto) ncolth--;
      } else {
         ncolth = bsize;
         firstcol = mythid*bsize + resto;
      }
      lastcol = firstcol + ncolth;
      int nt = 0;
      int nc = 0;
      int iend = iat_patt[firstcol];
      for (int i = firstcol; i < lastcol; i++){
         int istart = iend;
         iend = iat_patt[i+1];
         int len = iend - istart;
         int lenc = len-ntv;
         if ( lenc > 0 ){
            nc += lenc;
            nt += len*lenc;
         }
         nrmax_blk = max(nrmax_blk,len);
      }
      pt_THpart[mythid+1] = nt;
      pt_THcol[mythid+1] = nc;
   }
   // Reduce pt_THpart and pt_THcol
   pt_THpart[0] = 0;
   for (int ip = 0; ip < np; ip++) pt_THpart[ip+1] += pt_THpart[ip];
   pt_THcol[0] = 0;
   for (int ip = 0; ip < np; ip++) pt_THcol[ip+1] += pt_THcol[ip];

   // Allocate Z and vec_f
   int nrows_Z = iat_patt[nn];
   int ncols_Z = pt_THcol[np];
   int nterm_Z = pt_THpart[np];
   pt_Z = (int*) malloc( (nn+1)*sizeof(int) );
   pt_col_Z = (int*) malloc( (nn+1)*sizeof(int) );
   mat_Z = (double*) malloc( nterm_Z*sizeof(double) );
   vec_f = (double*) malloc( nrows_Z*sizeof(double) );
   if (mat_Z == nullptr || pt_Z == nullptr || pt_col_Z == nullptr ||
       vec_f == nullptr) return ierr = 3;

   // Allocate shared scratches
   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   double *vec_g = (double*) malloc( nrows_Z*sizeof(double) );
   if (c2glo == nullptr || vec_g == nullptr) return ierr = 1;

   #pragma omp parallel num_threads(np)
   {
      // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstcol, ncolth, lastcol;
      if (mythid <= resto) {
         ncolth = bsize+1;
         firstcol = mythid*ncolth;
         if (mythid == resto) ncolth--;
      } else {
         ncolth = bsize;
         firstcol = mythid*bsize + resto;
      }
      lastcol = firstcol + ncolth;

      // Estimate number of rows and first position for this chunk of columns
      int pos_g = iat_patt[firstcol];
      int ind_Z = pt_THpart[mythid];
      int ind_col_Z = pt_THcol[mythid];

      // Set proper position in g
      double *g_scr = &(vec_g[pos_g]);

      // Query work space for DGEQRF and DORGQR
      lapack_int ierr_lapack;
      lapack_int l_nn = static_cast<lapack_int>(max(ntv,nrmax_blk));
      lapack_int l_mm = static_cast<lapack_int>(ntv);
      lapack_int l_kk = static_cast<lapack_int>(ntv);
      double query_work_1;
      double query_work_2;
      double *BB_scr = nullptr;
      double *tau = nullptr;
      double *work = nullptr;

      lapack_int lwork = -1;
      dgeqrf(&l_nn,&l_mm,BB_scr,&l_nn,tau,&query_work_1,&lwork,&ierr_lapack);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      dorgqr(&l_nn,&l_nn,&l_kk,BB_scr,&l_nn,tau,&query_work_2,&lwork,&ierr_lapack);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;
      lwork = static_cast<lapack_int>(max(query_work_1,query_work_2));

      // Allocate private workspace
      BB_scr = (double*) malloc( (nrmax_blk*nrmax_blk)*sizeof(double) );
      tau = (double*) malloc( (ncolth*ntv)*sizeof(double) );
      work = (double*) malloc( lwork*sizeof(double) );
      if (BB_scr == nullptr || tau == nullptr || work == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      if (ierr > 0) goto exit_pragma;

      // Create mapping from coarse node numbering to global (original) numbering
      #pragma omp for
      for (int i = 0; i < nn; i++){
         int k = fcnode[i];
         if (k >= 0) c2glo[k] = i;
      }

      // Loop over the current chunk of columns
      int ind_g, ind_tau;
      ind_g = 0;
      ind_tau = 0;
      int istart_patt, iend_patt;
      iend_patt = iat_patt[firstcol];
      for (int icol = firstcol; icol < lastcol; icol++){

         // Set pt_Z
         pt_Z[icol] = ind_Z;
         pt_col_Z[icol] = ind_col_Z;

         // Check that this is a FINE node
         if (fcnode[icol] < 0){

            istart_patt = iend_patt;
            iend_patt = iat_patt[icol+1];
            int nr_BB_loc = iend_patt-istart_patt;

            // Check that the row is not empty
            if (nr_BB_loc > 0){
               // Copy TV into BB_scr
               int k = 0;
               for (int i = istart_patt; i < iend_patt; i++){
                  int i_F = c2glo[ja_patt[i]];
                  int kk = k;
                  for (int j = 0; j < ntv; j++){
                     BB_scr[kk] = TV[i_F][j];
                     kk += nr_BB_loc;
                  }
                  k++;
               }
               // Copy TV into g
               for (int j = 0; j < ntv; j++) g_scr[ind_g+j] = TV[icol][j];

               /* form of op(A) & op(B) to use in matrix vector multiplication */
               char const *chn = "N", *cht = "T";
               /* scalar values to use in dgemv */
               double const one = 1.0, mone = -1.0;
               lapack_int const oneint = 1;

               // Compute g = g - BB^T*coef_P0
               lapack_int const b_m = static_cast<lapack_int>( nr_BB_loc );
               lapack_int const b_n = static_cast<lapack_int>( ntv );
               dgemv( cht, &b_m, &b_n, &mone, BB_scr, &b_m, &(coef_P0[istart_patt]),
                      &oneint, &one, &(g_scr[ind_g]), &oneint );

               // Perform QR on BB
               l_nn = static_cast<lapack_int>(nr_BB_loc);
               lapack_int l_ll = l_nn;
               lapack_int i_lpk_1, i_lpk_2, i_lpk_3;
               dgeqrf(&l_nn,&l_mm,BB_scr,&l_ll,&(tau[ind_tau]),work,&lwork,&i_lpk_1);

               // Solve transposed triangular system
               char const *chu = "U";
               dtrtrs(chu,cht,chn,&l_mm,&oneint,BB_scr,&l_ll,&(g_scr[ind_g]),&l_mm,&i_lpk_2);

               // Transform QQ from Householder rotation to standard form
               dorgqr(&l_nn,&l_nn,&l_kk,BB_scr,&l_ll,&(tau[ind_tau]),work,&lwork,&i_lpk_3);
               if (i_lpk_1 || i_lpk_2 || i_lpk_3){
                  #pragma omp atomic write
                  ierr = 4;
                  goto exit_loop_icol;
               }

               // Update coef_P0 to ensure the TV constraint: coef_P0 += QQ*g
               dgemv( chn, &b_m, &b_n, &one, BB_scr, &b_m, &(g_scr[ind_g]),
                      &oneint, &one, &(coef_P0[istart_patt]), &oneint );

               // Copy the last columns of QQ (stored in BB_scr) into Z
               int nt_blk_Z = nr_BB_loc*(nr_BB_loc-ntv);
               memcpy(&(mat_Z[ind_Z]),&(BB_scr[nr_BB_loc*ntv]),nt_blk_Z*sizeof(double));

               // Update pointers in Z, g and tau
               ind_Z += nt_blk_Z;
               ind_col_Z += nr_BB_loc-ntv;
               ind_g += ntv;
               ind_tau += ntv;
            }

         }
      } // End loop over columns
      exit_loop_icol: ;

      // Free private workspace
      free(BB_scr);
      free(tau);
      free(work);

      // Exit point
      exit_pragma: ;

   } // End of parallel region

   // Gather vec_f
   #pragma omp parallel for num_threads(np)
   for (int icol = 0; icol < nn_C; icol++){
       int istart = iat_Tpatt[icol];
       int iend = iat_Tpatt[icol+1];
       int irow = c2glo[icol];
       int istart_A = iat_A[irow];
       int iend_A = iat_A[irow+1];
       gather_f(iend-istart,&(ja_Tpatt[istart]),iend_A-istart_A,&(ja_A[istart_A]),
                &(coef_A[istart_A]),&(vec_f[istart]));
   }

   // Set last pointer in pt_Z
   pt_Z[nn] = nterm_Z;
   pt_col_Z[nn] = ncols_Z;

   // Free shared scratches
   free(c2glo);
   free(vec_g);

   return ierr;

}
