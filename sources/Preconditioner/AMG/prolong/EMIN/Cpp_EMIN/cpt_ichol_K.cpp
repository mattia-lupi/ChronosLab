#include <stdlib.h>
#include <omp.h>
#include <iostream>
using namespace std;

#include "DebEnv.h"
#include "parm_EMIN.h"
#include "ICHOL_wrapper.h"

/*****************************************************************************************
 *
 * This function computes the incomplete factorization of the block diagonal matrix K.
 *
 * Parameters:
 *
 * np:            number of openMP threads.
 * itmax:         number of EnerMinCG iterations.
 * prec_type:     preconditioner for energy minimization
 * min_lfil:      minimum fill-in for Cholesky factorization
 * max_lfil:      maximum fill-in for Cholesky factorization
 * D_lfil:        fill-in increase for Cholesky factorization
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> allocation error for private scratches
 * 3 ---> allocation error for the final output
 * 4 ---> not enough memopry in incomplete Cholesky factorization
 * 5 ---> maximum fill-in reached in incomplete Cholesky factorization
 *
*****************************************************************************************/
int cpt_ichol_K(const int np, const int min_lfil, const int max_lfil, const int D_lfil,
                const int n_blk, const int *pt_blk, const int nnmax_blk,
                const int ntmax_blk, const int nn_K, const int *iat_K, const int *ja_K,
                const double *coef_K, double &avg_lfil, int *&it_U, int *&jcol_U,
                double *&coef_U, double *&D_inv){

   // Init error code
   int ierr = 0;

   // Allocate D_inv
   D_inv = (double*) malloc( nn_K*sizeof(double) );
   if (D_inv == nullptr) return ierr = 3;

    // Allocate scratches
   int *pt_scr = (int*) malloc( (np+1)*sizeof(int) );
   if (pt_scr == nullptr) return ierr = 1;

   // Compute incomplete factorization of each diagonal block
   int nt_U;
   int lfil_total = 0;
   #pragma omp parallel num_threads(np) reduction(+:lfil_total)
   {
      // Get thread ID and row partition
      int mythid = omp_get_thread_num();
      int bsize = n_blk/np;
      int resto = n_blk%np;
      int firstblk, nblkth, lastblk;
      if (mythid <= resto) {
         nblkth = bsize+1;
         firstblk = mythid*nblkth;
         if (mythid == resto) nblkth--;
      } else {
         nblkth = bsize;
         firstblk = mythid*bsize + resto;
      }
      lastblk = firstblk + nblkth;

      // Allocate buffer to store a stripe of factorization
      int firstrow_K = pt_blk[firstblk];
      int lastrow_K = pt_blk[lastblk];
      int nrows_scr = lastrow_K - firstrow_K;
      int nterm_scr = (iat_K[lastrow_K]-iat_K[firstrow_K]) +
                      nrows_scr*(min_lfil+max_lfil)/2;
      int *iat_scr = (int*) malloc( (nrows_scr+1)*sizeof(int) );
      int *ja_scr = (int*) malloc( nterm_scr*sizeof(int) );
      double *coef_scr = (double*) malloc( nterm_scr*sizeof(double) );
      if (iat_scr == nullptr || ja_scr == nullptr || coef_scr == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;

      // Allocate scratch to locally store a block
      int *iat_loc, *ja_loc;
      double *coef_loc;
      iat_loc = (int*) malloc( (nnmax_blk+1)*sizeof(int) ); 
      ja_loc = (int*) malloc( ntmax_blk*sizeof(int) );
      coef_loc = (double*) malloc( ntmax_blk*sizeof(double) );
      if (iat_loc == nullptr || ja_loc == nullptr || coef_loc == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;
      
      // Allocate scratch for ICHOL
      int iwk_U, ireg_scr_size, iext_scr_size;
      int *ireg_scr, *iext_scr;
      iwk_U = ntmax_blk + nnmax_blk * (max_lfil+1);
      ireg_scr_size = nnmax_blk + iwk_U;
      iext_scr_size = 4*nnmax_blk + iwk_U;
      ireg_scr = (int*) malloc( ireg_scr_size*sizeof(int) );
      iext_scr = (int*) malloc( iext_scr_size*sizeof(int) );
      if (ireg_scr == nullptr || iext_scr == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;

      // Loop over the blocks of this thread
      int ind_scr;
      int *iU_loc, *jU_loc;
      double *coefU_loc;
      ind_scr = 0;
      iU_loc = iat_scr;
      jU_loc = ja_scr;
      coefU_loc = coef_scr;
      for (int iblk = firstblk; iblk < lastblk; iblk++){

         // Copy the block into the local area (switching indices to FORTRAN style)
         int ind_loc = 0;
         iat_loc[0] = 1;
         int istart_loc = pt_blk[iblk];
         int iend_loc = pt_blk[iblk+1];
         int nn_loc = iend_loc - istart_loc;
         for (int i = istart_loc; i < iend_loc; i++){
            for (int j = iat_K[i]; j < iat_K[i+1]; j++){
               ja_loc[ind_loc] = ja_K[j] - istart_loc + 1;
               coef_loc[ind_loc] = coef_K[j];
               ind_loc++;
            }
            iat_loc[i-istart_loc+1] = ind_loc + 1;
         }
         int nt_loc = ind_loc;

         // Factorize the block with ICHOL
         int lfil = min_lfil;
         int jcol_offset = 0;
         bool FACT_flag = true;
         int ind_scr_old = ind_scr;
         while (FACT_flag){
            int ierr_ichol =
                ICHOL_wrapper(lfil,jcol_offset,nn_loc,nt_loc,ireg_scr_size,iext_scr_size,
                              iat_loc,ja_loc,coef_loc,iU_loc,jU_loc,coefU_loc,
                              &(D_inv[istart_loc]),ireg_scr,iext_scr);
            if (ierr_ichol == 0){
               // Successful factorization
               FACT_flag = false;
            } else if (ierr_ichol == 1){
               // Not enough scratch
               cout << "ICHOL ERROR: NOT ENOUGH SCRATCH FOR ICHOL" << endl;
               #pragma omp atomic write
               ierr = 4;
               goto exit_loop_icol;
            } else if (ierr_ichol == 2){
               // Increase fill-il if possible
               if (lfil == max_lfil){
                  // Maximum fill-in reached
                  cout << "ICHOL ERROR: MAXIMUM FILL-IN REACHED FOR BLOK " <<
                          firstblk+iblk << endl;
                  #pragma omp atomic write
                  ierr = 5;
                  goto exit_loop_icol;
               }
               cout << "INCREASING FILL-IN FOR BLOCK " << iblk << endl;
               lfil += D_lfil;
            }
         }
         int nt_U_loc = iU_loc[nn_loc] - 1;
         ind_scr += nt_U_loc;

         // Check size of ja_scr/coef_scr
         if (ind_scr > nterm_scr){
            /*****************************************************************************
            // DA SISTEMARE MEGLIO @@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // Evaluate memory needs (+20%)
            double exp_fact = 1.2 * static_cast<double>(lastblk-firstblk) /
                                    static_cast<double>(iblk);
            int new_nterm_scr = static_cast<int>(exp_fact*static_cast<double>(nterm_scr));
            // Try to reallocate scratch
            ja_scr   = (int*) realloc( ja_scr , new_nterm_scr*sizeof(int) );
            coef_scr = (double*) realloc( coef_scr, new_nterm_scr*sizeof(double) );
            if (ja_scr == nullptr || coef_scr == nullptr){
               #pragma omp atomic write
               ierr = 2;
               goto exit_loop_icol;
            }
            // Reset pointers
            jU_loc = ja_scr;
            coefU_loc = coef_scr;
            nterm_scr =new_nterm_scr;
            *******************************************************************************/ 
            // Not enough scratch
            cout << "ICHOL ERROR: NOT ENOUGH SCRATCH FOR ICHOL" << endl;
            #pragma omp atomic write
            ierr = 4;
            goto exit_loop_icol;
         }
         // Adjust indices (from Fortran as well)
         for (int i = 0; i < nn_loc+1; i++) iU_loc[i] += ind_scr_old - 1;
         for (int j = 0; j < nt_U_loc; j++) jU_loc[j] += istart_loc - 1;

         // Update local pointers
         iU_loc += nn_loc;
         jU_loc += nt_U_loc;
         coefU_loc += nt_U_loc;

         // Increase fill-in count
         lfil_total += lfil;

      } // End loop on blocks
      exit_loop_icol: ;

      // Free ICHOL scratches
      free(ireg_scr);
      free(iext_scr);
      free(iat_loc);
      free(ja_loc);
      free(coef_loc);

      // Store the number of entries of this stripe
      int myterms;
      myterms = ind_scr;
      pt_scr[mythid+1] = myterms;
      #pragma omp barrier

      #pragma omp single
      {  
         // Reduce pt_scr
         pt_scr[0] = 0;
         for (int ip = 0; ip < np; ip++) pt_scr[ip+1] += pt_scr[ip];
         nt_U = pt_scr[np];
         // Allocare final U
         it_U  = (int*) malloc( (nn_K+1)*sizeof(int) );
         jcol_U   = (int*) malloc( nt_U*sizeof(int) );
         coef_U = (double*) malloc( nt_U*sizeof(double) );
         if (it_U == nullptr || jcol_U == nullptr || coef_U == nullptr){
            ierr = 3;
         }
      }
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;

      // Each thread copies its chunk of rows
      int row_offset, nterm_offset;
      row_offset = firstrow_K;
      nterm_offset = pt_scr[mythid];
      for (int k = 0; k < nrows_scr; k++) it_U[row_offset+k] = iat_scr[k] + nterm_offset;
      if (mythid == np-1) it_U[nn_K] = nt_U;
      for (int k = 0; k < myterms; k++){
         jcol_U[nterm_offset+k] = ja_scr[k];
         coef_U[nterm_offset+k] = coef_scr[k];
      }

      // local buffers
      free(iat_scr);
      free(ja_scr);
      free(coef_scr);

      // Exit point
      exit_pragma: ;

   } // End parallel region

   // Compute average fill-in
   avg_lfil = static_cast<double>(lfil_total) / static_cast<double>(n_blk);

   // Free shared scratch
   free(pt_scr);

   return ierr;

}
