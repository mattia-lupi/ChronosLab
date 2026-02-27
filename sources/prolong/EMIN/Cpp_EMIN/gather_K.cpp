#include <stdlib.h>
#include <omp.h>
#include <algorithm>

#include "get_diagpos.h"
#include "gather_col.h"

/*****************************************************************************************
 *
 * This function gathers from the prolongation pattern (given by columns) and the operator
 * matrix A, the block diagonal matrix K that is used to minimize the energy of each
 * column of the prolongation.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> allocation error for private scratches
 * 3 ---> allocation error for the final output
 *
*****************************************************************************************/
int gather_K(const int np, const int nn, const int nn_C,
             const int *iat_A, const int *ja_A, const double *coef_A,
             const int *iat_Pcol, const int *ja_Pcol, int &max_nrows_blk,
             int &max_nterm_blk, int &nn_K, int *&iat_K, int *&ja_K,
             double *&coef_K){

   // Init error code
   int ierr = 0;

   // Estimate average non zeros per row in A
   int avg_nnzr_A = iat_A[nn] / nn + 1;

   // Allocate scratches
   int *indrow_scr = (int*) malloc( (np+1)*sizeof(int) );
   int *pt_scr = (int*) malloc( (np+1)*sizeof(int) );
   if (indrow_scr == nullptr || pt_scr == nullptr) return ierr = 1;

   int nt_K;
   max_nrows_blk = 0;
   max_nterm_blk = 0;
   #pragma omp parallel num_threads(np) reduction(max:max_nrows_blk,max_nterm_blk)
   {
      int n_entries,n_added;
      int k_scr, ind_scr;
      int nrows_prev;
      int myrows, myterms;
      int row_offset, nterm_offset;

      // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn_C/np;
      int resto = nn_C%np;
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

      // Estimate max non-zeroes per col
      int istart,iend;
      int max_nnzr = 0;
      iend = iat_Pcol[firstcol];
      for (int icol = firstcol; icol < lastcol; icol++){
         istart = iend;
         iend = iat_Pcol[icol+1];
         max_nnzr = std::max(max_nnzr,iend-istart);
      }
      max_nrows_blk = max_nnzr;

      // Allocate buffer
      int nrows_scr = iat_Pcol[lastcol] - iat_Pcol[firstcol];
      int scr_size = (avg_nnzr_A*nrows_scr + max_nnzr*max_nnzr);
      int *iat_scr = (int*) malloc( (nrows_scr+1)*sizeof(int) );
      int *ja_scr = (int*) malloc( scr_size*sizeof(int) );
      double *coef_scr = (double*) malloc( scr_size*sizeof(double) );
      if (iat_scr == nullptr || ja_scr == nullptr || coef_scr == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;

      // Loop over the chunk of columns of this thread
      k_scr = 0;
      ind_scr = 0;
      nrows_prev = iat_Pcol[firstcol];
      iat_scr[0] = 0;
      iend = iat_Pcol[firstcol];
      for (int icol = firstcol; icol < lastcol; icol++){
         istart = iend;
         iend = iat_Pcol[icol+1];
         n_entries = iend-istart;
         // Check there is enough room in scratch
         if ( ind_scr+n_entries*n_entries >= scr_size ){
            // Evaluate memory needs (+20%)
            double exp_fact = 1.2 * static_cast<double>(lastcol-firstcol) /
                                    static_cast<double>(icol);
            int new_scr_size = static_cast<int>(exp_fact*static_cast<double>(scr_size));
            // Try to reallocate scratch
            ja_scr   = (int*) realloc( ja_scr , new_scr_size*sizeof(int) );
            coef_scr = (double*) realloc( coef_scr, new_scr_size*sizeof(double) );
            if (ja_scr == nullptr || coef_scr == nullptr){
               #pragma omp atomic write
               ierr = 2;
               goto exit_loop_icol;
            }
            scr_size = new_scr_size;
         }
         int offset = k_scr + nrows_prev;
         int ind_scr_old = ind_scr;
         for (int j = 0; j < n_entries; j++){
            int irow = ja_Pcol[istart];
            // Gather the row indices of this column (upper part only)
            int ind_diag = get_diagpos(irow,iat_A,ja_A);
            int len_A = iat_A[irow+1] - ind_diag;
            gather_col(n_entries-j,offset+j,&(ja_Pcol[istart]),len_A,&(ja_A[ind_diag]),
                       &(coef_A[ind_diag]),n_added,&(ja_scr[ind_scr]),&(coef_scr[ind_scr]));
            ind_scr += n_added;
            k_scr++;
            iat_scr[k_scr] = ind_scr;
            istart++;
         }
         max_nterm_blk = std::max(max_nterm_blk,ind_scr-ind_scr_old);
      }
      exit_loop_icol: ;

      // Store the number of rows and entries of this stripe
      myrows = k_scr;
      indrow_scr[mythid+1] = myrows;
      myterms = ind_scr;
      pt_scr[mythid+1] = myterms;
      #pragma omp barrier
      if (ierr > 0) goto exit_pragma;

      #pragma omp single
      {
         // Reduce indrow_scr and pt_scr
         indrow_scr[0] = 0;
         pt_scr[0] = 0;
         for (int ip = 0; ip < np; ip++){
            indrow_scr[ip+1] += indrow_scr[ip];
            pt_scr[ip+1] += pt_scr[ip];
         }
         nn_K = indrow_scr[np];
         nt_K = pt_scr[np];
         // Allocare K
         iat_K  = (int*) malloc( (nn_K+1)*sizeof(int) );
         ja_K   = (int*) malloc( nt_K*sizeof(int) );
         coef_K = (double*) malloc( nt_K*sizeof(double) );
         if (iat_K == nullptr || ja_K == nullptr || coef_K == nullptr){
            ierr = 3;
         }
      }
      if (ierr > 0) goto exit_pragma;

      // Each thread copies its chunk of rows
      row_offset = indrow_scr[mythid];
      nterm_offset = pt_scr[mythid];
      for (int k = 0; k < myrows; k++) iat_K[row_offset+k] = iat_scr[k] + nterm_offset;
      if (mythid == np-1) iat_K[nn_K] = nt_K;
      for (int k = 0; k < myterms; k++){
         ja_K[nterm_offset+k] = ja_scr[k];
         coef_K[nterm_offset+k] = coef_scr[k];
      }

      // Free local buffers
      free(iat_scr);
      free(ja_scr);
      free(coef_scr);
 
      // Exit point
      exit_pragma: ;

   } // End parallel region

   // Free shared scratches
   free(indrow_scr);
   free(pt_scr);

   return ierr;

}
