#include <stdlib.h>
#include <omp.h>
#include <iostream>
using namespace std;

#include "ichol_mf.h"
#include "gather_col.h"
#include "get_diagpos.h"
#include "ICHOL_wrapper.h"

/*****************************************************************************************
 *
 * This function performs the matrix-free application of incomplete Cholesky
 *
 * ierr:   error code
 *         0 ---> successful run
 *         1 ---> not enough scratch for ichol
 *         2 ---> error in factorizing matrix
*****************************************************************************************/
int ichol_matfree(const int np, const int nn_C, struct ichol_mf ic_var,
                  const int *iat_A, const int *ja_A, const double *coef_A,
                  const int *iat_Pcol, const int *ja_Pcol, const double *vec_x,
                  double *vec_y){

   // Init error code
   int ierr = 0;

   // Extract some info from ic_var
   int min_lfil  = ic_var.min_lfil;
   int max_lfil  = ic_var.max_lfil;
   int D_lfil    = ic_var.D_lfil;
   int nn_max    = ic_var.max_nrows;
   int max_nnzr  = ic_var.max_nnzr;
   int I_dim     = ic_var.I_dim;
   int R_dim     = ic_var.R_dim;

   #pragma omp parallel num_threads(np)
   {

      // Get thread ID
      int mythid = omp_get_thread_num();
      
      // Get local part of the scratch
      int *myscr_I    = &(ic_var.I_scr[I_dim*mythid]);
      double *myscr_R = &(ic_var.R_scr[R_dim*mythid]);

      // Get max dimension of local arrays
      int nt_max = nn_max*max_nnzr;
      int iwk_U = nt_max + nn_max * (max_lfil+1);
      int ireg_scr_size = iwk_U + nn_max;
      int iext_scr_size = iwk_U + 4*nn_max;

      // Set pointers to int
      int *iat_K = myscr_I;
      int *ja_K = iat_K + nn_max + 1;
      int *iU = ja_K + nt_max;
      int *jU = iU + nn_max + 1;
      int *ireg_scr = jU + iwk_U;
      int *iext_scr = ireg_scr + ireg_scr_size;

      // Set pointers to double
      double *coef_K = myscr_R;
      double *coefU = coef_K + nt_max;
      double *D_inv = coefU + iwk_U;

      // Init first pointer in iat_K (FORTRAN style)
      iat_K[0] = 1;

      #pragma omp for
      for (int icol = 0; icol < nn_C; icol++){

         // Gather the K block (in FORTRAN style)
         int istart = iat_Pcol[icol];
         int jstart = istart;
         int iend   = iat_Pcol[icol+1];
         int nn_K =  iend - istart;
         int irow_K = 0;
         int ind_K = 0;
         for (int j = 0; j < nn_K; j++){
            int n_added;
            int irow = ja_Pcol[jstart];
            // Gather the row indices of this column (upper part only)
            int ind_diag = get_diagpos(irow,iat_A,ja_A);
            int len_A = iat_A[irow+1] - ind_diag;
            // The offset starts from 1 to be consistent with the FORTRAN style
            gather_col(nn_K-j,1+j,&(ja_Pcol[jstart]),len_A,&(ja_A[ind_diag]),
                       &(coef_A[ind_diag]),n_added,&(ja_K[ind_K]),&(coef_K[ind_K]));
            ind_K += n_added;
            irow_K++;
            iat_K[irow_K] = ind_K + 1;
            jstart++;
         }
         int nt_K = ind_K - 1;

         // Factorize the block K with ICHOL
         int lfil = min_lfil;
         int jcol_offset = 0;
         bool FACT_flag = true;
         while (FACT_flag){
            int ierr_ichol =
                ICHOL_wrapper(lfil,jcol_offset,nn_K,nt_K,ireg_scr_size,iext_scr_size,
                              iat_K,ja_K,coef_K,iU,jU,coefU,D_inv,ireg_scr,iext_scr);
            if (ierr_ichol == 0){
               // Successful factorization
               FACT_flag = false;
            } else if (ierr_ichol == 1){
               // Not enough scratch
               cout << "ICHOL ERROR: NOT ENOUGH SCRATCH FOR ICHOL" << endl;
               #pragma omp atomic write
               ierr = 1;
               break;
            } else if (ierr_ichol == 2){
               // Increase fill-il if possible
               if (lfil == max_lfil){
                  // Maximum fill-in reached
                  cout << "ICHOL ERROR: MAXIMUM FILL-IN REACHED FOR BLOK " << icol << endl;
                  #pragma omp atomic write
                  ierr = 2;
                  break;
               }
               cout << "INCREASING FILL-IN FOR BLOCK " << icol << endl;
               lfil += D_lfil;
            }
         }

         // Perform Forward and Backward substitution
         double *pvec = &(vec_y[istart]);
         const double *vec = &(vec_x[istart]);

         // Initialize pvec
         for (int i = 0; i < nn_K; i++) pvec[i] = 0.0;

         // Forward substitution
         int iend_U = iU[0]-1;
         for (int k = 0; k < nn_K; k++){
            int istart_U = iend_U;
            iend_U = iU[k+1]-1;
            pvec[k] = vec[k] - pvec[k];
            for (int m = istart_U; m < iend_U; m++) pvec[jU[m]-1] += coefU[m]*pvec[k];
         }

         // Scale by D_inv
         for (int i = 0; i < nn_K; i++) pvec[i] *= D_inv[i];

         // Backward substitution
         iend_U = iU[nn_K]-1;
         for (int k = nn_K-1; k >= 0; k--){
            int istart_U = iend_U-1;
            iend_U = iU[k]-1;
            double a = 0.0;
            for (int m = istart_U; m >= iend_U; m--) a += coefU[m]*pvec[jU[m]-1];
            pvec[k] -= a;
         }

      } // End loop on blocks

   }

   return ierr;

}
