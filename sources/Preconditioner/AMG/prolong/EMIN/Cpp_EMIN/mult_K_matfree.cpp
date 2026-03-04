#include <stdlib.h>
#include <omp.h>

#include "mult_K_col.h"

/*****************************************************************************************
 *
 * This function performs the matrix-free multiplication by K.
 *
*****************************************************************************************/
void mult_K_matfree(const int np, const int nn_C, const int *iat_A, const int *ja_A,
                    const double *coef_A, const int *iat_Pcol, const int *ja_Pcol,
                    const double *vec_x, double *vec_y){

   #pragma omp parallel for num_threads(np)
   for (int icol = 0; icol < nn_C; icol++){

      int istart = iat_Pcol[icol];;
      int iend   = iat_Pcol[icol+1];;
      int n_entries =  iend - istart;

      for (int j = istart; j < iend; j++){
         int irow_A = ja_Pcol[j];
         int ind_A = iat_A[irow_A];
         int len_A = iat_A[irow_A+1] - ind_A;
         vec_y[j] = mult_K_col(n_entries,&(ja_Pcol[istart]),len_A,&(ja_A[ind_A]),
                               &(coef_A[ind_A]),&(vec_x[istart]));
      }

   }

}
