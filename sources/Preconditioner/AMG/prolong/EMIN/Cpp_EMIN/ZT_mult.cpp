#include "omp.h"
#include "emin_blas.h"

void ZT_mult(const int np, const int nblk, const int ntv, const int *pt_blk,
             const int *pt_Z, const int *pt_col_Z, const double *mat_Z,
             const double *vec_in, double* vec_out){

   /* form of op(A) & op(B) to use in matrix vector multiplication */
   char const *cht = "T";
   /* scalar values to use in dgemv */
   double const one = 1.0, zero = 0.0;
   lapack_int const oneint = 1;

    #pragma omp parallel for num_threads(np)
    for (int iblk = 0; iblk < nblk; iblk++){
       int istart_in = pt_blk[iblk];
       int iend_in = pt_blk[iblk+1];
       int nc = iend_in - istart_in;
       if (nc > 0){
          int nr = nc-ntv;
          int istart_out = pt_col_Z[iblk];
          if (nr > 0){
             int ind_Z = pt_Z[iblk];
             lapack_int const b_m = static_cast<lapack_int>( nc );
             lapack_int const b_n = static_cast<lapack_int>( nr );
             dgemv(cht,&b_m,&b_n,&one,&(mat_Z[ind_Z]),&b_m,&(vec_in[istart_in]),
                   &oneint,&zero,&(vec_out[istart_out]),&oneint);
          }
       }
    }

}
