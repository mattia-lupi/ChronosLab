#include "omp.h"
#include "blas.h"
#include "lapacke.h"

void Z_mult(const int np, const int nblk, const int ntv, const int *pt_blk,
            const int *pt_Z, const int *pt_col_Z, const double *mat_Z,
            const double *vec_in, double* vec_out){

   /* form of op(A) & op(B) to use in matrix vector multiplication */
   char const *chn = "N", *cht = "T";
   /* scalar values to use in dgemv */
   double const one = 1.0, mone = -1.0, zero = 0.0;
   lapack_int const oneint = 1;

    #pragma omp parallel for num_threads(np)
    for (int iblk = 0; iblk < nblk; iblk++){
       int istart_out = pt_blk[iblk];
       int iend_out = pt_blk[iblk+1];
       int nr = iend_out - istart_out;
       if (nr > 0){
          int nc = nr-ntv;
          if (nc > 0){
             int ind_Z = pt_Z[iblk];
             int istart_in = pt_col_Z[iblk];
             lapack_int const b_m = static_cast<lapack_int>( nr );
             lapack_int const b_n = static_cast<lapack_int>( nc );
             dgemv(chn,&b_m,&b_n,&one,&(mat_Z[ind_Z]),&b_m,&(vec_in[istart_in]),
                   &oneint,&zero,&(vec_out[istart_out]),&oneint);
          } else {
             for (int k = istart_out; k < istart_out+nr; k++) vec_out[k] = 0.0;
          }
       }
    }

}
