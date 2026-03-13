#include "omp.h"
#include "emin_blas.h"

// Performs the orthogonalization v_out = v_in - Q*QT*v_in
void Orth_Q(const int np, const int nn, const int nn_K, const int ntv,
            const int* __restrict__ iat_patt, const int* __restrict__ ja_patt,
            const double* __restrict__ mat_Q, const double* __restrict__ v_in,
            double* __restrict__ v_ntv, double* __restrict__ v_out){

   /* form of op(A) & op(B) to use in matrix vector multiplication */
   char const *chn = "N", *cht = "T";
   /* scalar values to use in dgemv */
   double const one = 1.0, zero = 0.0;
   lapack_int const oneint = 1;

   #pragma omp parallel num_threads(np)
   {
      // Select the local scratch for v_ntv
      int myid = omp_get_thread_num();
      double *my_v_ntv = &(v_ntv[myid*ntv]);
      #pragma omp for
      for (int icol = 0; icol < nn; icol++){
         int istart = iat_patt[icol];
         int iend = iat_patt[icol+1];
         int nr_loc = iend - istart;
         if (nr_loc > 0){
            int ind_Q = istart*ntv;

            // Perform product by QT: my_v_ntv <-- QT*v_in;
            lapack_int const b_m = static_cast<lapack_int>( nr_loc );
            lapack_int const b_n = static_cast<lapack_int>( ntv );
            dgemv( cht, &b_m, &b_n, &one, &(mat_Q[ind_Q]), &b_m,
                   &(v_in[istart]), &oneint, &zero, my_v_ntv, &oneint );

            // Perform product by Q: v_out <-- Q*my_v_ntv;
            dgemv( chn, &b_m, &b_n, &one, &(mat_Q[ind_Q]), &b_m,
                   my_v_ntv, &oneint, &zero, &(v_out[istart]), &oneint );

            // Compute v_out <-- v_in - v_out
            for (int j = istart; j < iend; j++) v_out[j] = v_in[j] - v_out[j];
         }
      }
   }

}
