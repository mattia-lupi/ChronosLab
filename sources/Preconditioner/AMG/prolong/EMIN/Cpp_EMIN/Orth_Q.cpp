#include "omp.h"
#include "emin_blas.h"
#include "Orth_Q.h"

// Performs the orthogonalization v_out = v_in - Q*QT*v_in
void Orth_Q(const int np, const int nn, const int ntv, const int* RESTRICT iat_patt,
            const double* RESTRICT mat_Q, const double* RESTRICT v_in,
            double* RESTRICT v_ntv, double* RESTRICT v_out){

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
         int const istart = iat_patt[icol];
         int const iend = iat_patt[icol+1];
         int const nr_loc = iend - istart;
         if (nr_loc > 0){
            const double* RESTRICT q_blk = &(mat_Q[istart*ntv]);
            const double* RESTRICT x     = &(v_in[istart]);
            double* RESTRICT y           = &(v_out[istart]);

            if (ntv == 1) {
               double a0 = 0.0;
               for (int j = 0; j < nr_loc; j++) a0 += q_blk[j] * x[j];
               for (int j = 0; j < nr_loc; j++) y[j] = x[j] - a0 * q_blk[j];
            } else if (ntv == 6 && nr_loc <= 64) {
               const double* RESTRICT q0 = q_blk;
               const double* RESTRICT q1 = q0 + nr_loc;
               const double* RESTRICT q2 = q1 + nr_loc;
               const double* RESTRICT q3 = q2 + nr_loc;
               const double* RESTRICT q4 = q3 + nr_loc;
               const double* RESTRICT q5 = q4 + nr_loc;
               double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0;
               for (int j = 0; j < nr_loc; j++) {
                  double const xj = x[j];
                  a0 += q0[j] * xj;
                  a1 += q1[j] * xj;
                  a2 += q2[j] * xj;
                  a3 += q3[j] * xj;
                  a4 += q4[j] * xj;
                  a5 += q5[j] * xj;
               }
               for (int j = 0; j < nr_loc; j++) {
                  y[j] = x[j] - (a0 * q0[j] + a1 * q1[j] + a2 * q2[j]
                               + a3 * q3[j] + a4 * q4[j] + a5 * q5[j]);
               }
            } else {
               // Perform product by QT: my_v_ntv <-- QT*v_in;
               lapack_int const b_m = static_cast<lapack_int>( nr_loc );
               lapack_int const b_n = static_cast<lapack_int>( ntv );
               dgemv( cht, &b_m, &b_n, &one, q_blk, &b_m,
                      x, &oneint, &zero, my_v_ntv, &oneint );

               // Perform product by Q: v_out <-- Q*my_v_ntv;
               dgemv( chn, &b_m, &b_n, &one, q_blk, &b_m,
                      my_v_ntv, &oneint, &zero, y, &oneint );

               // Compute v_out <-- v_in - v_out
               for (int j = 0; j < nr_loc; j++) y[j] = x[j] - y[j];
            }
         }
      }
   }

}
