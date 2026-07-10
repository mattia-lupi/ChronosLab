#include "cpt_resRho.h"
#include "precision.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm> // For std::max and std::min_element
#include <iterator>  // For std::distance
#include <cmath>
#include "lapack.h"

// blas Fortran routines declarations
extern "C" {
   void dgemv_(const char* trans, const lapack_int* m, const lapack_int* n,
               const double* alpha, const double* a, const lapack_int* lda,
               const double* x, const lapack_int* incx, const double* beta,
               double* y, const lapack_int* incy);
               
   double ddot_(const lapack_int* n, const double* x, const lapack_int* incx, 
                const double* y, const lapack_int* incy);

   double dnrm2_(const lapack_int* n, const double* x, const lapack_int* incx);
   
   void daxpy_(const lapack_int* n, const double* alpha, const double* x, 
               const lapack_int* incx, double* y, const lapack_int* incy);
}

void cptRhoJ2(const iReg JtildeSize, double *normColJ, const iExt __restrict *jatAJtilde,  
              const iReg __restrict *iaAJtilde, const double __restrict *coefAJtilde,      
              const iExt nn_A, double *res, double *tmpRes, const double normRes) {
   const double normResSq = normRes * normRes;
   double dot, sum, val, rho;
   iExt colStart, colEnd;
   iReg row;

   // Loop through each column of the CSC matrix
   for (iReg i = 0; i < JtildeSize; ++i) {
      dot = 0.0;
      sum = 0.0;

      // Get the start and end boundaries for the current column 'i'
      colStart = jatAJtilde[i];
      colEnd   = jatAJtilde[i + 1];

      // Iterate only over the non-zero elements of this column
      for (iExt k = colStart; k < colEnd; ++k) {
         row   = iaAJtilde[k];
         val = coefAJtilde[k];

         dot += val * val;
         sum += val * res[row];
      }

      if (dot == 0.0) {
         dot = 1e25;
      }

      tmpRes[i] = sum;

      // Compute the rhoJ2 inline
      rho = normResSq - (sum * sum) / dot;
      normColJ[i] = std::max(rho, 0.0);
   }

   // Copy the result where needed
   std::memcpy(res, tmpRes, JtildeSize * sizeof(double));

   return;
}

// Find the index of the minimum 
iReg minIdx(double *rhoJ2, iReg JtildeSize){
    if (JtildeSize <= 0) return 0;
    
    double* min_element_ptr = std::min_element(rhoJ2, rhoJ2 + JtildeSize);
    return static_cast<iReg>(std::distance(rhoJ2, min_element_ptr));
}

void cptRes(iReg nn_A, iReg sizeJ, const double * __restrict A0k,
            const iExt * __restrict jatAJ, const iReg * __restrict iaAJ,
            const double * __restrict coefAJ, const double * __restrict mHat,
            double * __restrict res, double &resRelNorm, double &resNorm,
            int* __restrict ws_idx, double* __restrict ws_val) {

    // Initialize res and compute baseline sq_sum
    double sq_sum0 = 0.0, sq_sum1 = 0.0, sq_sum2 = 0.0, sq_sum3 = 0.0;
    iReg i = 0;

    // Unrolled from 0 to nn_A - 3, compute residual in a sparse accumulator
    for (; i <= nn_A - 4; i += 4) {
        double a0 = A0k[i];   double a1 = A0k[i+1];
        double a2 = A0k[i+2]; double a3 = A0k[i+3];

        res[i]   = -a0; res[i+1] = -a1;
        res[i+2] = -a2; res[i+3] = -a3;

        sq_sum0 += a0 * a0; sq_sum1 += a1 * a1;
        sq_sum2 += a2 * a2; sq_sum3 += a3 * a3;
    }
    // Loop from nn_A - 3 to the end
    for (; i < nn_A; ++i) {
        double a = A0k[i];
        res[i] = -a;
        sq_sum0 += a * a;
    }
    double sq_sum = sq_sum0 + sq_sum1 + sq_sum2 + sq_sum3;

    // Accumulate into dense ws_val, log unique rows in ws_idx, compute SpMV
    int ws_count = 0;
    for (iReg j = 0; j < sizeJ; ++j) {
        const double m_j = mHat[j];
        if (m_j == 0.0) continue;

        iExt start = jatAJ[j];
        iExt end   = jatAJ[j + 1];
        for (iExt k = start; k < end; ++k) {
            iReg row = iaAJ[k];
            double val = coefAJ[k] * m_j;

            // If this row hasn't been touched yet in this call, log its index
            if (ws_val[row] == 0.0) {
                ws_idx[ws_count++] = row;
            }
            ws_val[row] += val;
        }
    }

    // Compute normAjMh and dynamically adjust res and sq_sum
    double normAjMh_sq = 0.0;
    for (int n = 0; n < ws_count; ++n) {
        iReg row = ws_idx[n];
        double val = ws_val[row];

        // Handles both true zeros and skips duplicate indices safely
        if (val == 0.0) continue;

        normAjMh_sq += val * val;

        double old_res = res[row]; // Currently holds -A0k[row]
        double new_res = old_res + val;
        res[row] = new_res;

        // Compute the total sum of squares
        sq_sum += (new_res * new_res) - (old_res * old_res);

        // Reset workspace to 0.0 for the next function call
        ws_val[row] = 0.0;
    }

    // Compute the relative residual norm
    double normAjMh = std::sqrt(normAjMh_sq);
    resNorm = std::sqrt(sq_sum);

    resRelNorm = 2.0 * resNorm / (normAjMh + resRelNorm);
}
