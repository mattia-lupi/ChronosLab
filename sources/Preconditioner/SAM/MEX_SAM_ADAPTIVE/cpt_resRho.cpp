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

void cptRhoJ2(const iReg JtildeSize, double *normColJ, const iExt* RESTRICT jatAJtilde,
              const iReg* RESTRICT iaAJtilde, const double* RESTRICT coefAJtilde,
              double *res, double *tmpRes, const double normRes) {
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

      tmpRes[i] = sum;

      // Compute the rhoJ2 inline
      normColJ[i] = normResSq - (sum * sum) / dot;
   }

   return;
}

// Find the index of the minimum 
iReg minIdx(double *rhoJ2, iReg JtildeSize){
    if (JtildeSize <= 0) return 0;
    
    double* min_element_ptr = std::min_element(rhoJ2, rhoJ2 + JtildeSize);
    return static_cast<iReg>(std::distance(rhoJ2, min_element_ptr));
}

void cptRes(iReg sizeJ, const iExt * RESTRICT jatAJ, const iReg * RESTRICT iaAJ, 
            const double * RESTRICT coefAJ, const double * RESTRICT mHat,
            const iReg * RESTRICT A0k_idx, const double * RESTRICT A0k, iReg A0k_nnz,
            double * RESTRICT res, iReg * RESTRICT L, iReg &usedL,
            double &resRelNorm, double &resNorm,
            int* RESTRICT ws_idx, double* RESTRICT ws_val) {

    // Clear the elements of 'res' that were modified 
    // in the previous call, using the old L array.
    for (iReg n = 0; n < usedL; ++n) {
        res[L[n]] = 0.0;
    }

    int ws_count = 0;

    // Compute AJ * mHat and save into the workspace
    for (iReg j = 0; j < sizeJ; ++j) {
        const double m_j = mHat[j];
        if (m_j == 0.0) continue;

        iExt start = jatAJ[j];
        iExt end   = jatAJ[j + 1];
        for (iExt k = start; k < end; ++k) {
            iReg row = iaAJ[k];
            if (ws_val[row] == 0.0) {
                ws_idx[ws_count++] = row;
            }
            ws_val[row] += coefAJ[k] * m_j;
        }
    }

    // Compute the norm looking only at the full entries
    double normAjMh_sq = 0.0;
    for (int n = 0; n < ws_count; ++n) {
        double spmv_val = ws_val[ws_idx[n]];
        normAjMh_sq += spmv_val * spmv_val;
    }

    // Subtract A0k from the workspace
    for (iReg n = 0; n < A0k_nnz; ++n) {
        iReg row = A0k_idx[n];
        if (ws_val[row] == 0.0) {
            ws_idx[ws_count++] = row;
        }
        ws_val[row] -= A0k[row]; 
    }

    // Populate res, generate L, and compute total sq_sum
    double sq_sum = 0.0;
    iReg local_usedL = 0;

    for (int n = 0; n < ws_count; ++n) {
        iReg row = ws_idx[n];
        double final_res = ws_val[row];

        // If the values didn't perfectly cancel out to 0, it belongs in L
        if (final_res != 0.0) {
            res[row] = final_res;
            L[local_usedL++] = row;
            sq_sum += final_res * final_res;
        }

        // Reset workspace to 0.0 for the next call
        ws_val[row] = 0.0;
    }

    // Update the reference with the new count of non-zeros
    usedL = local_usedL;

    // Final Global Math
    double normAjMh = std::sqrt(normAjMh_sq);
    resNorm = std::sqrt(sq_sum);
    std::sort(L,L+usedL);

    resRelNorm = 2.0 * resNorm / (normAjMh + resRelNorm);
}
