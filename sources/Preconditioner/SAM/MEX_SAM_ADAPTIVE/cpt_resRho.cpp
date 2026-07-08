#include "cpt_resRho.h"
#include "precision.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm> // For std::max and std::min_element
#include <iterator>  // For std::distance
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

void cptRes(iReg nn_A, iReg sizeJ, double *A0k, double *AJ, double *mHat, double *res, double &resRelNorm, double &resNorm){
   // AJ is column-major with nn_A rows and sizeJ columns.
   // Matrix-vector multiplication computes: res = AJ * mHat
   char trans = 'N';
   lapack_int m_v = static_cast<lapack_int>(nn_A);
   lapack_int n_v = static_cast<lapack_int>(sizeJ);
   double alpha = 1.0;
   double beta = 0.0;
   lapack_int lda = m_v; 
   lapack_int incx = 1;
   lapack_int incy = 1;

   // Compute AJ * mHat
   dgemv_(&trans, &m_v, &n_v, &alpha, AJ, &lda, mHat, &incx, &beta, res, &incy);

   // Compute the norm of AJ * mHat
   lapack_int N = static_cast<lapack_int>(nn_A);
   double normAjMh = dnrm2_(&N, res, &incx);

   // Compute res = AJ * mHat - A0k
   double alpha_axpy = -1.0;
   daxpy_(&N, &alpha_axpy, A0k, &incx, res, &incy);

   // Compute the residual norm
   resNorm = dnrm2_(&N, res, &incx);

   // Compute the relative version of the norm
   resRelNorm = 2 * resNorm / (normAjMh + resRelNorm);

   return;
}

void cptRhoJ2(iReg JtildeSize, double *normColJ, double *AJtilde, iReg nn_A, double *res, double *tmpRes, double normRes){
   const double normResSq = normRes * normRes;

   // Fuse ddot, dgemv, and the final adjustments into a single pass
   for (iReg i = 0; i < JtildeSize; ++i) {
      double dot = 0.0;
      double sum = 0.0;
      const double* col = &AJtilde[i * nn_A];

      // This single inner loop streams a column into cache ONCE
      // and performs both operations simultaneously.
      for (iReg j = 0; j < nn_A; ++j) {
         double val = col[j];
         dot += val * val;
         sum += val * res[j];
      }

      if (dot == 0.0) {
         dot = 1e25;
      }

      tmpRes[i] = sum;
      
      // Compute the rhoJ2 inline
      double rho = normResSq - (sum * sum) / dot;
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