#include "cpt_resRho.h"
#include "precision.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm> // For std::max
#include "lapack.h"
#include "blas.h"

// // 1. Declare the Fortran BLAS functions using extern "C"
// extern "C" {
//     void dgemm_(const char* transa, const char* transb, const ptrdiff_t* m, const ptrdiff_t* n, const ptrdiff_t* k,
//                 const double* alpha, const double* a, const ptrdiff_t* lda, const double* b, const ptrdiff_t* ldb,
//                 const double* beta, double* c, const ptrdiff_t* ldc);
//     double dnrm2_(const ptrdiff_t* n, const double* x, const ptrdiff_t* incx);
//     void daxpy_(const ptrdiff_t* n, const double* alpha, const double* x, const ptrdiff_t* incx, double* y, const ptrdiff_t* incy);
// }

void cptRes(ptrdiff_t nn_A, ptrdiff_t sizeJ, double *A0k, double *AJ, double *mHat, double *res, double &resRelNorm, double &resNorm){
   // 2. RowMajor to ColumnMajor conversion for dgemm:
   // Original: cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nn_A, 1, sizeJ, 1.0, AJ, sizeJ, mHat, 1, 0.0, res, 1);
   // Swapping the matrices handles the row-major layout using Fortran's column-major logic.
   char transa = 'N';
   char transb = 'N';
   ptrdiff_t m = 1;
   ptrdiff_t n = nn_A;
   ptrdiff_t k = sizeJ;
   double alpha = 1.0;
   double beta = 0.0;
   ptrdiff_t lda = 1;
   ptrdiff_t ldb = sizeJ;
   ptrdiff_t ldc = 1;

   dgemm_(&transa, &transb, &m, &n, &k, &alpha, mHat, &lda, AJ, &ldb, &beta, res, &ldc);

   // 3. Compute the norm of Aj*mHat
   ptrdiff_t incx = 1;
   double normAjMh = dnrm2_(&n, res, &incx);

   // 4. Compute res = Aj*mHat - A0k and save it in res
   double alpha_axpy = -1.0;
   ptrdiff_t incy = 1;
   daxpy_(&n, &alpha_axpy, A0k, &incx, res, &incy);

   // 5. Compute the residual norm
   resNorm = dnrm2_(&n, res, &incx);

   // Compute the relative version of the norm
   resRelNorm = 2*resNorm/(normAjMh+resRelNorm);

   return;
}


void cptRhoJ2(ptrdiff_t JtildeSize, double *normColJ, double *AJtilde, ptrdiff_t nn_A, double *res, double *tmpRes, double normRes){
   ptrdiff_t incx = 1;
   ptrdiff_t k = nn_A;
   // Compute the norm for each column. 
   // The matrix is saved in column major so doing the norm is easy
   for (ptrdiff_t i = 0; i < JtildeSize; ++i){
      normColJ[i] = dnrm2_(&k, &(AJtilde[i*nn_A]), &incx);
      normColJ[i] *= normColJ[i];

      // If the norm is zero then discard this column
      if (normColJ[i] == 0.){
         normColJ[i] = 1e25;
      }
   }

   // Copy the residual inside the temporary vector
   std::memcpy(tmpRes, res, nn_A*sizeof(double));

   // 6. ColumnMajor dgemm maps directly, just switch to pointers
   // Original: cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, JtildeSize, 1, nn_A, 1.0, AJtilde, nn_A, tmpRes, nn_A, 0.0, res, nn_A);
   char transa = 'T';
   char transb = 'N';
   ptrdiff_t m = JtildeSize;
   ptrdiff_t n = 1;
   double alpha = 1.0;
   double beta = 0.0;
   ptrdiff_t lda = nn_A;
   ptrdiff_t ldb = nn_A;
   ptrdiff_t ldc = nn_A;

   dgemm_(&transa, &transb, &m, &n, &k, &alpha, AJtilde, &lda, tmpRes, &ldb, &beta, res, &ldc);

   double temp;
   // Compute the rhoJ2 and save it in normColJ
   for (ptrdiff_t i = 0; i < JtildeSize; ++i){
      temp = res[i];
      normColJ[i] = normRes*normRes - temp*temp/normColJ[i];
      // Avoid having possibly negative rhos 
      normColJ[i] = std::max(normColJ[i],0.);
   }

   return;
}

// Compute the index in which rhoJ2 is minimum
ptrdiff_t minIdx(double *rhoJ2, ptrdiff_t JtildeSize){
   ptrdiff_t idx = 0;
   // Initialize the minimum
   double minimum = rhoJ2[0];

   // Loop over the rhoJ2 to get the minimum
   for (ptrdiff_t i = 1; i < JtildeSize; ++i){
      if (rhoJ2[i] < minimum){
         idx = i;
         minimum = rhoJ2[i];
      }
   }

   return idx;
}