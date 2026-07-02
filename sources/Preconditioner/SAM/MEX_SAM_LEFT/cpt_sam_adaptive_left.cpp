#pragma once
#include "cpt_sam_adaptive_left.h"
#include <vector>
#include <iostream>
#include "lapack.h"
#include "cblas.h"
#include "qr_functions.h"
#include "find_stuff.h"
#include "cpt_resRho.h"

int checkCol = 1;

// LAPACK Fortran routines declarations
// extern "C" {
//    void dgeqrf_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
//    void dormqr_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* lwork, int* info);
//    void dtrtrs_(char* uplo, char* trans, char* diag, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);
// }



void print_matrix(const char* desc, int m, int n, double* mat, int lda) {
    printf("\n--- %s (%dx%d) ---\n", desc, m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.4f ", mat[i + j * lda]);
        }
        printf("\n");
    }
}
void print_matrix(const char* desc, int m, int n, int* mat, int lda) {
    printf("\n--- %s (%dx%d) ---\n", desc, m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", mat[i + j * lda]);
        }
        printf("\n");
    }
}

void print_vector(const char* desc, int n, double* vec) {
    printf("\n--- %s ---\n", desc);
    for (int i = 0; i < n; i++) {
        printf("%10.4f\n", vec[i]);
    }
}
void print_vector(const char* desc, int n,  int* vec) {
    printf("\n--- %s ---\n", desc);
    for (int i = 0; i < n; i++) {
        printf("%d\n", vec[i]);
    }
}

void cpt_sam_adaptive_left(int *iatk, int *jak,double *coefk,
                           int *iat0, int *ja0,double *coef0,
                           int nthread,int n_step,int step_size,double eps, int nn_A,
                           int *&iatN, int *&jaN,double *&coefN){


   // Suppose the pattern is symmetric
   // Count the maximum number of entries per row
   int max_nnz_r = 0;
   for (int i = 0; i < nn_A; ++i){
      max_nnz_r = std::max(max_nnz_r,iatk[i+1]-iatk[i]);
   }

   // Allocate the space for storing the outputs of the loop over the columns
   std::vector<int> storageJsizevec(nn_A);
   int *storageJsize = storageJsizevec.data();
   std::vector<int> storageJvec(n_step*nn_A);
   int *storageJ = storageJvec.data();
   std::vector<double> storageNvec(n_step*nn_A);
   double *storageN = storageNvec.data();

   int maxJsize = n_step*step_size;

   // Cycle over the columns
   for (int k = 0; k < nn_A; ++k){
      // All allocation could be done outside the cycle multiplying the size of the memory
      // needed by the size of the parallel thread pool. Then each thread can read and write
      // only on its part of the allocated memory. in this way the number of allocations is 
      // reduced to 1 per object instead of being nn_A while also reducing the actual max size
      // allocated

      double resRelNorm = 1.0, resNorm;
      int usedL, JtildeSize;

      // Allocate pointers to where the new factorization starts
      std::vector<int> qStartvec(n_step);
      int *qStart = qStartvec.data();
      qStart[0] = 0;

      // Allocate pointers to how big was I at each step
      std::vector<int> sizeIvec(n_step+1);
      int *sizeI = sizeIvec.data();
      sizeI[0] = 0;

      // Allocate pointers to how big was J at each step
      std::vector<int> sizeJvec(n_step+2);
      int *sizeJ = sizeJvec.data();
      sizeJ[0] = 0;
      sizeJ[1] = 1;

      // Allocate J
      std::vector<int> Jvec(maxJsize);
      int *J = Jvec.data();
      // Set first sparsity pattern to be diagonal
      J[0] = k;

      // Allocate Jtilde
      std::vector<int> Jtildevec(nn_A);
      int *Jtilde = Jtildevec.data();

      // Allocate I
      std::vector<int> Ivec(max_nnz_r);
      int *I = Ivec.data();

      // Allocate space for A(:,J)
      std::vector<double> aJvec(maxJsize);
      double *aJ = aJvec.data();

      // Allocate space for A(:,J)
      std::vector<double> normvec(maxJsize);
      double *normColJ = normvec.data();

      // Allocate space for A0(I,k)
      std::vector<double> a0kvec(max_nnz_r);
      double *a0k = a0kvec.data();

      // Allocate space for mHat(I,k)
      std::vector<double> mHatvec(max_nnz_r);
      double *mHat = mHatvec.data();

      // Allocate space for A0(:,k)
      std::vector<double> A0kvec(nn_A);
      double *A0k = A0kvec.data();

      // Allocate space for residuals
      std::vector<double> resvec(nn_A);
      double *res = resvec.data();

      // Allocate space for L
      std::vector<int> Lvec(nn_A);
      int *L = Lvec.data();

      // Fill the current column of the "old" matrix 
      // Do it once per column
      fullA0k(nn_A, iat0, ja0, coef0, k, A0k);
      // Compute the norm only once
      double normA0k = cblas_dnrm2(nn_A,A0k,1);

      // Allocate Ahat buffer for the max possible size
      std::vector<double> AhatBuffer(max_nnz_r*maxJsize);
      double *Ahat = AhatBuffer.data();

      // Allocate AJ buffer for the max possible size
      // Contains the householder vectors for all the updates
      std::vector<double> AJBuffer(nn_A*maxJsize);
      double *AJ = AJBuffer.data();

      // Allocate AJtilde buffer for the max possible size
      std::vector<double> AJtildeBuffer(nn_A*maxJsize);
      double *AJtilde = AJtildeBuffer.data();

      // Allocate tau for the max possible size
      std::vector<double> tauVec(maxJsize);
      double *tau = tauVec.data();

      // Allocate R for the max possible size
      std::vector<double> RVec(maxJsize*maxJsize);
      double *R = RVec.data();

      // Allocate R for the max possible size (only triangular values)
      std::vector<double> RtriangVec((maxJsize)*(maxJsize+1)/2);
      double *Rtriang = RtriangVec.data();

      // DGEQRF Workspace Query
      int lwork = -1;
      double work_query;
      int info;
      
      dgeqrf_(&max_nnz_r, &maxJsize, Ahat, &max_nnz_r, tau, &work_query, &lwork, &info);
      if (info != 0){
         printf("Error in allocating workspace, error %d\n",info);
      }

      lwork = static_cast<int>(work_query);
      // Allocate workspace for the max possible size
      std::vector<double> workVec(lwork);
      double *work = workVec.data();

      // Initialize the size of J to 1 as the number of entries
      int n2 = 1, n2old = 0;
      int oldSizeI, sizeIcurr = 0;
      int Astart = 0;

      for (int t = 0; t < n_step; ++t){
         if (k == checkCol) printf("-------------------------------------------------------\n");
         if (k == checkCol) printf("----------------------- %d -----------------------------\n",t);
         if (k == checkCol) printf("-------------------------------------------------------\n");
         oldSizeI = sizeIcurr;

         if (k == checkCol) print_vector("Vector J", n2, J);
         // print_vector("Vector iatk", nn_A+1, iatk);
         // print_vector("Vector jak", 6, jak);

         // Get the row entries of the matrix in columns J
         findNonZeroInColJ(J, iatk, jak, n2, I, sizeIcurr);
         sizeI[t+1] = sizeIcurr;
         // printf("%d -> %d\n", oldSizeI,sizeIcurr);

         // check what happens if I does not change
         if (k == checkCol) print_vector("Vector I", sizeIcurr, I);

         // Get the new piece of A0k and add it to ak0
         getA0k(a0k, I, sizeIcurr, oldSizeI, iat0, ja0, coef0, k);

         // Copy onto mHat which can be modified
         std::memcpy(mHat,a0k,sizeIcurr*sizeof(double));
         if (k == checkCol) print_vector("Vector A0k", sizeIcurr, a0k);
         // print_vector("Vector mhat", sizeIcurr, mHat);

         // Add to Ahat the new part of the matrix on which to work
         // printf("Astart %d, Jstart %d, Jend %d\n", Astart,n2old,n2);
         getAhat(I, sizeIcurr, J, n2old, n2, iatk, jak, coefk, Ahat, Astart);

         if (k == checkCol) print_vector("init Matrix Ahat QR", max_nnz_r*maxJsize, Ahat);
         // print_vector("Vector sizeI", t+2, sizeI);
         // print_vector("Vector sizeJ", t+1, sizeJ);

         if (t == 0){
            // Compute the QR factorization of the first matrix
            computeFirstQR(Ahat, sizeIcurr, n2, R, Rtriang, tau, work, lwork, info);
            qStart[t+1] = sizeIcurr * n2;
            if (k == checkCol) print_vector("cpt 1 Matrix Ahat QR", max_nnz_r*maxJsize, Ahat);
            if (k == checkCol) print_vector("tau", maxJsize, tau);

            // Apply Qt for the first matrix
            applyFirstQt(Ahat, sizeIcurr, n2, tau, mHat, work, lwork, info);
            if (k == checkCol) print_vector("Vector chat", sizeIcurr, mHat);
         } else{
            // Apply the previously computed Q to the new Ahat part
            applyQt(t, sizeJ, sizeI, qStart, Ahat, tau, Ahat + qStart[t], sizeI[t], sizeJ[t+1]-sizeJ[t], work, lwork, info);
            if (k == checkCol) print_vector("apply Matrix Ahat QR", max_nnz_r*maxJsize, Ahat);

            if (k == checkCol) printf("rowSizeB2 %d, colSizeB2 %d, startB2 %d, rowSizeB %d, oldSizeTau %d\n", 
                   sizeI[t]-sizeJ[t],sizeJ[t+1] - sizeJ[t],sizeI[t]*sizeJ[t] + sizeJ[t+1] - sizeJ[t],sizeI[t+1],sizeJ[t]);
            
            computeNewQR(t, sizeI, sizeJ, qStart, Ahat, tau, R, Rtriang, work, lwork, info);
            if (k == checkCol) print_vector("cpt n Matrix Ahat QR", max_nnz_r*maxJsize, Ahat);
            if (k == checkCol) print_vector("Rtriang", n2*(n2+1)/2, Rtriang);
            if (k == checkCol) print_vector("Vector a0k", n2+1, mHat);
            applyQt(t+1, sizeJ, sizeI, qStart, Ahat, tau, mHat, sizeI[t], sizeJ[t+1]-sizeJ[t], work, lwork, info);
            if (k == checkCol) print_vector("Vector chat", n2, mHat);
            // print_matrix("Matrix R", n2, n2, R, n2);
         }
         // print_vector("qStart", t+2, qStart);
         // Solve the triangolar system
         applyR(n2, R, mHat, info);
         if (k == checkCol) print_matrix("Matrix R", n2, n2, R, n2);
         if (k == checkCol) print_vector("Vector mhat", n2, mHat);

         // Get the matrix A(:,J)
         getAJ(J, n2, nn_A, iatk, jak, coefk, AJ);
         if (k == checkCol) print_vector("Matrix Aj", nn_A*n2, AJ);

         // Assign the norm of A0(:,k) to resRelNorm
         resRelNorm = normA0k;
         // Compute res = AJ*mHat-A0k and save it in AJ
         cptRes(nn_A, n2, A0k, AJ, mHat, res, resRelNorm, resNorm);
         if (k == checkCol) print_vector("vector res", nn_A, res);

         if (resRelNorm < eps){
            break;
         }

         if (t < n_step - 1){
            // Find the values for L
            fillL(L, res, nn_A, usedL);
            if (k == checkCol) print_vector("L", usedL, L);

            // Find the values for Jtilde
            findJtilde(Jtilde, JtildeSize, L, usedL, iatk, jak, J, n2);
            if (k == checkCol) print_vector("Vector Jtilde", JtildeSize, Jtilde);

            // Nothing new to add
            if (JtildeSize == 0){
               break;
            }

            // If there is more than one need to decide which to add
            if (JtildeSize != 1){
               // Find A(:,Jtilde)
               fullAJtilde(nn_A, iatk, jak, coefk, Jtilde, JtildeSize, AJtilde);
               if (k == checkCol) print_matrix("Matrix Ajtilde", nn_A, JtildeSize, AJtilde, nn_A);
   
               // Compute rhoJ2 = norm(res)^2 - (rTA.^2 ./ sum_A2) and save it in normColJ
               cptRhoJ2(JtildeSize, normColJ, AJtilde, nn_A, res, resNorm);
               if (k == checkCol) print_vector("Vector rhoJ2", JtildeSize, normColJ);
               // print_matrix("Matrix Ajtilde", JtildeSize, JtildeSize, AJtilde, nn_A);

               // Add a single one
               if (step_size == 1){
                  // Get the new J
                  J[n2] = Jtilde[minIdx(normColJ, JtildeSize)];

                  // Update the size of J
                  n2old = n2;
                  n2++;
                  sizeJ[t+2] = n2;
               }
               else{
                  printf("Error, not implemented yet\n");
                  return;
               }
            }
            else{
               // Add the index to J
               J[n2] = Jtilde[0];

               // Update the size of J
               n2old = n2;
               n2++;
               sizeJ[t+2] = n2;
            }
         }
      }
      memcpy(&(storageJ[k*nn_A]),J,n2*sizeof(int));
      memcpy(&(storageN[k*nn_A]),mHat,n2*sizeof(double));
      storageJsize[k] = n2;
      if (k == checkCol){
         break;
      }
   }

   print_matrix("storageJ",n_step, nn_A, storageJ, n_step);
   print_matrix("storageN",n_step, nn_A, storageN, n_step);

   return;
}







