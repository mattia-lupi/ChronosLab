
#include <vector>
#include <cstring>
#include <iostream>
#include <numeric>
#include "lapack.h"
#include "cpt_sam_adaptive_left.h"
#include "qr_functions.h"
#include "find_stuff.h"
#include "cpt_resRho.h"
#include "omp.h"

bool debug = false;
iReg checkCol = 0;
iReg checkLevel = 0;

// blas Fortran routines declarations
extern "C" {
    void dgemm_(const char* transa, const char* transb, const lapack_int* m, const lapack_int* n, const lapack_int* k,
                const double* alpha, const double* a, const lapack_int* lda, const double* b, const lapack_int* ldb,
                const double* beta, double* c, const lapack_int* ldc);
    double dnrm2_(const lapack_int* n, const double* x, const lapack_int* incx);
    void daxpy_(const lapack_int* n, const double* alpha, const double* x, const lapack_int* incx, double* y, const lapack_int* incy);
}



void print_matrix(const char* desc, iReg m, iReg n, double* mat, iReg lda) {
    printf("\n--- %s (%dx%d) ---\n", desc, m, n);
    for (iReg i = 0; i < m; i++) {
        for (iReg j = 0; j < n; j++) {
            printf("%10.4g ", mat[i + j * lda]);
        }
        printf("\n");
    }
}
void print_matrix(const char* desc, iReg m, iReg n, iReg* mat, iReg lda) {
    printf("\n--- %s (%dx%d) ---\n", desc, m, n);
    for (iReg i = 0; i < m; i++) {
        for (iReg j = 0; j < n; j++) {
            printf("%d ", mat[i + j * lda]);
        }
        printf("\n");
    }
}

void print_vector(const char* desc, iReg n, double* vec) {
    printf("\n--- %s ---\n", desc);
    for (iReg i = 0; i < n; i++) {
        printf("%10.4g\n", vec[i]);
    }
}
void print_vector(const char* desc, iReg n,  iReg* vec) {
    printf("\n--- %s ---\n", desc);
    for (iReg i = 0; i < n; i++) {
        printf("%d\n", vec[i]);
    }
}

void cpt_sam_adaptive_left(iExt *iatk, iReg *jak, double *coefk,
                           iExt *iat0, iReg *ja0, double *coef0,
                           iReg nthread, iReg n_step, iReg step_size, double eps, iExt nn_A,
                           iExt *&iatN, iReg *&jaN, double *&coefN, double &avg_resRelNorm){

   // Allocate the space for storing the outputs of the loop over the columns
   std::vector<iReg> storageJsizevec(nn_A);
   iReg *storageJsize = storageJsizevec.data();
   std::vector<iReg> storageJvec(n_step*nn_A);
   iReg *storageJ = storageJvec.data();
   std::vector<double> storageNvec(n_step*nn_A);
   double *storageN = storageNvec.data();
   avg_resRelNorm = 0;

   iReg maxJsize = n_step*step_size;

   // Cycle over the columns
   #pragma omp parallel for num_threads(nthread)
   for (iReg k = 0; k < nn_A; ++k){
      // All allocation could be done outside the cycle multiplying the size of the memory
      // needed by the size of the parallel thread pool. Then each thread can read and write
      // only on its part of the allocated memory. in this way the number of allocations is 
      // reduced to 1 per object instead of being nn_A while also reducing the actual max size
      // allocated

      double resRelNorm = 1.0, resNorm;
      iReg usedL, JtildeSize;

      // Allocate pointers to where the new factorization starts
      std::vector<lapack_int> qStartvec(n_step+1);
      lapack_int *qStart = qStartvec.data();
      qStart[0] = 0;

      // Allocate pointers to how big was I at each step
      std::vector<lapack_int> sizeIvec(n_step+1);
      lapack_int *sizeI = sizeIvec.data();
      sizeI[0] = 0;

      // Allocate pointers to how big was J at each step
      std::vector<lapack_int> sizeJvec(n_step+2);
      lapack_int *sizeJ = sizeJvec.data();
      sizeJ[0] = 0;
      sizeJ[1] = 1;

      // Allocate J
      std::vector<iReg> Jvec(maxJsize);
      iReg *J = Jvec.data();
      // Set first sparsity pattern to be diagonal
      J[0] = k;

      // Allocate Jtilde
      std::vector<iReg> Jtildevec(nn_A);
      iReg *Jtilde = Jtildevec.data();

      // Allocate I
      std::vector<iReg> Ivec(nn_A);
      iReg *I = Ivec.data();

      // Allocate space for the vector norms
      std::vector<double> normvec(nn_A);
      double *normColJ = normvec.data();

      // Allocate space for A0(I,k)
      std::vector<double> a0kvec(nn_A);
      double *a0k = a0kvec.data();

      // Allocate space for mHat(I,k)
      std::vector<double> mHatvec(nn_A);
      double *mHat = mHatvec.data();

      // Allocate space for A0(:,k)
      std::vector<double> A0kvec(nn_A);
      double *A0k = A0kvec.data();

      // Allocate space for residuals
      std::vector<double> resvec(nn_A);
      double *res = resvec.data();

      // Allocate space for L
      std::vector<iReg> Lvec(nn_A);
      iReg *L = Lvec.data();

      // Fill the current column of the "old" matrix 
      // Do it once per column
      fullA0k(nn_A, iat0, ja0, coef0, k, A0k);
      // Compute the norm only once
      lapack_int N = static_cast<lapack_int>(nn_A);
      lapack_int one = 1;
      double normA0k = dnrm2_(&N,A0k,&one);

      // Allocate Ahat buffer for the max possible size
      std::vector<double> AhatBuffer(nn_A*maxJsize);
      double *Ahat = AhatBuffer.data();

      // Allocate AJ buffer for the max possible size
      // Contains the householder vectors for all the updates
      std::vector<double> AJBuffer(nn_A*maxJsize);
      double *AJ = AJBuffer.data();

      // Allocate AJtilde buffer for the max possible size
      std::vector<double> AJtildeBuffer(nn_A*nn_A);
      double *AJtilde = AJtildeBuffer.data();

      // Allocate tau for the max possible size
      std::vector<double> tauVec(maxJsize);
      double *tau = tauVec.data();

      // Allocate R for the max possible size
      std::vector<double> RVec(maxJsize*maxJsize);
      double *R = RVec.data();

      // Allocate R for the max possible size (only triangular values)
      std::vector<double> RtriangVec(maxJsize*maxJsize);
      double *Rtriang = RtriangVec.data();

      // Manual workspace allocation
      char side = 'L';
      char trans = 'T';
      lapack_int lwork = -1;
      lapack_int mxJ = static_cast<lapack_int>(maxJsize);
      double work_query, work_query1;
      lapack_int info;
      
      // Query workspace size
      dormqr_(&side, &trans, &N, &one, &mxJ, Ahat, &N, tau, a0k, &N, &work_query1, &lwork, &info);
      if (info != 0){
         printf("Error in allocating workspace, error %d\n",static_cast<int>(info));
      }

      // Query workspace size
      dgeqrf_(&N, &mxJ, Ahat, &N, tau, &work_query, &lwork, &info);
      if (info != 0){
         printf("Error in allocating workspace, error %d\n",static_cast<int>(info));
      }

      lwork = static_cast<lapack_int>(work_query);
      lwork = std::max(static_cast<lapack_int>(work_query1),lwork);
      // Allocate workspace for the max possible size
      std::vector<double> workVec(lwork);
      double *work = workVec.data();

      // Initialize the size of J to 1 as the number of entries
      iReg n2 = 1, n2old = 0;
      iReg oldSizeI, sizeIcurr = 0;
      iReg Astart = 0;

      for (iReg t = 0; t < n_step; ++t){
         if (k == checkCol && debug) printf("-------------------------------------------------------\n");
         if (k == checkCol && debug) printf("----------------------- %d -----------------------------\n",t);
         if (k == checkCol && debug) printf("-------------------------------------------------------\n");
         oldSizeI = sizeIcurr;

         if (k == checkCol && debug) print_vector("Vector J", n2, J);
         // print_vector("Vector iatk", nn_A+1, iatk);
         // print_vector("Vector jak", 30, jak);

         // Get the row entries of the matrix in columns J
         findNonZeroInColJ(J, iatk, jak, n2, I, sizeIcurr);
         sizeI[t+1] = sizeIcurr;
         // printf("%d -> %d\n", oldSizeI,sizeIcurr);

         // check what happens if I does not change
         if (k == checkCol && debug) print_vector("Vector I", sizeIcurr, I);

         // Get the new piece of A0k and add it to ak0
         getA0k(a0k, I, sizeIcurr, oldSizeI, iat0, ja0, coef0, k);

         // Copy into mHat to be modified later
         std::memcpy(mHat,a0k,sizeIcurr*sizeof(double));
         if (k == checkCol && debug) print_vector("Vector A0k", sizeIcurr, a0k);

         // Add to Ahat the new part of the matrix on which to work
         // printf("Astart %d, Jstart %d, Jend %d\n", Astart,n2old,n2);
         getAhat(I, sizeIcurr, J, n2old, n2, iatk, jak, coefk, Ahat, Astart);

         if (k == checkCol && debug) print_vector("init Matrix Ahat QR", qStart[t] + sizeI[t+1]*(sizeJ[t+1]- sizeJ[t]), Ahat);
         // print_vector("Vector sizeI", t+2, sizeI);
         // print_vector("Vector sizeJ", t+1, sizeJ);

         if (t == 0){
            // Compute the QR factorization of the first matrix
            computeFirstQR(Ahat, sizeIcurr, n2, R, Rtriang, tau, work, lwork, info);
            qStart[t+1] = sizeIcurr * n2;
            if (k == checkCol && debug) print_vector("cpt 1 Matrix Ahat QR", sizeIcurr*n2, Ahat);
            if (k == checkCol && debug) print_vector("tau", n2, tau);

            // Apply Qt for the first matrix
            applyFirstQt(Ahat, sizeIcurr, n2, tau, mHat, work, lwork, info);
            if (k == checkCol && debug) print_vector("Vector chat", sizeIcurr, mHat);
         } else{
            // Apply the previously computed Q to the new Ahat part
            applyQt(t, sizeJ, sizeI, qStart, Ahat, tau, Ahat + qStart[t], sizeI[t], sizeJ[t+1]-sizeJ[t], work, lwork, info);
            if (k == checkCol && debug) print_vector("apply Matrix Ahat QR", qStart[t] + sizeI[t+1]*(sizeJ[t+1]- sizeJ[t]), Ahat);

            computeNewQR(t, sizeI, sizeJ, qStart, Ahat, tau, R, Rtriang, work, lwork, info);
            if (k == checkCol && debug) print_vector("cpt n Matrix Ahat QR", qStart[t+1], Ahat);
            if (k == checkCol && debug) print_vector("tau", n2, tau);
            if (k == checkCol && debug) print_vector("Rtriang", n2*(n2+1)/2, Rtriang);
            if (k == checkCol && debug) print_vector("Vector chat", sizeI[t+1], mHat);
            applyQt(t+1, sizeJ, sizeI, qStart, Ahat, tau, mHat, sizeIcurr, 1, work, lwork, info);
            if (k == checkCol && debug) print_vector("Vector chat", sizeI[t+1], mHat);
            // print_matrix("Matrix R", n2, n2, R, n2);
         }
         // print_vector("qStart", t+2, qStart);
         // Solve the triangolar system
         applyR(n2, R, mHat, info);
         if (k == checkCol && debug) print_matrix("Matrix R", n2, n2, R, n2);
         if (k == checkCol && debug) print_vector("Vector mhat", n2, mHat);

         // Get the matrix A(:,J)
         getAJ(J, n2, nn_A, iatk, jak, coefk, AJ);
         // if (k == checkCol && debug) print_vector("Matrix Aj", nn_A*n2, AJ);

         // Assign the norm of A0(:,k) to resRelNorm
         resRelNorm = normA0k;
         // Compute res = AJ*mHat-A0k and save it in AJ
         cptRes(nn_A, n2, A0k, AJ, mHat, res, resRelNorm, resNorm);
         // if (k == checkCol && debug) print_vector("vector res", std::min(nn_A,20), res);

         if (resRelNorm < eps){
            break;
         }

         if (t < n_step - 1){
            // Find the values for L
            fillL(L, res, nn_A, usedL);
            // if (k == checkCol && debug) print_vector("L", usedL, L);

            // Find the values for Jtilde
            findJtilde(Jtilde, JtildeSize, L, usedL, iatk, jak, J, n2);
            /// if (k == checkCol && debug) print_vector("Vector Jtilde", JtildeSize, Jtilde);

            // Nothing new to add
            if (JtildeSize == 0){
               break;
            }

            // If there is more than one need to decide which to add
            if (JtildeSize != 1){
               // Find A(:,Jtilde)
               fullAJtilde(nn_A, iatk, jak, coefk, Jtilde, JtildeSize, AJtilde);
               // if (k == checkCol) print_vector("Matrix Ajtilde", nn_A*JtildeSize, AJtilde);//print_matrix("Matrix Ajtilde", nn_A, JtildeSize, AJtilde, nn_A);
   
               // Compute rhoJ2 = norm(res)^2 - (rTA.^2 ./ sum_A2) and save it in normColJ
               cptRhoJ2(JtildeSize, normColJ, AJtilde, nn_A, res, mHat, resNorm);
               if (k == checkCol && debug) print_vector("Vector rhoJ2", JtildeSize, normColJ);
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
                  //return;
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
         if (checkLevel == t && k == checkCol && debug){
            return;
         }
      }

      // Compute the average relative residual norm
      #pragma omp atomic 
      avg_resRelNorm += resRelNorm;

      printf("col %d, avgRes %.2g, t %d\n",k,resRelNorm,n2);

      // Copy the results in the storage
      std::memcpy(&(storageJ[k*n_step]),J,n2*sizeof(iReg));
      std::memcpy(&(storageN[k*n_step]),mHat,n2*sizeof(double));
      storageJsize[k] = n2;

      if (k == checkCol && debug){
         break;
      }
   }

   // Take the average dividing by the number of rows once
   avg_resRelNorm /= nn_A;

   // Compute total nnz of the SAM matrix
   iReg total_nnz = std::accumulate(storageJsizevec.begin(), storageJsizevec.end(), 0.0);
   // printf("%d\n", total_nnz);

   // Allocate the space for the SAM matrix
   iatN = new iReg[nn_A+1];
   jaN = new iReg[total_nnz];
   coefN = new double[total_nnz];

   // print_vector("storageJsize",nn_A,storageJsize);
   // print_matrix("storageJ",n_step, nn_A, storageJ, n_step);
   // print_matrix("storageN",n_step, nn_A, storageN, n_step);

   // Set values
   iatN[0] = 0;
   iReg nEnt, count = 0, entryPos;

   // Loop over all columns
   for (iReg i = 0; i < nn_A; ++i){
      nEnt = storageJsize[i];
      iatN[i+1] = iatN[i] + nEnt;
      // Loop over the max number of entries for each column
      for (iReg j = 0; j < nEnt; ++j){
         // Get where to save the column index or the value
         entryPos = i*n_step + j;

         // Save column indices
         jaN[count] = storageJ[entryPos];

         // Save column values
         coefN[count] = storageN[entryPos];
         count++;
      }
   }

   // print_vector("iatN",nn_A+1,iatN);
   // print_vector("jaN",total_nnz,jaN);
   // print_vector("iatN",total_nnz,coefN);

   return;
}







