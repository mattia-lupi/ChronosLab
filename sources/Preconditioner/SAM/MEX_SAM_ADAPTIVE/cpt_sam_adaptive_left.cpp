
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
#include "transpose.h"

#define debug false
iReg checkCol = 0;
iReg checkLevel = 52;

// blas Fortran routines declarations
extern "C" {
    void dgemm_(const char* transa, const char* transb, const lapack_int* m, 
                const lapack_int* n, const lapack_int* k, const double* alpha, 
                const double* a, const lapack_int* lda, const double* b, 
                const lapack_int* ldb, const double* beta, double* c, 
                const lapack_int* ldc);

    double dnrm2_(const lapack_int* n, const double* x, const lapack_int* incx);
    
    void daxpy_(const lapack_int* n, const double* alpha, const double* x, 
                const lapack_int* incx, double* y, const lapack_int* incy);
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
        printf("%10.5g\n", vec[i]);
    }
}
void print_vector(const char* desc, iReg n,  iReg* vec) {
    printf("\n--- %s ---\n", desc);
    for (iReg i = 0; i < n; i++) {
        printf("%d\n", vec[i]);
    }
}

void print_spVec(const char* desc, iReg n, double* vec) {
    printf("\n--- %s (Nonzero Elements) ---\n", desc);
    for (iReg i = 0; i < n; i++) {
        if (vec[i] != 0.0) {
            printf("Id [%ld]: %10.4g\n", (long)i, vec[i]);
        }
    }
}

void cpt_sam_adaptive_left(iExt *iatk, iReg *jak, double *coefk,
                           iExt *iat0, iReg *ja0, double *coef0,
                           iReg nthread, iReg n_step, iReg step_size, double eps, iExt nn_A,
                           iExt *&iatN, iReg *&jaN, double *&coefN, double &avg_resRelNorm){

   // Compute the number of nonzeros
   iExt nnzAk = iatk[nn_A];

   // Get the transpose to do stuff much faster
   std::vector<iExt> iatkTvec(nn_A+1,0.);
   std::vector<iReg> jakTvec(nnzAk,0.);
   std::vector<double> coefkTvec(nnzAk,0.);
   iExt *iatkT = iatkTvec.data();
   iReg *jakT = jakTvec.data();
   double *coefkT = coefkTvec.data();

   // Transpose the matrix
   int info1 = transpose(nn_A,iatk,jak,coefk,iatkT,jakT,coefkT);
   if (info1 != 0 ){
      printf("error in tranposing matrix\n");
      return;
   }

   // Allocate the space for storing the outputs of the loop over the columns
   std::vector<iReg> storageJsizevec(nn_A);
   iReg *storageJsize = storageJsizevec.data();
   std::vector<iReg> storageJvec(n_step*nn_A);
   iReg *storageJ = storageJvec.data();
   std::vector<double> storageNvec(n_step*nn_A);
   double *storageN = storageNvec.data();
   avg_resRelNorm = 0;

   iReg maxJsize = n_step*step_size;

   // Compute the values most used to avoid useless computations
   iExt tmpVal = (n_step+1)*nthread;

   // Allocate pointers to where the new factorization starts
   std::vector<lapack_int> qStartvec(tmpVal);

   // Allocate pointers to how big was I at each step
   std::vector<lapack_int> sizeIvec(tmpVal);

   // Allocate pointers to how big was J at each step
   std::vector<lapack_int> sizeJvec((n_step+2)*nthread);

   // Allocate Jtilde
   tmpVal = nn_A*nthread;
   std::vector<iReg> Jtildevec(tmpVal);

   // Allocate I
   std::vector<iReg> Ivec(tmpVal);

   // Allocate J
   std::vector<iReg> Jvec(maxJsize*nthread);

   // Allocate space for the vector norms
   std::vector<double> normvec(tmpVal);

   // Allocate space for A0(I,k)
   std::vector<double> a0kvec(tmpVal);

   // Allocate space for mHat(I,k)
   std::vector<double> mHatvec(tmpVal);

   // Allocate space for A0(:,k)
   std::vector<double> A0kvec(tmpVal);

   // Allocate space for residuals
   std::vector<double> resvec(tmpVal);

   // Allocate space for L
   std::vector<iReg> Lvec(tmpVal);

   // Allocate Ahat buffer for the max possible size
   tmpVal *= maxJsize;
   std::vector<double> AhatBuffer(tmpVal);

   // Allocate AJ buffer for the max possible size
   // Allocate it for csc matrix
   std::vector<iExt> jatAJvec((maxJsize+1)*nthread);
   std::vector<iReg> iaAJvec(nnzAk*nthread);
   std::vector<double> coefAJvec(nnzAk*nthread);

   // Allocate AJtilde buffer for the max possible size
   // Allocate it for csc matrix
   std::vector<iExt> jatAtildeVec((nn_A+1)*nthread);
   std::vector<iReg> iaAtildeVec(nnzAk*nthread);
   std::vector<double> coefAtildeVec(nnzAk*nthread);

   // Allocate tau for the max possible size
   std::vector<double> tauVec(maxJsize*nthread);

   // Allocate R for the max possible size
   tmpVal = maxJsize*maxJsize*nthread;
   std::vector<double> RVec(tmpVal);

   // Allocate R for the max possible size (only triangular values)
   std::vector<double> RtriangVec(maxJsize*(maxJsize+1)*nthread/2);

   // Manual workspace allocation
   char side = 'L';
   char trans = 'T';
   lapack_int N = static_cast<lapack_int>(nn_A);
   lapack_int one = 1;
   lapack_int lwork = -1;
   lapack_int mxJ = static_cast<lapack_int>(maxJsize);
   double work_query, work_query1;
   lapack_int info;
   
   // Query workspace size
   #if defined(__clang__) && !defined(MATLAB_MEX_FILE)
      dormqr_(&side, &trans, &N, &one, &mxJ, AhatBuffer.data(), &N, 
              tauVec.data(), a0kvec.data(), &N, &work_query1, &lwork, &info,1,1);
   #else
      dormqr_(&side, &trans, &N, &one, &mxJ, AhatBuffer.data(), &N, 
              tauVec.data(), a0kvec.data(), &N, &work_query1, &lwork, &info);
   #endif
   if (info != 0){
      printf("Error in allocating workspace, error %d\n",static_cast<int>(info));
   }

   // Query workspace size
   dgeqrf_(&N, &mxJ, AhatBuffer.data(), &N, tauVec.data(), &work_query, &lwork, &info);
   if (info != 0){
      printf("Error in allocating workspace, error %d\n",static_cast<int>(info));
   }

   lwork = static_cast<lapack_int>(work_query);
   lwork = std::max(static_cast<lapack_int>(work_query1),lwork);
   // Allocate workspace for the max possible size
   std::vector<double> workVec(lwork*nthread);

   // Cycle over the columns
   #pragma omp parallel for num_threads(nthread)
   for (iReg k = 0; k < nn_A; ++k){
      double resRelNorm = 1.0, resNorm;
      iReg usedL, JtildeSize;

      // Get local id number
      iReg thId = omp_get_thread_num();

      // Compute most used values to avoid useless computations
      iExt tmpLocVal = (n_step+1)*thId;
      
      // Get the local point for the data
      lapack_int *qStart = qStartvec.data() + tmpLocVal;
      qStart[0] = 0;
      lapack_int *sizeI = sizeIvec.data() + tmpLocVal;
      sizeI[0] = 0;

      lapack_int *sizeJ = sizeJvec.data() + (n_step+2)*thId;
      sizeJ[0] = 0;
      sizeJ[1] = 1;

      // Get the local point for the data
      tmpLocVal = nn_A*thId;
      iReg *Jtilde = Jtildevec.data() + tmpLocVal;
      iReg *I = Ivec.data() + tmpLocVal;
      double *normColJ = normvec.data() + tmpLocVal;
      double *a0k = a0kvec.data() + tmpLocVal;
      double *mHat = mHatvec.data() + tmpLocVal;
      double *A0k = A0kvec.data() + tmpLocVal;
      double *res = resvec.data() + tmpLocVal;
      iReg *L = Lvec.data() + tmpLocVal;

      // Get the local point for the data
      tmpLocVal *= maxJsize;
      double *Ahat = AhatBuffer.data() + tmpLocVal;

      // Get local pointer for AJ csc
      iExt   *jatAJ  = jatAJvec.data() + (maxJsize+1)*thId;
      iReg   *iaAJ   =  iaAJvec.data() + nnzAk*thId;
      double *coefAJ = coefAJvec.data()+ nnzAk*thId;

      // Get local pointer for Atilde csc
      iExt   *jatAJtilde  = jatAtildeVec.data() + (nn_A+1)*thId;
      iReg   *iaAJtilde   =  iaAtildeVec.data() + nnzAk*thId;
      double *coefAJtilde = coefAtildeVec.data()+ nnzAk*thId;

      double *tau = tauVec.data() + maxJsize*thId;
      double *work = workVec.data() + lwork*thId;

      // Get the local point for the data
      tmpLocVal = maxJsize*maxJsize*thId;
      double *R = RVec.data() + tmpLocVal;
      double *Rtriang = RtriangVec.data() + tmpLocVal;

      iReg *J = Jvec.data() + maxJsize*thId;
      // Set first sparsity pattern to be diagonal
      J[0] = k;

      // Fill the current column of the "old" matrix 
      // Do it once per column
      fullA0k(nn_A, iat0, ja0, coef0, k, A0k);
      // Compute the norm only once
      double normA0k = dnrm2_(&N,A0k,&one);

      // Initialize the size of J to 1 as the number of entries
      iReg n2 = 1, n2old = 0;
      iReg oldSizeI, sizeIcurr = 0;
      iReg Astart = 0;

      for (iReg t = 0; t < n_step; ++t){
         oldSizeI = sizeIcurr;

         // Get the row entries of the matrix in columns J
         findNonZeroInColJ(J, iatk, jak, n2, I, sizeIcurr);
         sizeI[t+1] = sizeIcurr;

         // Debug prints
         #if debug
         if (k == checkCol) printf("-------------------------------------------------------\n");
         if (k == checkCol) printf("----------------------- %d -----------------------------\n",t);
         if (k == checkCol) printf("-------------------------------------------------------\n");
         if (k == checkCol) print_vector("Vector J", n2, J);
         if (k == checkCol) print_vector("Vector I", sizeIcurr, I);
         #endif

         // Get the new piece of A0k and add it to ak0
         getA0k(a0k, I, sizeIcurr, oldSizeI, iat0, ja0, coef0, k);

         // Copy into mHat to be modified later
         std::memcpy(mHat,a0k,sizeIcurr*sizeof(double));
         
         // Add to Ahat the new part of the matrix on which to work
         getAhat(I, sizeIcurr, J, n2old, n2, iatk, jak, coefk, Ahat, Astart);

         // Debug prints
         #if debug
         if (k == checkCol) print_vector("Vector A0k", sizeIcurr, a0k);
         if (k == checkCol) print_spVec("init Matrix Ahat QR", qStart[t] + sizeI[t+1]*(sizeJ[t+1]- sizeJ[t]), Ahat);
         #endif

         if (t == 0){
            // Compute the QR factorization of the first matrix
            computeFirstQR(Ahat, sizeIcurr, n2, R, Rtriang, tau, work, lwork, info);
            qStart[t+1] = sizeIcurr * n2;

            // Apply Qt for the first matrix
            applyFirstQt(Ahat, sizeIcurr, n2, tau, mHat, work, lwork, info);

            // Debug prints
            #if debug
            if (k == checkCol) print_spVec("cpt 1 Matrix Ahat QR", qStart[t+1], Ahat);
            if (k == checkCol) print_vector("tau", n2, tau);
            if (k == checkCol) print_vector("Rtriang", n2*(n2+1)/2, Rtriang);
            if (k == checkCol) print_vector("Vector chat", sizeI[t+1], mHat);
            #endif
         } else{

            // Apply the previously computed Q to the new Ahat part
            applyQt(t, sizeJ, sizeI, qStart, Ahat, tau, Ahat + qStart[t], sizeI[t], sizeJ[t+1]-sizeJ[t], work, lwork, info);

            // Debug prints
            #if debug
            if (k == checkCol) print_spVec("apply Matrix Ahat QR", qStart[t] + sizeI[t+1]*(sizeJ[t+1]- sizeJ[t]), Ahat);
            #endif

            // Compute the new part of the qr_factorization
            computeNewQR(t, sizeI, sizeJ, qStart, Ahat, tau, R, Rtriang, work, lwork, info);
            
            // Apply the Qt to the rhs to get chat
            applyQt(t+1, sizeJ, sizeI, qStart, Ahat, tau, mHat, sizeIcurr, 1, work, lwork, info);
            
            // Debug prints
            #if debug
            if (k == checkCol) print_spVec("cpt n Matrix Ahat QR", qStart[t+1], Ahat);
            if (k == checkCol) print_vector("tau", n2, tau);
            if (k == checkCol) print_vector("Rtriang", n2*(n2+1)/2, Rtriang);
            if (k == checkCol) print_vector("Vector chat", sizeI[t+1], mHat);
            #endif
         }

         // Solve the triangolar system
         applyR(n2, R, mHat, info);

         // Get the matrix A(:,J)
         getAJ(J, n2old, n2, iatkT, jakT, coefkT, jatAJ, iaAJ, coefAJ);

         // Debug prints
         #if debug
         if (k == checkCol) print_matrix("Matrix R", n2, n2, R, n2);
         if (k == checkCol) print_vector("Vector mhat", n2, mHat);
         if (k == checkCol) print_vector("vector jatAJ", n2+1, jatAJ);
         if (k == checkCol) print_vector("vector iaAJ", jatAJ[n2], iaAJ);
         if (k == checkCol) print_vector("vector coefAJ", jatAJ[n2], coefAJ);
         #endif

         // Assign the norm of A0(:,k) to resRelNorm
         resRelNorm = normA0k;

         // Compute res = AJ*mHat-A0k and save it in AJ
         cptRes(nn_A, n2, A0k, jatAJ, iaAJ, coefAJ, mHat, res, resRelNorm, resNorm,iaAJtilde,coefAJtilde);
         
         // Debug prints
         #if debug
         if (k == checkCol) print_spVec("vector res", nn_A, res);
         #endif

         if (resRelNorm < eps){
            break;
         }

         if (t < n_step - 1){
            // Find the values for L
            fillL(L, res, nn_A, usedL);

            // Find the values for Jtilde
            findJtilde(Jtilde, JtildeSize, L, usedL, iatk, jak, J, n2);
            std::sort(Jtilde, Jtilde + JtildeSize);

            // Debug prints
            #if debug
            if (k == checkCol) print_vector("L", usedL, L);
            if (k == checkCol) print_vector("Vector Jtilde", JtildeSize, Jtilde);
            #endif

            // Nothing new to add
            if (JtildeSize == 0){
               break;
            }

            // If there is more than one need to decide which to add
            if (JtildeSize != 1){
               // Find A(:,Jtilde)
               getAJ(Jtilde, 0, JtildeSize, iatk, jak, coefk, jatAJtilde, iaAJtilde, coefAJtilde);

               // Compute rhoJ2 = norm(res)^2 - (rTA.^2 ./ sum_A2) and save it in normColJ
               cptRhoJ2(JtildeSize, normColJ, jatAJtilde, iaAJtilde, coefAJtilde, 
                        res, mHat, resNorm);

               // Debug prints
               #if debug
               // if (k == checkCol) print_vector("vector jatAJtilde", n2+1, jatAJtilde);
               // if (k == checkCol) print_vector("vector iaAJtilde", jatAJtilde[n2], iaAJtilde);
               // if (k == checkCol) print_vector("vector coefAJtilde", jatAJtilde[n2], coefAJtilde);
               if (k == checkCol) print_vector("Vector rhoJ2", JtildeSize, normColJ);
               #endif

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

         #if debug
         if (t == checkLevel && k == checkCol){
            break;
         }
         #endif
      }

      #if debug
      if (k == checkCol){
         break;
      }
      #endif

      // Compute the average relative residual norm
      #pragma omp atomic 
      avg_resRelNorm += resRelNorm;

      printf("col %d, avgRes %.2g, t %d\n",k,resRelNorm,n2);

      // Copy the results in the storage
      std::memcpy(&(storageJ[k*n_step]),J,n2*sizeof(iReg));
      std::memcpy(&(storageN[k*n_step]),mHat,n2*sizeof(double));
      storageJsize[k] = n2;
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

   return;
}







