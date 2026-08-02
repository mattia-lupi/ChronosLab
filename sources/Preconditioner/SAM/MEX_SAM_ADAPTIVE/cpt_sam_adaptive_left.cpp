#include <vector>
#include <cstring>
#include <iostream>
#include <numeric>
#include <algorithm>
#include "lapack.h"
#include "cpt_sam_adaptive_left.h"
#include "qr_functions.h"
#include "find_stuff.h"
#include "cpt_resRho.h"
#include "omp.h"

#define debug false
iReg checkCol = 0;
iReg checkLevel = 2;

// BLAS Fortran routines declarations
extern "C" {
   double dnrm2_(const lapack_int* n, const double* x, const lapack_int* incx);
}
void print_matrix(const char* desc, iReg m, iReg n, double* mat, iReg lda);
void print_matrix(const char* desc, iReg m, iReg n, iReg* mat, iReg lda);
void print_vector(const char* desc, iReg n, double* vec);
void print_vector(const char* desc, iReg n, iReg* vec);
void print_spVec(const char* desc, iReg n, double* vec);
void print_spVec(const char* desc, iReg n, int* vec);

void cpt_sam_adaptive_left(iExt *iatk, iReg *jak, double *coefk, double *coefkT,
                           iExt *iat0, iReg *ja0, double *coef0,
                           iReg nthread, iReg n_step, iReg step_size, double eps, iExt nn_A,
                           iExt *&iatN, iReg *&jaN, double *&coefN, double &avg_resRelNorm) {

   iExt *iatkT = iatk;
   iReg *jakT = jak;
   iReg dCol = 0, dCol0 = 0;

   std::vector<double> colANormVec(nn_A);
   double *colANorm = colANormVec.data();

   #pragma omp parallel num_threads(nthread)
   {
      iReg local_d_col = 0;
      iReg local_d0_col = 0;
   
      #pragma omp for schedule(static)
      for (iReg i = 0; i < nn_A; i++) {
         // Get values for the current column
         int startCol = iatkT[i];
         lapack_int size = iatkT[i + 1] - startCol;
         if (size > local_d_col) {
            local_d_col = size;
         }

         // Compute column norms in parallel
         lapack_int one = 1;
         colANorm[i] = dnrm2_(&size, &(coefkT[startCol]), &one);
   
         // Compute column norms in parallel
         iReg size0 = static_cast<iReg>(iat0[i + 1] - iat0[i]);
         if (size0 > local_d0_col) {
            local_d0_col = size0;
         }
      }
   
      // Compute size useful for various vectors allocation
      #pragma omp critical
      {
         if (local_d_col > dCol) 
            dCol = local_d_col;
         if (local_d0_col > dCol0) 
            dCol0 = local_d0_col;
      }
   }

   // Allocate space for storing column output metadata
   std::vector<iReg> storageJsizeVec(nn_A);
   iReg *storageJsize = storageJsizeVec.data();
   std::vector<iReg> storageJVec(n_step * nn_A);
   iReg *storageJ = storageJVec.data();
   std::vector<double> storageNVec(n_step * nn_A);
   double *storageN = storageNVec.data();

   avg_resRelNorm = 0.0;
   iReg maxJsize = n_step * step_size;

   // LAPACK Workspace Query (performed once before parallel processing)
   iReg maxISize = std::min(nn_A,maxJsize*dCol);
   iReg maxLSize = std::min(nn_A,maxJsize*(dCol+dCol0));
   iReg maxJtildeSize = std::min(nn_A,maxLSize*dCol);
   char side = 'L';
   char trans = 'T';
   lapack_int N = static_cast<lapack_int>(maxISize);
   lapack_int one = 1;
   lapack_int lwork = -1;
   lapack_int mxJ = static_cast<lapack_int>(maxJsize);
   double work_query, work_query1;
   lapack_int info;

   double dummy = 0.0;
   #if defined(__clang__) && !defined(MATLAB_MEX_FILE)
   dormqr_(&side, &trans, &N, &one, &mxJ, &dummy, &N,
           &dummy, &dummy, &N, &work_query1, &lwork, &info, 1, 1);
   #else
   dormqr_(&side, &trans, &N, &one, &mxJ, &dummy, &N,
           &dummy, &dummy, &N, &work_query1, &lwork, &info);
   #endif
   if (info != 0) {
      printf("Error in allocating workspace, error %d\n", static_cast<int>(info));
   }

   dgeqrf_(&N, &mxJ, &dummy, &N, &dummy, &work_query, &lwork, &info);
   if (info != 0) {
      printf("Error in allocating workspace, error %d\n", static_cast<int>(info));
   }

   lwork = static_cast<lapack_int>(work_query);
   lwork = std::max(static_cast<lapack_int>(work_query1), lwork);

   #pragma omp parallel num_threads(nthread) reduction(+:avg_resRelNorm)
   {
      // Allocate the arrays
      std::vector<lapack_int> qStartVec(n_step + 1);
      std::vector<lapack_int> sizeIVec(n_step + 1);
      std::vector<lapack_int> sizeJVec(n_step + 2);

      std::vector<iReg> JtildeVec(maxJtildeSize);
      std::vector<iReg> IVec(maxISize);
      std::vector<iReg> JVec(maxJsize);
      std::vector<iReg> LVec(maxLSize);
      std::vector<double> AhatVec(maxISize * maxJsize);
      
      std::vector<double> normColJVec(nn_A);
      std::vector<double> a0kVec(nn_A);
      std::vector<double> mHatVec(nn_A);
      std::vector<double> A0kVec(nn_A);
      std::vector<iReg> A0k_idxVec(nn_A);
      std::vector<uint8_t> JtildeWorkVec(nn_A);
      std::vector<iReg> findIWorkVec(nn_A);
      std::vector<double> resVec(nn_A);
      std::vector<iReg> resIdxVec(nn_A);
      std::vector<double> resWSVec(nn_A);
      
      std::vector<double> tauVec(maxJsize);
      std::vector<double> RVec(maxJsize * maxJsize);
      std::vector<double> RtriangVec(maxJsize * (maxJsize + 1));
      std::vector<iExt> jatAJVec(maxJsize + 1);
      std::vector<const iReg*> iaAJVec(maxJsize);
      std::vector<const double*> coefAJVec(maxJsize);

      std::vector<iExt> jatAtildeVec(maxJtildeSize + 1);
      std::vector<const iReg*> iaAtildeVec(maxJtildeSize);
      std::vector<const double*> coefAtildeVec(maxJtildeSize);
      std::vector<double> workVec(lwork);

      // Get handles
      lapack_int *qStart = qStartVec.data();
      lapack_int *sizeI = sizeIVec.data();
      lapack_int *sizeJ = sizeJVec.data();
      iReg *Jtilde = JtildeVec.data();
      iReg *I = IVec.data();
      iReg *J = JVec.data();
      double *normColJ = normColJVec.data();
      double *a0k = a0kVec.data();
      double *mHat = mHatVec.data();
      double *A0k = A0kVec.data();

      iReg *A0k_idx = A0k_idxVec.data();
      iReg A0k_nnz = 0;

      uint8_t *JtildeWork = JtildeWorkVec.data();
      iReg *findIWork = findIWorkVec.data();
      double *res = resVec.data();
      iReg *resIdx = resIdxVec.data();
      double *resWS = resWSVec.data();
      iReg *L = LVec.data();
      double *Ahat = AhatVec.data();
      iExt *jatAJ = jatAJVec.data();
      const iReg **iaAJ = iaAJVec.data();
      const double **coefAJ = coefAJVec.data();
      iExt *jatAtilde = jatAtildeVec.data();
      const iReg **iaAtilde = iaAtildeVec.data();
      const double **coefAtilde = coefAtildeVec.data();
      double *tau = tauVec.data();
      double *R = RVec.data();
      double *Rtriang = RtriangVec.data();
      double *work = workVec.data();

      // Dynamically assign columns to available threads
      #pragma omp for schedule(dynamic)
      for (iReg k = 0; k < nn_A; ++k) {
         double resRelNorm = 1.0, resNorm;
         iReg usedL = 0, JtildeSize = 0;

         qStart[0] = 0;
         sizeI[0] = 0;
         sizeJ[0] = 0;
         sizeJ[1] = 1;

         J[0] = k;

         A0k_nnz = 0;
         double normA0k = fullA0k(iat0, ja0, coef0, k, A0k, A0k_idx, A0k_nnz);

         iReg n2 = 1, n2old = 0;
         iReg oldSizeI, sizeIcurr = 0;
         iReg Astart = 0;

         for (iReg t = 0; t < n_step; ++t) {
            oldSizeI = sizeIcurr;

            // Find nonzeros in the Jth column and save them in I appending on the old I
            findNonZeroInColJ(J, iatk, jak, n2, I, sizeIcurr, findIWork, k + 1);
            sizeI[t + 1] = sizeIcurr;

            #if debug
            if (k == checkCol) printf("-------------------------------------------------------\n");
            if (k == checkCol) printf("----------------------- %d -----------------------------\n", t);
            if (k == checkCol) printf("-------------------------------------------------------\n");
            if (k == checkCol) print_vector("Vector J", n2, J);
            #endif

            // Append the new entries of a0k in the vector. Then overwrite mHat
            getA0k(a0k, I, sizeIcurr, oldSizeI, iat0, ja0, coef0, k);
            std::memcpy(mHat, a0k, sizeIcurr * sizeof(double));

            // Get the Ahat matrix in the new position to prepare for the QR
            getAhat(I, sizeIcurr, J, n2old, n2, iatk, jak, coefk, Ahat, Astart);

            // Compute and apply correctly QR factorization for the current Ahat
            if (t == 0) {
               computeFirstQR(Ahat, sizeIcurr, n2, R, Rtriang, tau, work, lwork, info);
               qStart[t + 1] = sizeIcurr * n2;
               applyFirstQt(Ahat, sizeIcurr, n2, tau, mHat, work, lwork, info);

               #if debug
               if (k == checkCol) print_vector("tau", n2, tau);
               if (k == checkCol) print_vector("Rtriang", n2 * (n2 + 1) / 2, Rtriang);
               #endif
            } else {
               // Apply Q transposed in the correct place
               applyQt(t, sizeJ, sizeI, qStart, Ahat, tau, Ahat + qStart[t], sizeI[t],
                       sizeJ[t + 1] - sizeJ[t], work, lwork, info);
               // Compute new QR factor 
               computeNewQR(t, sizeI, sizeJ, qStart, Ahat, tau, R, Rtriang, work, lwork,
                            info);

               // Apply the new Q transpose
               applyQt(t + 1, sizeJ, sizeI, qStart, Ahat, tau, mHat, sizeIcurr, 1, work,
                       lwork, info);

               #if debug
               if (k == checkCol) print_vector("tau", n2, tau);
               if (k == checkCol) print_vector("Rtriang", n2 * (n2 + 1) / 2, Rtriang);
               #endif
            }

            // Apply R
            applyR(n2, R, mHat, info);

            // Get AJ pointers
            getAJ(J, n2old, n2, iatkT, jakT, coefkT, iaAJ, coefAJ, jatAJ);

            #if debug
            if (k == checkCol) print_matrix("Matrix R", n2, n2, R, n2);
            if (k == checkCol) print_vector("Vector mhat", n2, mHat);
            #endif

            // Initialize the norm
            resRelNorm = normA0k;

            // Compute residuals of the current iteration
            cptRes(n2, jatAJ, iaAJ, coefAJ, mHat, A0k_idx, A0k, A0k_nnz,
                   res, L, usedL, resRelNorm, resNorm, resIdx, resWS);

            // Check for convergence
            if (resRelNorm < eps) {
               break;
            }

            // If not the last step continue with the algorithm
            if (t < n_step - 1) {
               // Find possible candidates
               findJtilde(Jtilde, JtildeSize, L, usedL, iatk, jak, J, n2, JtildeWork);
               // Sort the candidates
               std::sort(Jtilde, Jtilde + JtildeSize);

               // If there are no candidates exit
               if (JtildeSize == 0) {
                  break;
               }

               if (JtildeSize != 1) {
                  // If there is more than one candidate get the pointers to the columns
                  getAJ(Jtilde, 0, JtildeSize, iatkT, jakT, coefkT, iaAtilde, coefAtilde,
                        jatAtilde);

                  // Compute the residual reductions
                  cptRhoJ2(JtildeSize, normColJ, jatAtilde, iaAtilde, coefAtilde,
                           res, resNorm, colANorm, Jtilde);

                  #if debug
                  if (k == checkCol) print_spVec("Vector rhoJ2", JtildeSize, normColJ);
                  #endif

                  if (step_size == 1) {
                     // If there is only one to add get the minimum index
                     J[n2] = Jtilde[minIdx(normColJ, JtildeSize)];

                     // Update sizes
                     n2old = n2;
                     n2++;
                     sizeJ[t + 2] = n2;
                  } else {
                     printf("Error, not implemented yet\n");
                  }
               } else {
                  // There is only one entry
                  J[n2] = Jtilde[0];

                  // Update sizes
                  n2old = n2;
                  n2++;
                  sizeJ[t + 2] = n2;
               }
            }

            #if debug
            if (t == checkLevel && k == checkCol) {
                break;
            }
            #endif
         }

         #if debug
         if (k == checkCol) {
             break;
         }
         #endif

         // Column average norm
         #pragma omp atomic
         avg_resRelNorm += resRelNorm;

         printf("col %d, avgRes %.2g, t %d\n", k, resRelNorm, n2);

         // Copy results to storage
         std::memcpy(&(storageJ[k * n_step]), J, n2 * sizeof(iReg));
         std::memcpy(&(storageN[k * n_step]), mHat, n2 * sizeof(double));
         storageJsize[k] = n2;
      }
   }

   // Mean residual
   avg_resRelNorm /= nn_A;

   // Parallel construction of the SAM matrix outputs (iatN, jaN, coefN)
   iatN = new iReg[nn_A + 1];
   iatN[0] = 0;
   for (iReg i = 0; i < nn_A; ++i) {
      iatN[i + 1] = iatN[i] + storageJsize[i];
   }

   iReg total_nnz = iatN[nn_A];
   printf("nnzr = %4.6f, nnz = %d, avgResNorm = %f\n",static_cast<double>(total_nnz)/static_cast<double>(nn_A),total_nnz,avg_resRelNorm);
   jaN = new iReg[total_nnz];
   coefN = new double[total_nnz];

   // Parallel copy using column offsets
   #pragma omp parallel for num_threads(nthread) schedule(static)
   for (iReg i = 0; i < nn_A; ++i) {
      iReg nEnt = storageJsize[i];
      iReg dstOffset = iatN[i];
      iReg srcOffset = i * n_step;

      std::memcpy(&jaN[dstOffset], &storageJ[srcOffset], nEnt * sizeof(iReg));
      std::memcpy(&coefN[dstOffset], &storageN[srcOffset], nEnt * sizeof(double));
   }
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
void print_vector(const char* desc, iReg n, iReg* vec) {
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
void print_spVec(const char* desc, iReg n, int* vec) {
    printf("\n--- %s (Nonzero Elements) ---\n", desc);
    for (iReg i = 0; i < n; i++) {
        if (vec[i] != 0.0) {
            printf("Id [%ld]: %d\n", (long)i, vec[i]);
        }
    }
}
