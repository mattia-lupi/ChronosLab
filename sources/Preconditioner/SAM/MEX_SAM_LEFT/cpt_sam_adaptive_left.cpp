// #pragma once
#include "cpt_sam_adaptive_left.h"
#include <vector>
#include <iostream>
#include "lapack.h"
#include "cblas.h"

// LAPACK Fortran routines declarations
// extern "C" {
//    void dgeqrf_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info);
//    void dormqr_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* lwork, int* info);
//    void dtrtrs_(char* uplo, char* trans, char* diag, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);
// }

void fullA0k(int nn_A, int *iat0, int *ja0, double *coef0, int k, double *A0k);
void findNonZeroInColJ(int *J, int *iatk, int *jak, int n2, int *I, int &sizeI);
void getA0k(double *a0k,  int *I, int sizeI, int oldSizeI, int *iat0,  int *ja0, double *coef0,  int k);
void getAhat( int *I,  int sizeI,  int *J,  int Jstart,  int Jend,
             int *iatk,  int *jak, double *coefk, double *Ahat,  int &Astart);
void computeFirstQR(double *Ahat, int sizeI, int sizeJ, double *R, double *tau, double *work, int lwork, int &info);
void computeUpdateQR();
void applyFirstQt(double *Ahat, int sizeI, int sizeJ, double *tau, double *a0k, double *work, int lwork, int &info);
void applyR(int sizeJ, double *R, double *a0k, int &info);
void getAJ(int *J, int Jsize, int nn_A, int *iatk, int *jak, double *coefk, double *AJ);
void cptRes(int nn_A, int sizeJ, double *A0k, double *AJ, double *mHat, double &resRelNorm);
void fillL(double *L, double *res, int nn_A, int &usedL);


void print_matrix(const char* desc, int m, int n, double* mat, int lda) {
    printf("\n--- %s (%dx%d) ---\n", desc, m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.4f ", mat[i + j * lda]);
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
                     int *iatN, int *jaN,double *coefN){

   // Suppose the pattern is symmetric
   // Count the maximum number of entries per row
   int max_nnz_r = 0;
   for (int i = 0; i < nn_A; ++i){
      max_nnz_r = std::max(max_nnz_r,iatk[i+1]-iatk[i]);
   }

    int maxJsize = n_step*step_size;

   // Cycle over the columns
   for ( int k = 0; k < nn_A; ++k){
      // All allocation could be done outside the cycle multiplying the size of the memory
      // needed by the size of the parallel thread pool. Then each thread can read and write
      // only on its part of the allocated memory. in this way the number of allocations is 
      // reduced to 1 per object instead of being nn_A while also reducing the actual max size
      // allocated

      double resRelNorm = 1.0;
      int usedL;

      // Allocate J
      std::vector< int> Jvec(maxJsize);
      int *J = Jvec.data();
      // Set first sparsity pattern to be diagonal
      J[0] = k;

      // Allocate Jtilde
      std::vector< int> Jtildevec(nn_A);
      int *Jtilde = Jtildevec.data();

      // Allocate I
      std::vector< int> Ivec(max_nnz_r);
      int *I = Ivec.data();

      // Allocate space for A(:,J)
      std::vector<double> aJvec(maxJsize);
      double *aJ = aJvec.data();

      // Allocate space for A0(I,k)
      std::vector<double> a0kvec(max_nnz_r);
      double *a0k = a0kvec.data();

      // Allocate space for A0(:,k)
      std::vector<double> A0kvec(nn_A);
      double *A0k = A0kvec.data();

      // Allocate space for the residuals
      std::vector<double> Lvec(nn_A);
      double *L = Lvec.data();

      // Fill the current column of the "old" matrix 
      // Do it once per column
      fullA0k(nn_A, iat0, ja0, coef0, k, A0k);
      // Compute the norm only once
      double normA0k = cblas_dnrm2(nn_A,A0k,1);

      // Allocate Ahat buffer for the max possible size
      // Contains the householder vectors for all the updates
      std::vector<double> AhatBuffer(max_nnz_r*maxJsize);
      double *Ahat = AhatBuffer.data();

      // Allocate AJ buffer for the max possible size
      // Contains the householder vectors for all the updates
      std::vector<double> AJBuffer(nn_A*maxJsize);
      double *AJ = AJBuffer.data();

      // Allocate tau for the max possible size
      std::vector<double> tauVec(maxJsize);
      double *tau = tauVec.data();

      // Allocate R for the max possible size
      std::vector<double> RVec((maxJsize)*(maxJsize+1)/2);
      double *R = RVec.data();
      int indexR = 0;

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
      int oldSizeI, sizeI = 0;
      int Astart = 0;

      for (int t = 0; t < n_step; ++t){

         oldSizeI = sizeI;

         // Get the row entries of the matrix in columns J
         findNonZeroInColJ(J, iatk, jak, n2, I, sizeI);

         // check what happens if I does not change
         print_vector("Vector I", sizeI, I);

         // Get the new piece of A0k
         getA0k(a0k, I, sizeI, oldSizeI, iat0, ja0, coef0, k);

         print_vector("Vector A0k", sizeI, a0k);

         // Add to Ahat the new part of the matrix on which to work
         getAhat(I, sizeI, J, n2old, n2, iatk, jak, coefk, Ahat, Astart);
         printf("Astart %d\n", Astart);

         print_matrix("Matrix Ahat", sizeI, n2, Ahat, sizeI);

         if (t == 0){
            // Compute the QR factorization of the first matrix
            computeFirstQR(Ahat, sizeI, n2, R, tau, work, lwork, info);
            print_matrix("Matrix Ahat QR", sizeI, n2, Ahat, sizeI);

            // Apply Qt for the first matrix
            applyFirstQt(Ahat, sizeI, n2, tau, a0k, work, lwork, info);
            print_vector("Vector chat", sizeI, a0k);
         } else{
            computeUpdateQR();
         }

         // Solve the triangolar system
         applyR(n2, R, a0k, info);
         print_vector("Vector mhat", sizeI, a0k);

         // Get the matrix A(:,J)
         getAJ(J, n2, nn_A, iatk, jak, coefk, AJ);
         print_matrix("Matrix Aj", nn_A, n2, AJ, nn_A);

         // Assign the norm of A0(:,k) to resRelNorm
         resRelNorm = normA0k;
         cptRes(nn_A, n2, A0k, AJ, a0k, resRelNorm);
         print_matrix("Matrix Aj", nn_A, n2, AJ, nn_A);

         if (resRelNorm < eps){
            break;
         }

         if (t < n_step){
            fillL(L, AJ, nn_A, usedL);
            print_vector("L", usedL, L);
         }
      }
      return;
   }

   return;
}

void fullA0k(int nn_A, int *iat0, int *ja0, double *coef0, int k, double *A0k){
   // Loop over the rows
   for (int row = 0; row < nn_A; ++row){
      // Loop over the single row entries
      for (int i = iat0[row]; i < iat0[row+1]; ++i){
         // If the column is the same of the current column k
         // Get the nonzero index
         if (ja0[i] == k){
            A0k[row] = coef0[i];
            break;
         }

         // No entry was found, set to zero
         if (i == iat0[i+1] - 1){
            A0k[row] = 0;
         }
      }
   }
   return;
}
void findNonZeroInColJ( int *J, int *iatk,  int *jak, int n2,  int *I, int &sizeI){

   // Assume the pattern is symmetric
   // the nonzero row entries in column J correspond to 
   // the nonzero column entries in row J

   int initial_count = sizeI;

   // Flag for fast exit in case of repeated index
   bool exitt = false;
   // Cycle over J to select the rows to get, need to check for no repeated indices
   for (int i = 0; i < n2; ++i){
      for (int j = iatk[J[i]]; j < iatk[J[i] + 1]; ++j){
         // If the sizeI is 0 then avoid checks and just copy all the first row
         if(initial_count != 0){
            // Check for duplicate value
            for (int q = 0; q < sizeI; ++q){
               if(jak[j] == I[q]){
                  exitt = true;
                  break;
               }
            }
   
            // This index is repeated, fast exit
            if (exitt){
               // Reset flag to false
               exitt = false;
               break;
            }
         }

         // This index is new
         // Copy this nonzero column entry inside I
         I[sizeI] = jak[j];
         sizeI++;
      }
   }
   return;
}


void getA0k(double *a0k,  int *I, int sizeI, int oldSizeI, int *iat0,  int *ja0, double *coef0,  int k){
   // Cycle over I to get which columns to seach for
   // Cycle over only the new entries of I
   for (int i = oldSizeI; i < sizeI; ++i){
      for (int j = iat0[I[i]]; j < iat0[I[i] + 1]; ++j){
         if (ja0[j] == k){
            // Get the entry of column k
            a0k[i] = coef0[j];
            break;
         }

         if(j == iat0[I[i] + 1] - 1){
            // If entered here there is no column entry equal to k 
            // set to zero
            a0k[i] = 0;
         }
      }
   }
   return;
}

void getAhat(int *I, int sizeI, int *J, int Jstart, int Jend,
             int *iatk, int *jak, double *coefk, double *Ahat, int &Astart){
   // Cycle over the columns in J that have been added
   for (int j = Jstart; j < Jend; ++j){
      // Cycle over the rows in I
      for (int i = 0; i < sizeI; ++i){
         // Cycle over the chosen row
         for(int k = iatk[I[i]]; k < iatk[I[i]+1]; ++k){
            // If the column in the row coincides with the column added then get the nonzero value
            // printf("%d %d %d ", k, jak[k], J[j]);
            if(jak[k] == J[j]){
               Ahat[Astart] = coefk[k];
               // printf("%f\n", coefk[k]);
               Astart++;
               break;
            }

            if(j == iatk[I[i] + 1] - 1){
               // If entered here there is no column entry equal to k 
               // set to zero
               // printf("0\n");
               Ahat[Astart] = 0;
               Astart++;
            }
            else{
               // printf("\n");
            }

         }
      }
   }
   return;
}

void computeFirstQR(double *Ahat, int sizeI, int sizeJ, double *R, double *tau, double *work, int lwork, int &info){
   // Compute QR Factorization
   dgeqrf_(&sizeI, &sizeJ, Ahat, &sizeI, tau, work, &lwork, &info);
   if (info != 0){
      printf("Exit at first QR due to error %d\n", info);
      return;
   }

   // Copy R into a separate vector to store and use also the updated R later
   R[0] = Ahat[0];

   return;
}

void applyFirstQt(double *Ahat, int sizeI, int sizeJ, double *tau, double *a0k, double *work, int lwork, int &info){
   char side = 'L';
   char trans = 'T';

   // Check why I needed to add the fortran string length sizes
   dormqr_(&side, &trans, &sizeI, &sizeJ, &sizeJ, Ahat, &sizeI, tau, a0k, &sizeI, work, &lwork, &info,3,3);
   if (info != 0){
      printf("Exit at first Qt apply due to error %d\n", info);
      return;
   }

   return;
}

void applyR(int sizeJ, double *R, double *a0k, int &info){
   char uplo = 'U';
   char trans = 'N';
   char diag = 'N';
   int nrhs = 1;

   // Check why I needed to add the fortran string length sizes
   dtrtrs_(&uplo, &trans, &diag, &sizeJ, &nrhs, R, &sizeJ, a0k, &sizeJ, &info,3,3,3);
   if (info != 0){
      printf("Exit at R apply due to error %d\n", info);
      return;
   }

   return;
}


// AJ is saved rowwise then used as transposed in the blas gemm to have better memory access in creating it
void getAJ(int *J, int Jsize, int nn_A, int *iatk, int *jak, double *coefk, double *AJ){
   int Astart = 0;
   int currJ = 0;

   // printf("%d\n",Jsize);

   // Cycle over the rows in I
   for (int i = 0; i < nn_A; ++i){
      currJ = 0;

      // Cycle over the chosen row
      for(int k = iatk[i]; k < iatk[i+1]; ++k){
         // printf("k %d, jak %d, Jcurr %d\n",k,jak[k],J[currJ]);
         if(jak[k] > J[currJ]){
            // If entered here there is no column entry equal to k 
            // set to zero
            // printf("0\n");
            // printf("entered when k %d, jak %d, Jcurr %d\n",k,jak[k],J[currJ]);
            AJ[Astart] = 0;
            // printf("A[%d] = %f\n", Astart,0);
            Astart++;
            currJ++;
         }

         // Check if reached the max size for J
         if(currJ >= Jsize){
            // printf("entered break at k = %d\n",k);
            break;
         }

         // If the column in the row coincides with the column needed then add
         if(jak[k] == J[currJ]){
            AJ[Astart] = coefk[k];
            // printf("A[%d] = %f\n", Astart,coefk[k]);
            Astart++;
            currJ++;

         }

         // Check if reached the max size for J
         if(currJ >= Jsize){
            // printf("entered break at k = %d\n",k);
            break;
         }
      }
   }
   return;
}

void cptRes(int nn_A, int sizeJ, double *A0k, double *AJ, double *mHat, double &resRelNorm){
   // Compute the product AJ*mHat and save it in AJ
   cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nn_A, 1, sizeJ, 1.0, AJ, sizeJ, mHat, 1, 0.0, AJ, sizeJ);

   // Compute the norm of Aj*mHat
   double normAjMh = cblas_dnrm2(nn_A,AJ,1);

   // Compute res = Aj*mHat - A0k and save it in AJ
   cblas_daxpy(nn_A, -1.0, A0k, 1, AJ, 1);

   // Compute the residual norm
   resRelNorm = 2*cblas_dnrm2(nn_A,AJ,1)/(normAjMh+resRelNorm);
   // printf("resNorm %f\n", resRelNorm);

   return;
}

void fillL(double *L, double *res, int nn_A, int &usedL){
   usedL = 0;
   // Loop over all residual entries
   for (int i = 0; i < nn_A; ++i){
      // The residual is not numerically zero
      // Add it as a possible column to be computed
      if (std::abs(res[i]) > 1e-13){
         L[usedL] = i;
         usedL++;
      }
   }
   return;
}

void findJtilde(){
   // Assume the pattern is symmetric
   // the nonzero row entries in column J correspond to 
   // the nonzero column entries in row J

   int initial_count = sizeI;

   // Flag for fast exit in case of repeated index
   bool exitt = false;
   // Cycle over J to select the rows to get, need to check for no repeated indices
   for (int i = 0; i < n2; ++i){
      for (int j = iatk[J[i]]; j < iatk[J[i] + 1]; ++j){
         // If the sizeI is 0 then avoid checks and just copy all the first row
         if(initial_count != 0){
            // Check for duplicate value
            for (int q = 0; q < sizeI; ++q){
               if(jak[j] == I[q]){
                  exitt = true;
                  break;
               }
            }
   
            // This index is repeated, fast exit
            if (exitt){
               // Reset flag to false
               exitt = false;
               break;
            }
         }

         // This index is new
         // Copy this nonzero column entry inside I
         I[sizeI] = jak[j];
         sizeI++;
      }
   }
   return;
}

void computeUpdateQR(){

   return;
}




