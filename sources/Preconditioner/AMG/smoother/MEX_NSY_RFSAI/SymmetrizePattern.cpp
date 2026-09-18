#include "SymmetrizePattern.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include "transpose.h"
#include "merge_row_patt.h"

// Symmetrizes the pattern of a non-symmetric matrix by padding with zeroes missing enries
int SymmetrizePattern(const int nrows, int *& iat, int *& ja, double *& coef,
                      double *&coef_T, const int num_threads){

   // Transpose input pattern
   int *iat_T = nullptr;
   int *ja_T = nullptr;
   int ierr = transpose(nrows,nrows,iat,ja,coef,iat_T,ja_T,coef_T,num_threads);
   if (ierr != 0) return 1;

   int *iat_tmp = (int*) malloc((nrows + 1) * sizeof(int));
   if (iat_tmp == nullptr) {
      free(iat_T); free(ja_T); free(coef_T);
      return 1;
   }

   // Pass 1: Parallel symbolic merge counting
   #pragma omp parallel for schedule(static) num_threads(num_threads)
   for (int i = 0; i < nrows; i++){
      int istart_A = iat[i];
      int len_A = iat[i+1] - istart_A;
      int istart_T = iat_T[i];
      int len_T = iat_T[i+1] - istart_T;

      iat_tmp[i+1] = count_merged_row_patt(len_A, &ja[istart_A], len_T, &ja_T[istart_T]);
   }

   // Prefix sum to compute row pointers
   iat_tmp[0] = 0;
   for (int i = 0; i < nrows; i++){
      iat_tmp[i+1] += iat_tmp[i];
   }

   int total_nnz = iat_tmp[nrows];
   int *ja_tmp = (int*) malloc(total_nnz * sizeof(int));
   double *coef_tmp = (double*) malloc(total_nnz * sizeof(double));

   if (ja_tmp == nullptr || coef_tmp == nullptr) {
      free(iat_T); free(ja_T); free(coef_T);
      free(iat_tmp);
      if (ja_tmp) free(ja_tmp);
      if (coef_tmp) free(coef_tmp);
      return 1;
   }

   // Pass 2: Parallel numeric merge directly into pre-allocated memory
   #pragma omp parallel for schedule(static) num_threads(num_threads)
   for (int i = 0; i < nrows; i++){
      int istart_A = iat[i];
      int len_A = iat[i+1] - istart_A;
      int istart_T = iat_T[i];
      int len_T = iat_T[i+1] - istart_T;

      int len_out = 0;
      int out_offset = iat_tmp[i];

      merge_row_patt(len_A, &ja[istart_A], &coef[istart_A],
                     len_T, &ja_T[istart_T], &coef_T[istart_T],
                     len_out, &ja_tmp[out_offset], &coef_tmp[out_offset]);
   }

   // Free temporary transposed matrix and previous arrays
   free(iat_T);
   free(ja_T);
   free(coef_T);

   // Free original matrix storage
   free(iat);
   free(ja);
   free(coef);

   // Replace with merged matrix
   iat = iat_tmp;
   ja = ja_tmp;
   coef = coef_tmp;

   // Transpose again to compute coef_T
   ierr = transpose(nrows, nrows, iat, ja, coef, iat_T, ja_T, coef_T, num_threads);
   if (ierr != 0) return 1;

   // Free unused topology of transposed matrix
   free(iat_T);
   free(ja_T);

   return 0;
}
