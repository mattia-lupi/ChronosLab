#include "transpose.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <omp.h>

int transpose(const int nrows, const int ncols, const int *const iat, const int *const ja,
              const double *const coef, int *&iat_T, int *&ja_T, double *&coef_T,
              const int num_threads){

   // Allocate output and scracth
   iat_T = (int*) malloc((ncols+1) * sizeof(int));
   int nterm = iat[nrows];
   ja_T = (int*) malloc((nterm) * sizeof(int));
   coef_T = (double*) malloc((nterm) * sizeof(double));
   if (iat_T == nullptr || ja_T == nullptr || coef_T == nullptr){
      // Allocation error
      std::cout << "Allocation Error in transpose" << std::endl;
      return 1;
   }

   if (num_threads == 1) {
      // Clean sequential transpose
      int *ISCR = (int*) malloc((ncols + 1) * sizeof(int));
      if (!ISCR) return 1;

      // Initialize pointers
      std::fill_n(iat_T, ncols + 1, 0);

      // Count non-zeroes for each column of the input matrix
      for (int i = 0; i < nrows; i++){
         for (int j = iat[i]; j < iat[i+1]; j++) iat_T[ja[j]]++;
      }

      // Set pointers
      ISCR[0] = 0;
      for (int i = 1; i <= ncols; i++) ISCR[i] = ISCR[i-1] + iat_T[i-1];
      for (int i = 0; i <= ncols; i++) iat_T[i] = ISCR[i];

      // Transpose column indices and coefficients
      for (int i = 0; i < nrows; i++){
         for (int j = iat[i]; j < iat[i+1]; j++){
            int ind = ISCR[ja[j]]++;
            ja_T[ind] = i;
            coef_T[ind] = coef[j];
         }
      }

      // Deallocate scratch
      free(ISCR);
      return 0;
   }

   // Parallel multi-thread transpose using num_threads threads
   std::vector<int> thread_counts((size_t)num_threads * ncols, 0);
   std::vector<int> thread_offsets((size_t)num_threads * ncols, 0);
   std::vector<int> col_totals(ncols, 0);

   #pragma omp parallel num_threads(num_threads)
   {
      int tid = omp_get_thread_num();
      int actual_nth = omp_get_num_threads();

      int r_start = (int)(((long long)tid * nrows) / actual_nth);
      int r_end   = (int)(((long long)(tid + 1) * nrows) / actual_nth);

      int *my_counts = &thread_counts[(size_t)tid * ncols];

      // Pass 1: Thread-local column histogram
      for (int i = r_start; i < r_end; ++i) {
         int row_end = iat[i+1];
         for (int j = iat[i]; j < row_end; ++j) {
            my_counts[ja[j]]++;
         }
      }

      #pragma omp barrier

      // Pass 2: Prefix sum across threads and columns
      #pragma omp for schedule(static)
      for (int c = 0; c < ncols; ++c) {
         int running_sum = 0;
         for (int t = 0; t < actual_nth; ++t) {
            thread_offsets[(size_t)t * ncols + c] = running_sum;
            running_sum += thread_counts[(size_t)t * ncols + c];
         }
         col_totals[c] = running_sum;
      }
   }

   // Build iat_T
   iat_T[0] = 0;
   for (int c = 0; c < ncols; ++c) {
      iat_T[c + 1] = iat_T[c] + col_totals[c];
   }

   #pragma omp parallel num_threads(num_threads)
   {
      int tid = omp_get_thread_num();
      int actual_nth = omp_get_num_threads();

      int r_start = (int)(((long long)tid * nrows) / actual_nth);
      int r_end   = (int)(((long long)(tid + 1) * nrows) / actual_nth);

      int *my_offsets = &thread_offsets[(size_t)tid * ncols];

      for (int c = 0; c < ncols; ++c) {
         my_offsets[c] += iat_T[c];
      }

      // Pass 3: Scatter transposed entries
      for (int i = r_start; i < r_end; ++i) {
         int row_end = iat[i+1];
         for (int j = iat[i]; j < row_end; ++j) {
            int col = ja[j];
            int dest = my_offsets[col]++;
            ja_T[dest] = i;
            coef_T[dest] = coef[j];
         }
      }
   }

   return 0;
}
