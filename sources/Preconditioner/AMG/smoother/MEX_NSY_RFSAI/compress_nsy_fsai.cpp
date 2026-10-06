#include "compress_nsy_fsai.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <omp.h>

using namespace std;

int compress_nsy_fsai(const int nn, int &nt_FL, int &nt_FU,
                      int *&iat_FL, int *&ja_FL, int *&iat_FU, int *&ja_FU,
                      double *&coef_FL, double *&coef_FU,
                      const int num_threads){

   // Allocate new arrays
   int *iat_FL_new = (int*) malloc((nn + 1) * sizeof(int));
   if (iat_FL_new == nullptr){
      cout << "Allocation Error in compress_nsy_fsai" << endl;
      return 1;
   }
   iat_FU = (int*) malloc((nn + 1) * sizeof(int));
   if (iat_FU == nullptr){
      cout << "Allocation Error in compress_nsy_fsai" << endl;
      free(iat_FL_new);
      return 1;
   }

   // Count true non-zero entries
   #pragma omp parallel for schedule(static) num_threads(num_threads)
   for (int i = 0; i < nn; i++){
      int count = 0;
      int istart = iat_FL[i];
      int iend = iat_FL[i+1];
      for (int j = istart; j < iend; j++){
         if (coef_FL[j] != 0.0) count++;
      }
      iat_FL_new[i+1] = count;
   }

   iat_FL_new[0] = 0;
   for (int i = 0; i < nn; i++){
      iat_FL_new[i+1] += iat_FL_new[i];
   }
   nt_FL = iat_FL_new[nn];

   // Allocate new arrays
   int *ja_FL_new = (int*) malloc(nt_FL * sizeof(int));
   double *coef_FL_new = (double*) malloc(nt_FL * sizeof(double));
   if (ja_FL_new == nullptr || coef_FL_new == nullptr){
      cout << "Allocation Error in compress_nsy_fsai" << endl;
      return 1;
   }

   // Pass 2: Extract non-zeros of FL in parallel
   #pragma omp parallel for schedule(guided) num_threads(num_threads)
   for (int i = 0; i < nn; i++){
      int out_idx = iat_FL_new[i];
      int istart = iat_FL[i];
      int iend = iat_FL[i+1];
      for (int j = istart; j < iend; j++){
         if (coef_FL[j] != 0.0){
            ja_FL_new[out_idx] = ja_FL[j];
            coef_FL_new[out_idx] = coef_FL[j];
            out_idx++;
         }
      }
   }

   // --- 2. Transpose FUT into FU in Compressed Format ---
   std::vector<int> thread_counts((size_t)num_threads * nn, 0);
   std::vector<int> thread_offsets((size_t)num_threads * nn, 0);
   std::vector<int> col_totals(nn, 0);

   #pragma omp parallel num_threads(num_threads)
   {
      int tid = omp_get_thread_num();
      int actual_nth = omp_get_num_threads();

      int r_start = (int)(((long long)tid * nn) / actual_nth);
      int r_end   = (int)(((long long)(tid + 1) * nn) / actual_nth);

      int *my_counts = &thread_counts[(size_t)tid * nn];

      // Pass 1: Thread-local column histogram of FUT
      for (int i = r_start; i < r_end; ++i) {
         int row_end = iat_FL[i+1];
         for (int j = iat_FL[i]; j < row_end; ++j) {
            if (coef_FU[j] != 0.0) {
               my_counts[ja_FL[j]]++;
            }
         }
      }

      #pragma omp barrier

      // Pass 2: Prefix sum across threads and columns
      #pragma omp for schedule(static)
      for (int c = 0; c < nn; ++c) {
         int running_sum = 0;
         for (int t = 0; t < actual_nth; ++t) {
            thread_offsets[(size_t)t * nn + c] = running_sum;
            running_sum += thread_counts[(size_t)t * nn + c];
         }
         col_totals[c] = running_sum;
      }
   }

   // Build iat_FU
   iat_FU[0] = 0;
   for (int c = 0; c < nn; ++c) {
      iat_FU[c + 1] = iat_FU[c] + col_totals[c];
   }
   nt_FU = iat_FU[nn];

   // Allocate FU arrays
   ja_FU = (int*) malloc(nt_FU * sizeof(int));
   double *coef_FU_new = (double*) malloc(nt_FU * sizeof(double));
   if (!ja_FU || !coef_FU_new) return 1;

   #pragma omp parallel num_threads(num_threads)
   {
      int tid = omp_get_thread_num();
      int actual_nth = omp_get_num_threads();

      int r_start = (int)(((long long)tid * nn) / actual_nth);
      int r_end   = (int)(((long long)(tid + 1) * nn) / actual_nth);

      int *my_offsets = &thread_offsets[(size_t)tid * nn];

      for (int c = 0; c < nn; ++c) {
         my_offsets[c] += iat_FU[c];
      }

      // Pass 3: Scatter transposed entries
      for (int i = r_start; i < r_end; ++i) {
         int row_end = iat_FL[i+1];
         for (int j = iat_FL[i]; j < row_end; ++j) {
            if (coef_FU[j] != 0.0) {
               int col = ja_FL[j];
               int dest = my_offsets[col]++;
               ja_FU[dest] = i;
               coef_FU_new[dest] = coef_FU[j];
            }
         }
      }
   }

   // Swap pointers for the entries of FU
   free(coef_FU);
   coef_FU = coef_FU_new;

   // Update FL pointers
   free(coef_FL);
   coef_FL = coef_FL_new;
   free(iat_FL);
   free(ja_FL);
   iat_FL = iat_FL_new;
   ja_FL = ja_FL_new;

   return 0;

}
