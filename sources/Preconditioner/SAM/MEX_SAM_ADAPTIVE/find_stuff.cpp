#include "find_stuff.h"
#include <algorithm>
#include <vector>



bool mattia = false;
void fullA0k(ptrdiff_t nn_A, ptrdiff_t *iat0, ptrdiff_t *ja0, double *coef0, ptrdiff_t k, double *A0k){
   // Zero all the vector then focus on finding the nonzeros
   std::fill_n(A0k, nn_A, 0.0);

   // Loop over the rows
   for (ptrdiff_t row = 0; row < nn_A; ++row) {
      ptrdiff_t row_start = iat0[row];
      ptrdiff_t row_end = iat0[row+1];

      // Loop over the single row entries
      for (ptrdiff_t i = row_start; i < row_end; ++i) {
      	// If the column is the same of the current column k
         // Get the nonzero index
         if (ja0[i] == k) {
            A0k[row] = coef0[i];
            break;
         }
      }
   }
   return;
}

void findNonZeroInColJ(ptrdiff_t *J, ptrdiff_t *iatk, ptrdiff_t *jak, ptrdiff_t n2, ptrdiff_t *I, ptrdiff_t &sizeI){
   // Assume the pattern is symmetric
   // the nonzero row entries in column J correspond to
   // the nonzero column entries in row J

   const ptrdiff_t initial_count = sizeI;
   ptrdiff_t start, end;

   // First time entering, copy all
   if (initial_count == 0) {
      // Direct copy without any duplicate checks
      for (ptrdiff_t i = 0; i < n2; ++i) {
         start = iatk[J[i]];
         end = iatk[J[i] + 1];
         for (ptrdiff_t j = start; j < end; ++j) {
            I[sizeI] = jak[j];
            sizeI++;
         }
      }
   }
   else {
      bool duplicate;
      ptrdiff_t val;

      // Cycle over J to select the rows to get, need to check for no repeated indices
      for (ptrdiff_t i = 0; i < n2; ++i) {
         start = iatk[J[i]];
         end = iatk[J[i] + 1];
         for (ptrdiff_t j = start; j < end; ++j) {
            // Initialize the values for this jak
            val = jak[j];
            duplicate = false;

            // Check for duplicate value
            for (ptrdiff_t q = 0; q < sizeI; ++q) {
               if (I[q] == val) {
                  duplicate = true;
                  break;
               }
            }

            // This index is not repeated, copy it
            if (!duplicate) {
               I[sizeI++] = val;
            }
         }
      }
   }
}

void getA0k(double *a0k, ptrdiff_t *I, ptrdiff_t sizeI, ptrdiff_t oldSizeI, ptrdiff_t *iat0, ptrdiff_t *ja0, double *coef0, ptrdiff_t k){
   ptrdiff_t row, row_start, row_end;

   // Cycle over I to get which columns to seach for
   // Cycle over only the new entries of I
   for (ptrdiff_t i = oldSizeI; i < sizeI; ++i) {
      row = I[i];
      row_start = iat0[row];
      row_end = iat0[row + 1];

      // Set default value to 0.0
      a0k[i] = 0.0;

      // Loop over the single row entries
      for (ptrdiff_t j = row_start; j < row_end; ++j) {
      	// If entered here there is no column entry equal to k 
      	// set to zero
         if (ja0[j] == k) {
            a0k[i] = coef0[j];
            break; 
         }
      }
   }
   return;
}


void getAhat(ptrdiff_t *I, ptrdiff_t sizeI, ptrdiff_t *J, ptrdiff_t Jstart, ptrdiff_t Jend,
             ptrdiff_t *iatk, ptrdiff_t *jak, double *coefk, double *Ahat, ptrdiff_t &Astart) {

   if(mattia == true){
      // Cycle over the columns in J that have been added
      for (ptrdiff_t j = Jstart; j < Jend; ++j){
         // Cycle over the rows in I
         for (ptrdiff_t i = 0; i < sizeI; ++i){
            // Cycle over the chosen row
            for(ptrdiff_t k = iatk[I[i]]; k < iatk[I[i]+1]; ++k){
               // If the column in the row coincides with the column added then get the nonzero value
               // printf("%d %d %d ", k, jak[k], J[j]);
               if(jak[k] == J[j]){
                  Ahat[Astart] = coefk[k];
                  // printf("%f\n", coefk[k]);
                  Astart++;
                  break;
               }

               if(k == iatk[I[i] + 1] - 1){
                  // If entered here there is no column entry equal to k
                  // set to zero
                  // printf("0\n");
                  Ahat[Astart] = 0;
                  Astart++;
               }
               // else{
               //    // printf("\n");
               // }

            }
         }
      }
      return;
   }else{
      ptrdiff_t num_cols = Jend - Jstart;
      ptrdiff_t row, row_start, row_end, row_len;
      ptrdiff_t target_col, dest_idx, k;
      ptrdiff_t *row_cols, *it;

      // Loop over rows
      for (ptrdiff_t i = 0; i < sizeI; ++i) {
         row = I[i];
         row_start = iatk[row];
         row_end = iatk[row + 1];
         row_len = row_end - row_start;

         // Loop over the added columns
         for (ptrdiff_t j = Jstart; j < Jend; ++j) {
            target_col = J[j];

            // Calculate the exact 1D destination index in Ahat
            dest_idx = Astart + (j - Jstart) * sizeI + i;

            // Binary search assuming 'jak' is sorted per row
            row_cols = &jak[row_start];
            it = std::lower_bound(row_cols, row_cols + row_len, target_col);

            // Assign value if found, otherwise explicitly assign 0 (fixes the empty row bug)
            if (it != row_cols + row_len && *it == target_col) {
               k = row_start + (it - row_cols);
               Ahat[dest_idx] = coefk[k];
            } else {
               Ahat[dest_idx] = 0.0;
            }
         }
      }

      // Update Astart once at the very end
      Astart += num_cols * sizeI;
      return;
   }
}


// AJ is saved rowwise then used as transposed in the blas gemm to have better memory access in creating it
void getAJ(ptrdiff_t *J, ptrdiff_t Jsize, ptrdiff_t nn_A, ptrdiff_t *iatk, ptrdiff_t *jak, double *coefk, double *AJ){
   ptrdiff_t Astart = 0;

   // Cycle over the rows in I
   for (ptrdiff_t i = 0; i < nn_A; ++i){
      // currJ = 0;

      // Loop over all possible J
      for (ptrdiff_t currJ = 0; currJ < Jsize; ++currJ){
      	// Cycle over the chosen row
      	for(ptrdiff_t k = iatk[i]; k < iatk[i+1]; ++k){
      		// If column index is equal to the chosen J column then add it
      		if(jak[k] == J[currJ]){
      			// printf("A[%d] = %f\n", Astart,coefk[k]);
      			AJ[Astart] = coefk[k];
         		Astart++;
         		break;
      		}

      		if(k == iatk[i+1] - 1){
				   // If entered here there is no column entry equal to this J
				   // set to zero
				   // printf("A[%d] = %f\n", Astart,0);
				   AJ[Astart] = 0;
				   Astart++;
				   break;
				}
      	}
      }
   }
   return;
}


void fillL(ptrdiff_t *L, double *res, ptrdiff_t nn_A, ptrdiff_t &usedL){
   usedL = 0;
   // Loop over all residual entries
   for (ptrdiff_t i = 0; i < nn_A; ++i){
      // The residual is not numerically zero
      // Add it as a possible column to be computed
      if (std::abs(res[i]) > 1e-13){
         L[usedL] = i;
         usedL++;
      }
   }
   return;
}

void findJtilde(ptrdiff_t *Jtilde, ptrdiff_t &JtildeSize, ptrdiff_t *L, ptrdiff_t sizeL, ptrdiff_t *iatk, ptrdiff_t *jak, ptrdiff_t *J, ptrdiff_t sizeJ){

   // Sanity check for sizeL == 0
   if (sizeL == 0){
      printf("sizeL == 0, check\n");
      return;
   }

   // Initialize counter to 0
   JtildeSize = 0;

   // Flag for fast exit in case of repeated index
   bool skip = false;
   ptrdiff_t maxJsize;
   for (ptrdiff_t i = 0; i < sizeL; ++i){
      for (ptrdiff_t j = iatk[L[i]]; j < iatk[L[i] + 1]; ++j){
         // Get the max size to search for duplicates
         maxJsize = std::max(sizeJ,JtildeSize);
         // std::cout << "maxJsize = " << maxJsize << std::endl;
         // std::cout << "ja["<< j << "] = " << jak[j] << std::endl;
         // Check for duplicate value
         for (ptrdiff_t q = 0; q < maxJsize; ++q){
            // Check if it is already in Jtilde
            if(q < JtildeSize){
               if(jak[j] == Jtilde[q]){
                  skip = true;
                  break;
               }
            }

            // Check if it is already in J
            if (q < sizeJ){
               if(jak[j] == J[q]){
                  skip = true;
                  break;
               }
            }
         }
   
         // This index is repeated, fast exit
         if (skip){
            // Reset flag to false
            skip = false;
         }
         else{
            // Index is not repeated, add it
            Jtilde[JtildeSize] = jak[j];
            JtildeSize++;
         }
      }
   }
   return;
}

// Compute A(:,Jtilde) in colmajor
void fullAJtilde(ptrdiff_t nn_A, ptrdiff_t *iatk, ptrdiff_t *jak, double *coefk, ptrdiff_t *Jtilde, ptrdiff_t JtildeSize, double *AJtilde){
   
   std::fill_n(AJtilde, nn_A * JtildeSize, 0.0);

   // Find both the minimum and maximum column indices in Jtilde
   auto minmax = std::minmax_element(Jtilde, Jtilde + JtildeSize);
   ptrdiff_t min_col = *minmax.first;
   ptrdiff_t max_col = *minmax.second;
   ptrdiff_t range = max_col - min_col + 1;

   // Fill the lookup array based purely on the range span
   std::vector<ptrdiff_t> lookup(range, -1);
   for (ptrdiff_t c = 0; c < JtildeSize; ++c) {
       lookup[Jtilde[c] - min_col] = c; // Apply the negative offset
   }

   // Iterate through the CSR matrix
   for (ptrdiff_t r = 0; r < nn_A; ++r) {
       ptrdiff_t row_end = iatk[r + 1];
       for (ptrdiff_t k = iatk[r]; k < row_end; ++k) {
           ptrdiff_t col = jak[k];
           
           // Quick bounds check using the min/max cluster boundaries
           if (col >= min_col && col <= max_col) {
               ptrdiff_t c = lookup[col - min_col]; // Apply same offset to query
               
               if (c != -1) {
                   // Compute column-major index: row + (col_index * total_rows)
                   AJtilde[r + c * nn_A] = coefk[k];
               }
           }
       }
   }
   // // Loop over possible columns for Jtilde
   // for (ptrdiff_t j = 0; j < JtildeSize; ++j){
   //    // Loop over the rows
   //    for (ptrdiff_t row = 0; row < nn_A; ++row){
   //       // Loop over the single row entries
   //       for (ptrdiff_t i = iatk[row]; i < iatk[row+1]; ++i){
         
   //          // If the column is the same of the current column k
   //          // Get the nonzero index
   //          if (jak[i] == Jtilde[j]){
   //             AJtilde[row + j*nn_A] = coefk[i];
   //             break;
   //          }
   
   //          // No entry was found, set to zero
   //          if (i == iatk[i+1] - 1){
   //             AJtilde[row + j*nn_A] = 0;
   //          }
   //       }
   //    }
   // }
   return;
}
