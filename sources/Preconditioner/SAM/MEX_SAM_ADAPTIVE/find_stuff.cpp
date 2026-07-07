#include "find_stuff.h"
#include <algorithm>
#include <vector>
#include <iostream>

bool mattia = false;
void fullA0k(iExt nn_A, iExt *iat0, iReg *ja0, double *coef0, iExt k, double *A0k){
   // Zero all the vector then focus on finding the nonzeros
   std::fill_n(A0k, nn_A, 0.0);

   // Loop over the rows
   for (iReg row = 0; row < nn_A; ++row) {
      iReg row_start = iat0[row];
      iReg row_end = iat0[row+1];

      // Loop over the single row entries
      for (iReg i = row_start; i < row_end; ++i) {
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

void findNonZeroInColJ(iReg *J, iExt *iatk, iReg *jak, iReg n2, iReg *I, iReg &sizeI){
   // Assume the pattern is symmetric
   // the nonzero row entries in column J correspond to
   // the nonzero column entries in row J

   const iReg initial_count = sizeI;
   iReg start, end;

   // First time entering, copy all
   if (initial_count == 0) {
      // Direct copy without any duplicate checks
      for (iReg i = 0; i < n2; ++i) {
         start = iatk[J[i]];
         end = iatk[J[i] + 1];
         for (iReg j = start; j < end; ++j) {
            I[sizeI] = jak[j];
            sizeI++;
         }
      }
   }
   else {
      bool duplicate;
      iReg val;

      // Cycle over J to select the rows to get, need to check for no repeated indices
      for (iReg i = 0; i < n2; ++i) {
         start = iatk[J[i]];
         end = iatk[J[i] + 1];
         for (iReg j = start; j < end; ++j) {
            // Initialize the values for this jak
            val = jak[j];
            duplicate = false;

            // Check for duplicate value
            for (iReg q = 0; q < sizeI; ++q) {
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

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, iExt *iat0, iReg *ja0, double *coef0, iExt k){
   iReg row, row_start, row_end;

   // Cycle over I to get which columns to seach for
   // Cycle over only the new entries of I
   for (iReg i = oldSizeI; i < sizeI; ++i) {
      row = I[i];
      row_start = iat0[row];
      row_end = iat0[row + 1];

      // Set default value to 0.0
      a0k[i] = 0.0;

      // Loop over the single row entries
      for (iReg j = row_start; j < row_end; ++j) {
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


void getAhat(iReg *I, iReg sizeI, iReg *J, iReg Jstart, iReg Jend,
             iExt *iatk, iReg *jak, double *coefk, double *Ahat, iReg &Astart) {

   if(mattia == true){
      // Cycle over the columns in J that have been added
      for (iReg j = Jstart; j < Jend; ++j){
         // Cycle over the rows in I
         for (iReg i = 0; i < sizeI; ++i){
            // Cycle over the chosen row
            for(iExt k = iatk[I[i]]; k < iatk[I[i]+1]; ++k){
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
      iReg num_cols = Jend - Jstart;
      iReg row, row_start, row_end, row_len;
      iReg target_col, dest_idx, k;
      iReg *row_cols, *it;

      // Loop over rows
      for (iReg i = 0; i < sizeI; ++i) {
         row = I[i];
         row_start = iatk[row];
         row_end = iatk[row + 1];
         row_len = row_end - row_start;

         // Loop over the added columns
         for (iReg j = Jstart; j < Jend; ++j) {
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
void getAJ(iReg *J, iReg Jsize, iExt nn_A, iExt *iatk, iReg *jak, double *coefk, double *AJ){
   iReg Astart = 0;

   // Cycle over the rows in I
   for (iReg i = 0; i < nn_A; ++i){
      // currJ = 0;

      // Loop over all possible J
      for (iReg currJ = 0; currJ < Jsize; ++currJ){
      	// Cycle over the chosen row
      	for(iExt k = iatk[i]; k < iatk[i+1]; ++k){
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


void fillL(iReg *L, double *res, iExt nn_A, iReg &usedL){
   usedL = 0;
   // Loop over all residual entries
   for (iReg i = 0; i < nn_A; ++i){
      // The residual is not numerically zero
      // Add it as a possible column to be computed
      if (std::abs(res[i]) > 1e-13){
         L[usedL] = i;
         usedL++;
      }
   }
   return;
}

void findJtilde(iReg *Jtilde, iReg &JtildeSize, iReg *L, iReg sizeL, iExt *iatk, iReg *jak, iReg *J, iReg sizeJ){

   // Sanity check for sizeL == 0
   if (sizeL == 0){
      printf("sizeL == 0, check\n");
      return;
   }

   // Initialize counter to 0
   JtildeSize = 0;

   // Flag for fast exit in case of repeated index
   bool skip = false;
   iReg maxJsize;
   for (iReg i = 0; i < sizeL; ++i){
      for (iReg j = iatk[L[i]]; j < iatk[L[i] + 1]; ++j){
         // Get the max size to search for duplicates
         maxJsize = std::max(sizeJ,JtildeSize);
         // std::cout << "maxJsize = " << maxJsize << std::endl;
         // std::cout << "ja["<< j << "] = " << jak[j] << std::endl;
         // Check for duplicate value
         for (iReg q = 0; q < maxJsize; ++q){
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
void fullAJtilde(iExt nn_A, iExt *iatk, iReg *jak, double *coefk, iReg *Jtilde, iReg JtildeSize, double *AJtilde){
   
   std::fill_n(AJtilde, nn_A * JtildeSize, 0.0);

   // Find both the minimum and maximum column indices in Jtilde
   auto minmax = std::minmax_element(Jtilde, Jtilde + JtildeSize);
   iReg min_col = *minmax.first;
   iReg max_col = *minmax.second;
   iReg range = max_col - min_col + 1;

   // Fill the lookup array based purely on the range span
   std::vector<iReg> lookup(range, -1);
   for (iReg c = 0; c < JtildeSize; ++c) {
       lookup[Jtilde[c] - min_col] = c; // Apply the negative offset
   }

   // Iterate through the CSR matrix
   for (iReg r = 0; r < nn_A; ++r) {
       iReg row_end = iatk[r + 1];
       for (iExt k = iatk[r]; k < row_end; ++k) {
           iReg col = jak[k];
           
           // Quick bounds check using the min/max cluster boundaries
           if (col >= min_col && col <= max_col) {
               iReg c = lookup[col - min_col]; // Apply same offset to query
               
               if (c != -1) {
                   // Compute column-major index: row + (col_index * total_rows)
                   AJtilde[r + c * nn_A] = coefk[k];
               }
           }
       }
   }
   // // Loop over possible columns for Jtilde
   // for (iReg j = 0; j < JtildeSize; ++j){
   //    // Loop over the rows
   //    for (iReg row = 0; row < nn_A; ++row){
   //       // Loop over the single row entries
   //       for (iReg i = iatk[row]; i < iatk[row+1]; ++i){
         
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
