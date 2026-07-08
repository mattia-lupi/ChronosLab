#include "find_stuff.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <unordered_set>

bool mattia = true;
void fullA0k(iExt nn_A, iExt *iat0, iReg *ja0, double *coef0, iExt k, double *A0k){
   std::fill_n(A0k, nn_A, 0.0);

   for (iReg row = 0; row < nn_A; ++row) {
      iReg row_start = iat0[row];
      iReg row_end = iat0[row+1];
      iReg row_len = row_end - row_start;

      if (row_len > 0) {
         iReg *row_cols = &ja0[row_start];
         auto it = std::lower_bound(row_cols, row_cols + row_len, k);
         if (it != row_cols + row_len && *it == k) {
            A0k[row] = coef0[row_start + (it - row_cols)];
         }
      }
   }
}

void findNonZeroInColJ(iReg *J, iExt *iatk, iReg *jak, iReg n2, iReg *I, iReg &sizeI){
   // Hash set to track unique items in O(1) time
   std::unordered_set<iReg> seen(I, I + sizeI);

   for (iReg i = 0; i < n2; ++i) {
      iReg start = iatk[J[i]];
      iReg end = iatk[J[i] + 1];
      for (iReg j = start; j < end; ++j) {
         iReg val = jak[j];
         // insert().second is true only if the item didn't already exist in the set
         if (seen.insert(val).second) {
            I[sizeI++] = val;
         }
      }
   }
}

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, iExt *iat0, iReg *ja0, double *coef0, iExt k){
   for (iReg i = oldSizeI; i < sizeI; ++i) {
      iReg row = I[i];
      iReg row_start = iat0[row];
      iReg row_end = iat0[row + 1];
      iReg row_len = row_end - row_start;

      a0k[i] = 0.0; // Default value

      if (row_len > 0) {
         iReg *row_cols = &ja0[row_start];
         auto it = std::lower_bound(row_cols, row_cols + row_len, k);
         if (it != row_cols + row_len && *it == k) {
            a0k[i] = coef0[row_start + (it - row_cols)];
         }
      }
   }
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

// Get the new columns and add them to AJ
void getAJ(iReg *J, iReg Jstart, iReg Jend, iExt nn_A, iExt *iatk, iReg *jak, double *coefk, double *AJ) {
    // Initialize the newly added columns in AJ to zero
    for (iReg j = Jstart; j < Jend; ++j) {
        for (iExt i = 0; i < nn_A; ++i) {
            AJ[j * nn_A + i] = 0.0;
        }
    }

    // Populate AJ with values from the CSR matrix A
    for (iExt i = 0; i < nn_A; ++i) {
        iExt row_start = iatk[i];
        iExt row_end = iatk[i + 1];
        
        for (iExt k = row_start; k < row_end; ++k) {
            iReg col_A = jak[k];
            
            // Search if the current column index exists in the active J window
            for (iReg j = Jstart; j < Jend; ++j) {
                if (J[j] == col_A) {
                    // Column-major indexing: column * structural_rows + row
                    AJ[j * nn_A + i] = coefk[k];
                }
            }
        }
    }
}


void fillL(iReg *L, double *res, iExt nn_A, iReg &usedL){
   usedL = 0;
   // Loop over all residual entries
   for (iReg i = 0; i < nn_A; ++i){
      // The residual is not numerically zero
      // Add it as a possible column to be computed
      if (std::abs(res[i]) > 0){
         L[usedL] = i;
         usedL++;
      }
   }
   return;
}

void findJtilde(iReg *Jtilde, iReg &JtildeSize, iReg *L, iReg sizeL, iExt *iatk, iReg *jak, iReg *J, iReg sizeJ){
   if (sizeL == 0){
      printf("sizeL == 0, check\n");
      return;
   }

   JtildeSize = 0;
   // Seed the set with existing J elements so we don't include them in Jtilde
   std::unordered_set<iReg> seen(J, J + sizeJ);

   for (iReg i = 0; i < sizeL; ++i){
      iReg start = iatk[L[i]];
      iReg end = iatk[L[i] + 1];
      for (iReg j = start; j < end; ++j){
         iReg val = jak[j];
         if (seen.insert(val).second) {
            Jtilde[JtildeSize++] = val;
         }
      }
   }
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
   return;
}
