#include "find_stuff.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <iostream>
#include <unordered_set>

void fullA0k(const iExt nn_A, const iExt * RESTRICT iat0,
             const iReg * RESTRICT ja0, const double * RESTRICT coef0,
             const iExt k, double * RESTRICT A0k,
             iReg * RESTRICT A0k_idx, iReg &A0k_nnz) {

   // Clear the elements modified in the previous call
   for (iReg n = 0; n < A0k_nnz; ++n) {
      A0k[A0k_idx[n]] = 0.0;
   }

   // Direct lookup of the start and end bounds for column k
   iExt col_start = iat0[k];
   iExt col_end   = iat0[k+1];

   iReg local_nnz = 0;

   // Populate the active rows and log their indices
   for (iExt idx = col_start; idx < col_end; ++idx) {
      iReg row = ja0[idx];
      A0k[row] = coef0[idx];
      A0k_idx[local_nnz++] = row;
   }

   // Record how many non-zero elements we found in this column
   A0k_nnz = local_nnz;
}

void findNonZeroInColJ(const iReg *RESTRICT J, const iExt *RESTRICT iatk,
                       const iReg *RESTRICT jak, const iReg n2, 
                       iReg *RESTRICT I, iReg &sizeI,
                       iReg *RESTRICT visited, const int t) {
    
    // Mark pre-existing elements of I as visited in the current t
    for (iReg k = 0; k < sizeI; ++k) {
        visited[I[k]] = t;
    }

    for (iReg i = 0; i < n2; ++i) {
        const iReg start = iatk[J[i]];
        const iReg end = iatk[J[i] + 1];
        
        for (iReg j = start; j < end; ++j) {
            const iReg val = jak[j];
            if (visited[val] != t) {
                visited[val] = t;
                I[sizeI++] = val;
            }
        }
    }
}

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, iExt *iat0, iReg *ja0, double *coef0, iExt k){
   // Look up column k boundaries
   iExt col_start = iat0[k];
   iExt col_end   = iat0[k+1];
   iReg col_len   = col_end - col_start;
   iReg *col_rows = &ja0[col_start];

   // Loop over the new rows
   for (iReg i = oldSizeI; i < sizeI; ++i) {
      iReg row = I[i];
      a0k[i] = 0.0; // Default value

      // Binary search for the 'row' within the contiguous rows of column k
      if (col_len > 0) {
         auto it = std::lower_bound(col_rows, col_rows + col_len, row);
         if (it != col_rows + col_len && *it == row) {
            a0k[i] = coef0[col_start + (it - col_rows)];
         }
      }
   }
}

void getAhat(iReg * RESTRICT I, iReg sizeI, iReg * RESTRICT J, iReg Jstart, 
             iReg Jend, iExt * RESTRICT iatk, iReg * RESTRICT jak, 
             double * RESTRICT coefk, double * RESTRICT Ahat, iReg &Astart) {

   // Cycle over the columns in J that have been added
   for (iReg j = Jstart; j < Jend; ++j){
      iReg colJ = J[j];
      // Cycle over the rows in I
      for (iReg i = 0; i < sizeI; ++i){
         iExt k, row = I[i];
         const iExt row_start = iatk[row];
         const iExt row_end = iatk[row+1];

         // Cycle over the chosen row
         for(k = row_start; k < row_end; ++k){
            // If the column in the row coincides with the column added then get the nonzero value
            if(jak[k] == colJ){
               Ahat[Astart] = coefk[k];
               Astart++;
               break;
            }
         }

         if(k == row_end){
            // If entered here there is no column entry equal to k
            // set to zero
            Ahat[Astart] = 0;
            Astart++;
         }
      }
   }
   return;
}

// Get the new columns by storing pointers to the values inside matrix Ak
void getAJ(iReg *J, iReg Jstart, iReg Jend, iExt *jatk, iReg *iak, double *coefk,
           const iReg **iaAJ, const double **coefAJ, iExt *jatAJ) {

   // Find starting point (virtual NNZ tracking)
   iExt current_nnz = jatAJ[Jstart];

   // Loop through the requested column indices in J
   for (iReg i = Jstart; i < Jend; ++i) {
      // Get the source column index from matrix Ak
      iReg col_A = J[i];
      
      // Find the start and end of this column inside Ak
      iExt start_A = jatk[col_A];
      iExt end_A   = jatk[col_A + 1];
      iExt num_elements = end_A - start_A;

      // ZERO COPY: Store the pointer directly to the slice inside Ak
      // Note: We index these by the column offset 'i' instead of 'current_nnz'
      iaAJ[i]   = &iak[start_A];
      coefAJ[i] = &coefk[start_A];
    
      // Keep track of the "virtual" cumulative NNZ
      current_nnz += num_elements;

      // Update the column pointer for the next column of AJ
      jatAJ[i + 1] = current_nnz;
   }
}


void findJtilde(iReg *Jtilde, iReg &JtildeSize,
                const iReg* RESTRICT L, const iReg sizeL,
                const iExt* RESTRICT iatk, const iReg* RESTRICT jak,
                const iReg* RESTRICT J, const iReg sizeJ,
                uint8_t* RESTRICT seen) {
    if (sizeL == 0) {
        printf("sizeL == 0, check\n");
        return;
    }

    JtildeSize = 0;

    // Fill the array with existing J elements
    for (iReg i = 0; i < sizeJ; ++i) {
        seen[J[i]] = 1;
    }

    // Main processing loop
    for (iReg i = 0; i < sizeL; ++i) {
        const iReg l_idx = L[i];
        const iExt start = iatk[l_idx];
        const iExt end   = iatk[l_idx + 1];

        for (iExt j = start; j < end; ++j) {
            const iReg val = jak[j];
            if (!seen[val]) {
                seen[val] = 1;
                Jtilde[JtildeSize++] = val;
            }
        }
    }

    // Selective zeroing
    for (iReg i = 0; i < sizeJ; ++i) {
        seen[J[i]] = 0;
    }
    for (iReg i = 0; i < JtildeSize; ++i) {
        seen[Jtilde[i]] = 0;
    }
}
