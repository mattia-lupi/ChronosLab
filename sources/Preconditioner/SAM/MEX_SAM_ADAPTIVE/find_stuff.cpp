#include "find_stuff.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <iostream>
#include <unordered_set>

void fullA0k(const iExt nn_A, const iExt* RESTRICT iat0,
             const iReg* RESTRICT ja0, const double* RESTRICT coef0,
             const iExt k, double *A0k){
   // Initialize the dense output vector with zeros
   std::fill_n(A0k, nn_A, 0.0);

   // Direct lookup of the start and end bounds for column k
   iExt col_start = iat0[k];
   iExt col_end   = iat0[k+1];

   // Populate only the rows that have non-zero entries in this column
   for (iExt idx = col_start; idx < col_end; ++idx) {
      iReg row = ja0[idx];
      A0k[row] = coef0[idx];
   }
}

void findNonZeroInColJ(const iReg* RESTRICT J, const iExt* RESTRICT iatk, 
                       const iReg* RESTRICT jak, const iReg n2, iReg *I, iReg &sizeI){
   std::unordered_set<iReg> seen(I, I + sizeI);

   for (iReg i = 0; i < n2; ++i) {
      iReg start = iatk[J[i]];
      iReg end = iatk[J[i] + 1];
      for (iReg j = start; j < end; ++j) {
         iReg val = jak[j];
         if (seen.insert(val).second) {
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

void getAhat(iReg *I, iReg sizeI, iReg *J, iReg Jstart, iReg Jend,
             iExt *iatk, iReg *jak, double *coefk, double *Ahat, iReg &Astart) {

   // Cycle over the columns in J that have been added
   for (iReg j = Jstart; j < Jend; ++j){
      // Cycle over the rows in I
      for (iReg i = 0; i < sizeI; ++i){
         // Cycle over the chosen row
         for(iExt k = iatk[I[i]]; k < iatk[I[i]+1]; ++k){
            // If the column in the row coincides with the column added then get the nonzero value
            if(jak[k] == J[j]){
               Ahat[Astart] = coefk[k];
               Astart++;
               break;
            }

            if(k == iatk[I[i] + 1] - 1){
               // If entered here there is no column entry equal to k
               // set to zero
               Ahat[Astart] = 0;
               Astart++;
            }
         }
      }
   }
   return;
}

// Get the new columns and add them to AJ
void getAJ(iReg *J, iReg Jstart, iReg Jend, iExt *jatk,iReg *iak,double *coefk,
           iExt *jatAJ, iReg *iaAJ, double *coefAJ) {

   // Find starting point
   iExt current_nnz = jatAJ[Jstart];

   // Loop through the requested column indices in J
   for (iReg i = Jstart; i < Jend; ++i) {
      // Get the source column index from matrix Ak
      iReg col_A = J[i];
      
      // Find the start and end of this column inside Ak
      iExt start_A = jatk[col_A];
      iExt end_A   = jatk[col_A + 1];
      iExt num_elements = end_A - start_A;

      // Copy row indices
      std::memcpy(&iaAJ[current_nnz], &iak[start_A], num_elements * sizeof(iReg));
    
      // Copy coefficients
      std::memcpy(&coefAJ[current_nnz], &coefk[start_A], num_elements * sizeof(double));
    
      // Get new nnz
      current_nnz += num_elements;

      // Update the column pointer for the next column of AJ
      jatAJ[i + 1] = current_nnz;
   }
}


void fillL(iReg * RESTRICT L, const double * RESTRICT res, iExt nn_A, iReg &usedL) {
   iReg local_usedL = 0;

   // Loop over all the rows
   for (iExt i = 0; i < nn_A; ++i) {
      if (res[i] != 0.0) { 
         L[local_usedL] = static_cast<iReg>(i);
         local_usedL++;
      }
   }

   // Write back to the reference exactly once
   usedL = local_usedL; 
}

void findJtilde(iReg *Jtilde, iReg &JtildeSize, const iReg* RESTRICT L,
                const iReg sizeL, const iExt* RESTRICT iatk, const iReg* RESTRICT jak, 
                iReg *J, iReg sizeJ) {
   if (sizeL == 0) {
      printf("sizeL == 0, check\n");
      return;
   }

   JtildeSize = 0;
   iReg l_idx, start, end, val;

   // Find the maximum ID to size our tracking array
   iReg max_val = 0;
   iReg jj;
   for (iReg i = 0; i < sizeJ; ++i) {
      jj = J[i];
      if (jj > max_val) max_val = jj;
   }
   for (iReg i = 0; i < sizeL; ++i) {
      l_idx = L[i];
      start = iatk[l_idx];
      end = iatk[l_idx + 1];
      for (iReg j = start; j < end; ++j) {
         if (jak[j] > max_val) max_val = jak[j];
      }
   }

   // Use a flat vector lookups
   std::vector<bool> seen(max_val + 1, false);

   // Seed the array with existing J elements
   for (iReg i = 0; i < sizeJ; ++i) {
      seen[J[i]] = true;
   }

   // Loop over the L
   for (iReg i = 0; i < sizeL; ++i) {
      l_idx = L[i]; 
      start = iatk[l_idx];
      end = iatk[l_idx + 1];

      // Loop over the row L[i]
      for (iReg j = start; j < end; ++j) {
         val = jak[j];

         if (!seen[val]) {
            seen[val] = true;
            Jtilde[JtildeSize++] = val;
         }
      }
   }
}

