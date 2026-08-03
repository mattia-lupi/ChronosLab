#include "find_stuff.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <cmath>

double fullA0k(const iExt * RESTRICT iat0,
             const iReg * RESTRICT ja0, const double * RESTRICT coef0,
             const iExt k, double * RESTRICT A0k,
             iReg * RESTRICT A0k_idx, iReg &A0k_nnz) {

   // Clear the elements modified in the previous call
   for (iReg n = 0; n < A0k_nnz; ++n) {
      A0k[A0k_idx[n]] = 0.0;
      A0k_idx[n] = 0;
   }

   // Direct lookup of the start and end bounds for column k
   iExt col_start = iat0[k];
   iExt col_end   = iat0[k+1];

   iReg local_nnz = 0;
   double sq_sum  = 0.0;

   // Populate the active rows and log their indices
   for (iExt idx = col_start; idx < col_end; ++idx) {
      // Get "handles"
      iReg row = ja0[idx];
      const double val = coef0[idx];

      // Fill dense vector
      A0k[row] = val;
      A0k_idx[local_nnz++] = row;

      // Compute the norm
      sq_sum += val * val;
   }

   // Record how many non-zero elements we found in this column
   A0k_nnz = local_nnz;
   return std::sqrt(sq_sum);
}

void findNonZeroInColJ(const iReg *RESTRICT J, const iExt *RESTRICT iatk,
                       const iReg *RESTRICT jak, const iReg n2, 
                       iReg *RESTRICT I, iReg &sizeI,
                       iReg *RESTRICT visited, const int t) {
    
   // Mark pre-existing elements of I as visited in the current t
   for (iReg k = 0; k < sizeI; ++k) {
      visited[I[k]] = t;
   }

   // Loop over all J
   for (iReg i = 0; i < n2; ++i) {
      // Get handles
      iReg row = J[i];
      const iReg start = iatk[row];
      const iReg end = iatk[row + 1];
      
      // Loop over row J[i]
      for (iReg j = start; j < end; ++j) {
         const iReg val = jak[j];
         // If not seen before, add it
         if (visited[val] != t) {
            visited[val] = t;
            I[sizeI++] = val;
         }
      }
   }
}

void getA0k(double *a0k, iReg *I, iReg sizeI, iReg oldSizeI, double* fullA0k) {

    // Loop over the new additions
    for (iReg i = oldSizeI; i < sizeI; ++i) {
        // Get the value from the full vector
        a0k[i] = fullA0k[I[i]];
    }
}

void getAhat(iReg * RESTRICT I, iReg sizeI, 
             iReg * RESTRICT J, iReg Jstart, iReg Jend, 
             iExt * RESTRICT iatk, iReg * RESTRICT jak,
             double * RESTRICT coefk, double * RESTRICT Ahat, 
             iReg &Astart) 
{
    const iReg nJ = Jend - Jstart;
    const iReg base_Astart = Astart;

    // Cycle over rows I on the outer loop
    for (iReg i = 0; i < sizeI; ++i) {
        const iExt row = I[i];
        const iExt row_start = iatk[row];
        const iExt row_end   = iatk[row + 1];

        // Cycle over columns J on the inner loop
        for (iReg j = Jstart; j < Jend; ++j) {
            const iReg colJ = J[j];

            // Compute the target index in Ahat to preserve colmajor ordering
            const iReg dest_idx = base_Astart + (j - Jstart) * sizeI + i;

            iExt k;
            // Check if the current row has entry in the current column
            for (k = row_start; k < row_end; ++k) {
                if (jak[k] == colJ) {
                    Ahat[dest_idx] = coefk[k];
                    break;
                }
            }

            if (k == row_end) {
                Ahat[dest_idx] = 0.0;
            }
        }
    }

    // Advance Astart by total entries written
    Astart += sizeI * nJ;
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

      // Store the pointer directly to the slice inside Ak
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
