#include "cpt_resRho.h"
#include "precision.h"
#include <vector>
#include <cstring>
#include <iostream>
#include <algorithm> // For std::max and std::min_element
#include <iterator>  // For std::distance
#include <cmath>
#include "lapack.h"


void cptRhoJ2(const iReg JtildeSize,
              double * RESTRICT normColJ,
              const iExt* RESTRICT jatAJtilde,
              const iReg * const * RESTRICT iaAJtilde,
              const double * const * RESTRICT coefAJtilde,
              const double * RESTRICT res,
              const double normRes, const double * RESTRICT colANorm,
              const iReg * RESTRICT Jtilde) {

    const double normResSq = normRes * normRes;

    for (iReg i = 0; i < JtildeSize; ++i) {
        const iReg* RESTRICT ptr_ia = iaAJtilde[i];
        const double* RESTRICT ptr_coef = coefAJtilde[i];

        const iExt num_elements = jatAJtilde[i + 1] - jatAJtilde[i];

        // 8 independent accumulators to break latency chains
        double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
        double sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;

        // Process blocks of 8 using pointer offsets
        const iReg* end_block = ptr_ia + (num_elements & ~7);

        while (ptr_ia < end_block) {
            // Fetch 8 indices
            iReg r0 = ptr_ia[0];  iReg r1 = ptr_ia[1];
            iReg r2 = ptr_ia[2];  iReg r3 = ptr_ia[3];
            iReg r4 = ptr_ia[4];  iReg r5 = ptr_ia[5];
            iReg r6 = ptr_ia[6];  iReg r7 = ptr_ia[7];

            // Issue 8 parallel loads & multiplications
            sum0 += ptr_coef[0] * res[r0];
            sum1 += ptr_coef[1] * res[r1];
            sum2 += ptr_coef[2] * res[r2];
            sum3 += ptr_coef[3] * res[r3];
            sum4 += ptr_coef[4] * res[r4];
            sum5 += ptr_coef[5] * res[r5];
            sum6 += ptr_coef[6] * res[r6];
            sum7 += ptr_coef[7] * res[r7];

            // Single pointer increment per 8 elements
            ptr_ia += 8;
            ptr_coef += 8;
        }

        // Tree reduction for minimum floating-point addition latency
        double sum = ((sum0 + sum1) + (sum2 + sum3)) +
                     ((sum4 + sum5) + (sum6 + sum7));

        // Tail loop for remaining elements (0 to 7)
        const iReg* end_ia = iaAJtilde[i] + num_elements;
        while (ptr_ia < end_ia) {
            sum += ptr_coef[0] * res[ptr_ia[0]];
            ptr_ia++;
            ptr_coef++;
        }

        normColJ[i] = normResSq - (sum * sum) / colANorm[Jtilde[i]];
    }
}



// Find the index of the minimum 
iReg minIdx(double *rhoJ2, iReg JtildeSize){
    if (JtildeSize <= 0) return 0;
    
    double* min_element_ptr = std::min_element(rhoJ2, rhoJ2 + JtildeSize);
    return static_cast<iReg>(std::distance(rhoJ2, min_element_ptr));
}

void cptRes(iReg sizeJ, const iExt * RESTRICT jatAJ, 
            const iReg * const * RESTRICT iaAJ, 
            const double * const * RESTRICT coefAJ, 
            const double * RESTRICT mHat,
            const iReg * RESTRICT A0k_idx, const double * RESTRICT A0k, iReg A0k_nnz,
            double * RESTRICT res, iReg * RESTRICT L, iReg &usedL,
            double &resRelNorm, double &resNorm,
            int* RESTRICT ws_idx, double* RESTRICT ws_val) {

    // Clear the elements of 'res' that were modified 
    // in the previous call, using the old L array.
    for (iReg n = 0; n < usedL; ++n) {
        res[L[n]] = 0.0;
        L[n] = 0;
    }

    int ws_count = 0;

    // Compute AJ * mHat and save into the workspace
    for (iReg j = 0; j < sizeJ; ++j) {
        const double m_j = mHat[j];
        if (m_j == 0.0) continue;

        // Get the direct pointers to this column's row indices and coefficients
        const iReg* RESTRICT col_ia = iaAJ[j];
        const double* RESTRICT col_coef = coefAJ[j];

        // Determine how many non-zero elements are in this column
        iExt start = jatAJ[j];
        iExt end   = jatAJ[j + 1];
        iExt num_elements = end - start;

        // Loop over the elements locally (0 to num_elements)
        for (iExt k = 0; k < num_elements; ++k) {
            iReg row = col_ia[k];
            if (ws_val[row] == 0.0) {
                ws_idx[ws_count++] = row;
            }
            ws_val[row] += col_coef[k] * m_j;
        }
    }

    // Compute the norm looking only at the full entries
    double normAjMh_sq = 0.0;
    for (int n = 0; n < ws_count; ++n) {
        double spmv_val = ws_val[ws_idx[n]];
        normAjMh_sq += spmv_val * spmv_val;
    }

    // Subtract A0k from the workspace
    for (iReg n = 0; n < A0k_nnz; ++n) {
        iReg row = A0k_idx[n];
        if (ws_val[row] == 0.0) {
            ws_idx[ws_count++] = row;
        }
        ws_val[row] -= A0k[row]; 
    }

    // Populate res, generate L, and compute total sq_sum
    double sq_sum = 0.0;
    iReg local_usedL = 0;

    for (int n = 0; n < ws_count; ++n) {
        iReg row = ws_idx[n];
        double final_res = ws_val[row];

        // If the values didn't perfectly cancel out to 0, it belongs in L
        if (final_res != 0.0) {
            res[row] = final_res;
            L[local_usedL++] = row;
            sq_sum += final_res * final_res;
        }

        // Reset workspace to 0.0 for the next call
        ws_val[row] = 0.0;
    }

    // Update the reference with the new count of non-zeros
    usedL = local_usedL;

    // Final Global Math
    double normAjMh = std::sqrt(normAjMh_sq);
    resNorm = std::sqrt(sq_sum);
    std::sort(L, L + usedL);

    resRelNorm = 2.0 * resNorm / (normAjMh + resRelNorm);
}


// Sifts down the root of a max-heap of size k
inline void replace_max_heap_top(MinCandidate* heap, iReg k, double new_val, iReg new_j) {
    iReg i = 0;
    while (true) {
        iReg left = 2 * i + 1;
        iReg right = 2 * i + 2;
        iReg largest = i;

        if (left < k && heap[largest].val < heap[left].val) {
            largest = left;
        }
        if (right < k && heap[largest].val < heap[right].val) {
            largest = right;
        }
        if (largest == i) break;

        heap[i] = heap[largest];
        i = largest;
    }
    heap[i] = { new_val, new_j };
}

void multiMinIdx(iReg step_size, iReg &JtildeSize, const iReg *Jtilde, const double *rhoJ2, iReg *Jstart) {
   step_size = std::min(step_size, JtildeSize);
   if (step_size <= 0) return;

   // Fixed stack buffer to prevent dynamic memory allocation
   constexpr size_t STACK_CAPACITY = 64;
   MinCandidate local_buffer[STACK_CAPACITY];

   // Safety fallback: if step_size ever exceeds stack capacity
   MinCandidate* heap = local_buffer;
   std::vector<MinCandidate> fallback_vec;
   if (static_cast<size_t>(step_size) > STACK_CAPACITY) {
      fallback_vec.resize(step_size);
      heap = fallback_vec.data();
   }

   // Initialize max-heap with the first step_size elements
   for (iReg i = 0; i < step_size; ++i) {
      heap[i] = { rhoJ2[i], Jtilde[i] };
   }
   std::make_heap(heap, heap + step_size);

   // Single pass over the remaining elements (N - k iterations)
   for (iReg i = step_size; i < JtildeSize; ++i) {
      double current_val = rhoJ2[i];
      // Fast rejection filter: only update if current_val is smaller than maximum in heap
      if (current_val < heap[0].val) {
         replace_max_heap_top(heap, step_size, current_val, Jtilde[i]);
      }
   }
 
   // Sort the selected candidates in ascending order by column index
   std::sort(heap, heap + step_size, [](const MinCandidate& a, const MinCandidate& b) {   
      return a.j_val < b.j_val; 
   });

   // Output values into Jstart
   for (iReg i = 0; i < step_size; ++i) {
      Jstart[i] = heap[i].j_val;
   }
   JtildeSize = step_size;
}
