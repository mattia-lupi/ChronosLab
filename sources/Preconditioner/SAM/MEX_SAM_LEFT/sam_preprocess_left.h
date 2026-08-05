#ifndef SAM_PREPROCESS_LEFT_H
#define SAM_PREPROCESS_LEFT_H

#include <vector>
#include <cstdint>

// Performs preprocessing to assemble the sparsity pattern R for SAM.
// Inputs are 0-based CSR matrices A and S.
void sam_preprocess_left(int n,
                         const int* a_ptr, const int* a_col,
                         const int* s_ptr, const int* s_data,
                         std::vector<int>& r_ptr,
                         std::vector<int>& r_data,
                         int num_threads = 0);

#endif // SAM_PREPROCESS_LEFT_H
