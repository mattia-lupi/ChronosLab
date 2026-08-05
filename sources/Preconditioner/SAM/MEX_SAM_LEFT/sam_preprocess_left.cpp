#include "sam_preprocess_left.h"
#include <algorithm>
#include <cstdint>
#include <omp.h>

void sam_preprocess_left(int n,
                         const int* a_ptr, const int* a_col,
                         const int* s_ptr, const int* s_data,
                         std::vector<int>& r_ptr,
                         std::vector<int>& r_data,
                         int num_threads)
{
    if (n <= 0) return;

    int total_a_nnz = a_ptr[n];
    int max_col = n - 1;
    for (int p = 0; p < total_a_nnz; ++p) {
        if (a_col[p] > max_col) max_col = a_col[p];
    }
    int marker_size = max_col + 1;

    r_ptr.resize(n + 1, 0);
    std::vector<int> r_counts(n, 0);

    int active_threads = (num_threads > 0) ? num_threads : omp_get_max_threads();

    // Pass 1: Count non-zeros per row in R
    #pragma omp parallel num_threads(active_threads)
    {
        std::vector<uint8_t> marker(marker_size, 0);
        std::vector<int> local_cols;
        local_cols.reserve(256);

        #pragma omp for schedule(guided)
        for (int i = 0; i < n; ++i) {
            local_cols.clear();
            int s0 = s_ptr[i], s1 = s_ptr[i + 1];

            for (int ji = s0; ji < s1; ++ji) {
                int j = s_data[ji];
                if (j < 0 || j >= n) continue;

                int a0 = a_ptr[j], a1 = a_ptr[j + 1];
                for (int p = a0; p < a1; ++p) {
                    int c = a_col[p];
                    // Bug Fix 1: Guard against out-of-bounds column indices
                    if (c < 0 || c >= marker_size) continue;

                    if (!marker[c]) {
                        marker[c] = 1;
                        local_cols.push_back(c);
                    }
                }
            }

            for (int c : local_cols) marker[c] = 0;
            r_counts[i] = static_cast<int>(local_cols.size());
        }
    }

    r_ptr[0] = 0;
    for (int i = 0; i < n; ++i) {
        r_ptr[i + 1] = r_ptr[i] + r_counts[i];
    }

    int total_r_nnz = r_ptr[n];
    r_data.resize(total_r_nnz);

    // Pass 2: Populate column indices
    #pragma omp parallel num_threads(active_threads)
    {
        std::vector<uint8_t> marker(marker_size, 0);
        std::vector<int> local_cols;
        local_cols.reserve(256);

        #pragma omp for schedule(guided)
        for (int i = 0; i < n; ++i) {
            local_cols.clear();
            int s0 = s_ptr[i], s1 = s_ptr[i + 1];

            for (int ji = s0; ji < s1; ++ji) {
                int j = s_data[ji];
                if (j < 0 || j >= n) continue;

                int a0 = a_ptr[j], a1 = a_ptr[j + 1];
                for (int p = a0; p < a1; ++p) {
                    int c = a_col[p];
                    if (c < 0 || c >= marker_size) continue;

                    if (!marker[c]) {
                        marker[c] = 1;
                        local_cols.push_back(c);
                    }
                }
            }

            for (int c : local_cols) marker[c] = 0;

            std::sort(local_cols.begin(), local_cols.end());
            int offset = r_ptr[i];
            for (size_t k = 0; k < local_cols.size(); ++k) {
                r_data[offset + k] = local_cols[k];
            }
        }
    }
}
