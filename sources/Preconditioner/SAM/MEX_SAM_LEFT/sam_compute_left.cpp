#include "sam_compute_left.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>

// Optimized Cholesky solver with Tikhonov Regularization
inline bool fast_cholesky_solve_reg(int n, double* B, const double* g, double* out_z,
                                    double* L, double* y)
{
    // Inline factorize step eliminates code duplication between standard and regularized passes
    auto factorize = [&](double min_pivot) -> bool {
        for (int i = 0; i < n; ++i) {
            const int i_n = i * n;
            for (int j = 0; j <= i; ++j) {
                const int j_n = j * n;
                double dot = 0.0;
                #pragma omp simd reduction(+:dot)
                for (int k = 0; k < j; ++k) {
                    dot += L[i_n + k] * L[j_n + k];
                }
                const double sum = B[i_n + j] - dot;

                if (i == j) {
                    if (sum <= min_pivot || !std::isfinite(sum)) return false;
                    L[i_n + j] = std::sqrt(sum);
                } else {
                    L[i_n + j] = sum / L[j_n + j];
                }
            }
        }
        return true;
    };

    // Pass 1: Standard Cholesky decomposition (B = L * L^T)
    if (!factorize(1e-15)) {
        // Pass 2: Tikhonov Regularization fallback if matrix is ill-conditioned
        double trace = 0.0;
        for (int i = 0; i < n; ++i) trace += std::abs(B[i * n + i]);
        const double eps = std::max(1e-12, 1e-8 * (trace / std::max(1, n)));

        for (int i = 0; i < n; ++i) B[i * n + i] += eps;

        if (!factorize(0.0)) return false;
    }

    // Forward substitution: L * y = g
    for (int i = 0; i < n; ++i) {
        const int i_n = i * n;
        double dot = 0.0;
        #pragma omp simd reduction(+:dot)
        for (int k = 0; k < i; ++k) {
            dot += L[i_n + k] * y[k];
        }
        y[i] = (g[i] - dot) / L[i_n + i];
    }

    // Backward substitution directly into out_z: L^T * out_z = y (eliminates z workspace array copy)
    for (int i = n - 1; i >= 0; --i) {
        const int i_n = i * n;
        double sum = y[i];
        for (int k = i + 1; k < n; ++k) {
            sum -= L[k * n + i] * out_z[k];
        }
        const double val = sum / L[i_n + i];
        if (!std::isfinite(val)) return false;
        out_z[i] = val;
    }

    return true;
}

void sam_compute_left(int n,
                      const int* ak_ptr, const int* ak_col, const double* ak_val,
                      const int* a0_ptr, const int* a0_col, const double* a0_val,
                      const int* s_ptr,  const int* s_data,
                      const int* r_ptr,  const int* r_data,
                      double* out_val,
                      int num_threads)
{
    if (n <= 0) return;
    const int active_threads = (num_threads > 0) ? num_threads : omp_get_max_threads();

    #pragma omp parallel num_threads(active_threads)
    {
        // Thread-local persistent workspace buffers (allocated ONCE per thread outside OpenMP loop)
        std::vector<int> col_ptr;
        std::vector<int> nz_row;
        std::vector<double> nz_val;
        std::vector<int> f_nz_idx;
        std::vector<double> f_dense;
        std::vector<double> w_dense;

        std::vector<double> B_vec;
        std::vector<double> g_vec;
        std::vector<double> L_vec;
        std::vector<double> y_vec;

        // Pre-reserve capacity to eliminate dynamic heap reallocations during push_back
        nz_row.reserve(512);
        nz_val.reserve(512);
        f_nz_idx.reserve(128);

        #pragma omp for schedule(guided)
        for (int i = 0; i < n; ++i) {
            const int s0 = s_ptr[i], s1 = s_ptr[i + 1], nnzi = s1 - s0;
            if (nnzi == 0) continue;

            const int r0 = r_ptr[i], r1 = r_ptr[i + 1], nri = r1 - r0;
            if (nri == 0) {
                std::fill_n(&out_val[s0], nnzi, 0.0);
                continue;
            }

            const int m_ls = nri, n_ls = nnzi;

            // Expand buffers if current subproblem exceeds thread workspace size
            if (static_cast<int>(w_dense.size()) < m_ls) {
                w_dense.resize(m_ls, 0.0);
                f_dense.resize(m_ls, 0.0);
            }
            if (static_cast<int>(col_ptr.size()) < n_ls + 1) {
                col_ptr.resize(n_ls + 1);
            }
            if (static_cast<int>(g_vec.size()) < n_ls) {
                B_vec.resize(n_ls * n_ls);
                L_vec.resize(n_ls * n_ls);
                g_vec.resize(n_ls);
                y_vec.resize(n_ls);
            }

            const int* si = s_data + s0;
            const int* ri = r_data + r0;

            // 1. Sparse Extraction of Vector f
            f_nz_idx.clear();
            const int as = a0_ptr[i], ae = a0_ptr[i + 1];
            int p_a0 = as, p_ri = 0;
            while (p_a0 < ae && p_ri < m_ls) {
                const int col_a0 = a0_col[p_a0];
                const int col_r  = ri[p_ri];
                if (col_a0 == col_r) {
                    f_dense[p_ri] = a0_val[p_a0];
                    f_nz_idx.push_back(p_ri);
                    p_a0++; p_ri++;
                } else if (col_a0 < col_r) {
                    p_a0++;
                } else {
                    p_ri++;
                }
            }

            // 2. Sparse Extraction of Matrix Atmp (CSC Format)
            nz_row.clear();
            nz_val.clear();
            col_ptr[0] = 0;

            for (int ji = 0; ji < n_ls; ++ji) {
                const int row_k = si[ji];
                if (row_k >= 0 && row_k < n) {
                    const int ks = ak_ptr[row_k], ke = ak_ptr[row_k + 1];
                    int p_ak = ks, p_r = 0;
                    while (p_ak < ke && p_r < m_ls) {
                        const int col_ak = ak_col[p_ak];
                        const int col_r  = ri[p_r];
                        if (col_ak == col_r) {
                            nz_row.push_back(p_r);
                            nz_val.push_back(ak_val[p_ak]);
                            p_ak++; p_r++;
                        } else if (col_ak < col_r) {
                            p_ak++;
                        } else {
                            p_r++;
                        }
                    }
                }
                col_ptr[ji + 1] = static_cast<int>(nz_row.size());
            }

            // 3. Normal Equations via Sparse Scatter
            double* B = B_vec.data();
            double* g = g_vec.data();

            for (int c1 = 0; c1 < n_ls; ++c1) {
                const int c1_start = col_ptr[c1];
                const int c1_end   = col_ptr[c1 + 1];

                double sum_g = 0.0;
                for (int k = c1_start; k < c1_end; ++k) {
                    const int r_idx = nz_row[k];
                    const double val = nz_val[k];
                    w_dense[r_idx] = val;
                    sum_g += val * f_dense[r_idx];
                }
                g[c1] = sum_g;

                // Write lower triangle of B only
                const int c1_n = c1 * n_ls;
                for (int c2 = 0; c2 <= c1; ++c2) {
                    const int c2_start = col_ptr[c2];
                    const int c2_end   = col_ptr[c2 + 1];

                    double sum_B = 0.0;
                    for (int k = c2_start; k < c2_end; ++k) {
                        sum_B += nz_val[k] * w_dense[nz_row[k]];
                    }
                    B[c1_n + c2] = sum_B;
                }

                // Unscatter Column c1
                for (int k = c1_start; k < c1_end; ++k) {
                    w_dense[nz_row[k]] = 0.0;
                }
            }

            // Clear touched entries of f_dense
            for (int idx : f_nz_idx) {
                f_dense[idx] = 0.0;
            }

            // 4. Solve System B * z = g
            const bool ok = fast_cholesky_solve_reg(n_ls, B, g, &out_val[s0],
                                                    L_vec.data(), y_vec.data());
            if (!ok) {
                std::fill_n(&out_val[s0], nnzi, 0.0);
            }
        }
    }
}
