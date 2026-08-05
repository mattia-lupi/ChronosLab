#include "sam_compute_left.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>

extern "C" {
void dposv_(const char* uplo, const ptrdiff_t* n, const ptrdiff_t* nrhs,
            double* a, const ptrdiff_t* lda, double* b, const ptrdiff_t* ldb, ptrdiff_t* info);

void dgels_(const char* trans, const ptrdiff_t* m, const ptrdiff_t* n, const ptrdiff_t* nrhs,
            double* a, const ptrdiff_t* lda, double* b, const ptrdiff_t* ldb,
            double* work, const ptrdiff_t* lwork, ptrdiff_t* info);
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
        std::vector<double> B_vec;
        std::vector<double> g_vec;
        std::vector<double> Atmp_vec;
        std::vector<double> f_vec;
        std::vector<double> work_vec(512);

        #pragma omp for schedule(guided) nowait
        for (int i = 0; i < n; ++i) {
            const int s0 = s_ptr[i], s1 = s_ptr[i + 1], nnzi = s1 - s0;
            if (nnzi == 0) continue;

            const int r0 = r_ptr[i], r1 = r_ptr[i + 1], nri = r1 - r0;
            const int* si = s_data + s0;
            const int* ri = r_data + r0;

            std::fill_n(&out_val[s0], nnzi, 0.0);
            if (nri == 0) continue;

            const int m_ls = nri, n_ls = nnzi;
            const int max_dim = std::max(n_ls, m_ls);
            bool solved = false;

            if (static_cast<int>(g_vec.size()) < max_dim) {
                B_vec.resize(max_dim * max_dim);
                g_vec.resize(max_dim);
            }

            // Sparse Assembly Speedup Path (Construct B and g via sparse intersection)
            if (n_ls <= 16 && m_ls >= n_ls) {
                double* B = B_vec.data();
                double* g = g_vec.data();

                for (int c1 = 0; c1 < n_ls; ++c1) {
                    const int row1 = si[c1];
                    double sum_g = 0.0;
                    int pa = ak_ptr[row1], ea = ak_ptr[row1 + 1];
                    int p0 = a0_ptr[i],    e0 = a0_ptr[i + 1];
                    while (pa < ea && p0 < e0) {
                        if (ak_col[pa] == a0_col[p0]) {
                            sum_g += ak_val[pa] * a0_val[p0];
                            pa++; p0++;
                        } else if (ak_col[pa] < a0_col[p0]) {
                            pa++;
                        } else {
                            p0++;
                        }
                    }
                    g[c1] = sum_g;

                    for (int c2 = 0; c2 <= c1; ++c2) {
                        const int row2 = si[c2];
                        double sum_B = 0.0;
                        int p1 = ak_ptr[row1], e1 = ak_ptr[row1 + 1];
                        int p2 = ak_ptr[row2], e2 = ak_ptr[row2 + 1];
                        while (p1 < e1 && p2 < e2) {
                            if (ak_col[p1] == ak_col[p2]) {
                                sum_B += ak_val[p1] * ak_val[p2];
                                p1++; p2++;
                            } else if (ak_col[p1] < ak_col[p2]) {
                                p1++;
                            } else {
                                p2++;
                            }
                        }
                        B[c1 * n_ls + c2] = sum_B;
                        B[c2 * n_ls + c1] = sum_B;
                    }
                }

                // LAPACK Accuracy: Call dposv_ directly on assembled system
                char uplo = 'U';
                ptrdiff_t pn = n_ls, pnrhs = 1, plda = n_ls, pldb = n_ls, info = 0;
                dposv_(&uplo, &pn, &pnrhs, B, &plda, g, &pldb, &info);

                if (info == 0) {
                    for (int ji = 0; ji < n_ls; ++ji) {
                        out_val[s0 + ji] = std::isfinite(g[ji]) ? g[ji] : 0.0;
                    }
                    solved = true;
                }
            }

            // Fallback Extraction Path: Use LAPACK dgels_ for ill-conditioned or non-square systems
            if (!solved) {
                const int ldb = std::max(m_ls, n_ls);
                if (static_cast<int>(Atmp_vec.size()) < m_ls * n_ls) {
                    Atmp_vec.resize(m_ls * n_ls);
                }
                if (static_cast<int>(f_vec.size()) < ldb) {
                    f_vec.resize(ldb);
                }

                double* Atmp = Atmp_vec.data();
                double* f    = f_vec.data();

                std::fill_n(Atmp, m_ls * n_ls, 0.0);
                std::fill_n(f, ldb, 0.0);

                for (int ji = 0; ji < nnzi; ++ji) {
                    const int row_k = si[ji];
                    const int ks = ak_ptr[row_k], ke = ak_ptr[row_k + 1];
                    int p_ak = ks, p_ri = 0;
                    while (p_ak < ke && p_ri < nri) {
                        if (ak_col[p_ak] == ri[p_ri]) {
                            Atmp[ji * m_ls + p_ri] = ak_val[p_ak];
                            p_ak++; p_ri++;
                        } else if (ak_col[p_ak] < ri[p_ri]) {
                            p_ak++;
                        } else {
                            p_ri++;
                        }
                    }
                }

                {
                    const int as = a0_ptr[i], ae = a0_ptr[i + 1];
                    int p_a0 = as, p_ri = 0;
                    while (p_a0 < ae && p_ri < nri) {
                        if (a0_col[p_a0] == ri[p_ri]) {
                            f[p_ri] = a0_val[p_a0];
                            p_a0++; p_ri++;
                        } else if (a0_col[p_a0] < ri[p_ri]) {
                            p_a0++;
                        } else {
                            p_ri++;
                        }
                    }
                }

                char trans = 'N';
                ptrdiff_t pm = m_ls, pn = n_ls, pnrhs = 1, plda = m_ls, pldb = ldb;
                ptrdiff_t plw = static_cast<ptrdiff_t>(work_vec.size()), info = 0;

                dgels_(&trans, &pm, &pn, &pnrhs, Atmp, &plda, f, &pldb, work_vec.data(), &plw, &info);

                if (info == 0) {
                    for (int ji = 0; ji < n_ls; ++ji) {
                        out_val[s0 + ji] = std::isfinite(f[ji]) ? f[ji] : 0.0;
                    }
                }
            }
        }
    }
}
