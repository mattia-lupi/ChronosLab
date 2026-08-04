// sam_compute_left_mex.cpp
#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
#include <omp.h>

using namespace matlab::data;
using namespace matlab::mex;

extern "C" {
void dgels_(const char* trans, const ptrdiff_t* m, const ptrdiff_t* n, const ptrdiff_t* nrhs,
            double* a, const ptrdiff_t* lda, double* b, const ptrdiff_t* ldb,
            double* work, const ptrdiff_t* lwork, ptrdiff_t* info);

void dposv_(const char* uplo, const ptrdiff_t* n, const ptrdiff_t* nrhs,
            double* a, const ptrdiff_t* lda, double* b, const ptrdiff_t* ldb, ptrdiff_t* info);

void dsyrk_(const char* uplo, const char* trans, const ptrdiff_t* n, const ptrdiff_t* k,
            const double* alpha, const double* a, const ptrdiff_t* lda,
            const double* beta, double* c, const ptrdiff_t* ldc);

void dgemv_(const char* trans, const ptrdiff_t* m, const ptrdiff_t* n,
            const double* alpha, const double* a, const ptrdiff_t* lda,
            const double* x, const ptrdiff_t* incx,
            const double* beta, double* y, const ptrdiff_t* incy);
}

struct CSR {
    int n;
    std::vector<int>    ptr;
    std::vector<int>    col;
    std::vector<double> val;
};

static CSR buildCSR(const TypedArray<double>& iat, const TypedArray<double>& ja, const TypedArray<double>& coef) {
    CSR csr;
    csr.n   = static_cast<int>(iat.getNumberOfElements()) - 1;
    int nnz = static_cast<int>(ja.getNumberOfElements());
    csr.ptr.resize(csr.n + 1);
    csr.col.resize(nnz);
    csr.val.resize(nnz);
    for (int i = 0; i <= csr.n; ++i) csr.ptr[i] = static_cast<int>(iat[i]) - 1;
    for (int p = 0; p < nnz; ++p) {
        csr.col[p] = static_cast<int>(ja[p]) - 1;
        csr.val[p] = coef[p];
    }
    return csr;
}

static int queryLwork(int max_m, int max_n) {
    if (max_m <= 0 || max_n <= 0) return 128;
    std::vector<double> Atmp(max_m * max_n, 0.0);
    int mn = std::min(max_m, max_n);
    for (int k = 0; k < mn; ++k) Atmp[k * max_m + k] = 1.0;

    char trans = 'N';
    ptrdiff_t pm = max_m, pn = max_n, pnrhs = 1;
    ptrdiff_t plda = pm, pldb = std::max(pm, pn), plw = -1, info = 0;
    double wq = 0.0;
    std::vector<double> ftmp(std::max(max_m, max_n), 0.0);

    dgels_(&trans, &pm, &pn, &pnrhs, Atmp.data(), &plda, ftmp.data(), &pldb, &wq, &plw, &info);
    return std::max(128, static_cast<int>(wq));
}

inline bool solve_micro_kernel(int n, const double* B, const double* g, double* z) {
    if (n == 1) {
        if (B[0] <= 0.0) return false;
        z[0] = g[0] / B[0];
        return true;
    } else if (n == 2) {
        double det = B[0]*B[3] - B[1]*B[2];
        if (det <= 0.0) return false;
        z[0] = (B[3]*g[0] - B[2]*g[1]) / det;
        z[1] = (B[0]*g[1] - B[1]*g[0]) / det;
        return true;
    } else if (n == 3) {
        double det = B[0]*(B[4]*B[8] - B[5]*B[7]) - B[3]*(B[1]*B[8] - B[2]*B[7]) + B[6]*(B[1]*B[5] - B[2]*B[4]);
        if (det <= 0.0) return false;
        double invDet = 1.0 / det;
        z[0] = ((B[4]*B[8] - B[5]*B[7])*g[0] + (B[6]*B[5] - B[3]*B[8])*g[1] + (B[3]*B[7] - B[6]*B[4])*g[2]) * invDet;
        z[1] = ((B[7]*B[2] - B[1]*B[8])*g[0] + (B[0]*B[8] - B[6]*B[2])*g[1] + (B[1]*B[6] - B[0]*B[7])*g[2]) * invDet;
        z[2] = ((B[1]*B[5] - B[4]*B[2])*g[0] + (B[3]*B[2] - B[0]*B[5])*g[1] + (B[0]*B[4] - B[3]*B[1])*g[2]) * invDet;
        return true;
    }
    return false;
}

class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        ArrayFactory factory;

        if (inputs.size() != 11) {
            getEngine()->feval(u"error", 0, std::vector<Array>({factory.createScalar("sam_compute_left_mex: requires iatk, jak, coefk, iat0, ja0, coef0, s_ptr, s_data, r_ptr, r_data, nnz_total.")}));
            return;
        }

        TypedArray<double> iatk  = std::move(inputs[0]); TypedArray<double> jak   = std::move(inputs[1]); TypedArray<double> coefk = std::move(inputs[2]);
        TypedArray<double> iat0  = std::move(inputs[3]); TypedArray<double> ja0   = std::move(inputs[4]); TypedArray<double> coef0 = std::move(inputs[5]);
        
        TypedArray<int32_t> s_ptr_arr = std::move(inputs[6]);
        TypedArray<int32_t> s_dat_arr = std::move(inputs[7]);
        TypedArray<int32_t> r_ptr_arr = std::move(inputs[8]);
        TypedArray<int32_t> r_dat_arr = std::move(inputs[9]);
        
        int total_nnz;
        if (inputs[10].getType() == ArrayType::INT32) {
            total_nnz = static_cast<TypedArray<int32_t>>(inputs[10])[0];
        } else {
            total_nnz = static_cast<int>(static_cast<TypedArray<double>>(inputs[10])[0]);
        }

        const CSR Ak = buildCSR(iatk, jak, coefk);
        const CSR A0 = buildCSR(iat0, ja0, coef0);
        const int n = Ak.n;

        // Zero-copy pointers to MATLAB memory
        const int32_t* s_ptr = &(*s_ptr_arr.begin());
        const int32_t* s_data = &(*s_dat_arr.begin());
        const int32_t* r_ptr = &(*r_ptr_arr.begin());
        const int32_t* r_data = &(*r_dat_arr.begin());

        int max_si = 0, max_ri = 0;
        for (int i = 0; i < n; ++i) {
            max_si = std::max(max_si, s_ptr[i+1] - s_ptr[i]);
            max_ri = std::max(max_ri, r_ptr[i+1] - r_ptr[i]);
        }
        const int LWORK = queryLwork(max_ri, max_si);
        const int max_ldb = std::max(max_ri, max_si);

        const size_t sz = static_cast<size_t>(total_nnz);
        TypedArray<double> out_row = factory.createArray<double>({sz, 1});
        TypedArray<double> out_col = factory.createArray<double>({sz, 1});
        TypedArray<double> out_val = factory.createArray<double>({sz, 1});

        double* p_row = &(*out_row.begin()); double* p_col = &(*out_col.begin()); double* p_val = &(*out_val.begin());

        for (int i = 0; i < n; ++i) {
            for (int p = s_ptr[i]; p < s_ptr[i+1]; ++p) {
                p_row[p] = static_cast<double>(i + 1);
                p_col[p] = static_cast<double>(s_data[p] + 1);
            }
        }

        #pragma omp parallel
        {
            std::vector<double> Atmp(max_ri * max_si);
            std::vector<double> f(max_ldb);
            std::vector<double> work(LWORK);
            std::vector<double> B(64 * 64);
            std::vector<double> g(64);

            #pragma omp for schedule(guided) nowait
            for (int i = 0; i < n; ++i) {
                const int s0 = s_ptr[i], s1 = s_ptr[i + 1], nnzi = s1 - s0;
                if (nnzi == 0) continue;

                const int r0 = r_ptr[i], r1 = r_ptr[i + 1], nri = r1 - r0;
                const int32_t* si = s_data + s0;
                const int32_t* ri = r_data + r0;

                for (int p = s0; p < s1; ++p) p_val[p] = 0.0;
                if (nri == 0) continue;

                const int m_ls = nri, n_ls = nnzi;
                const int ldb = std::max(m_ls, n_ls);
                bool solved = false;

                if (n_ls <= 16 && m_ls >= n_ls) {
                    for (int c1 = 0; c1 < n_ls; ++c1) {
                        int row1 = si[c1];
                        double sum_g = 0.0;
                        int pa = Ak.ptr[row1], ea = Ak.ptr[row1+1];
                        int p0 = A0.ptr[i], e0 = A0.ptr[i+1];
                        while(pa < ea && p0 < e0) {
                            if (Ak.col[pa] == A0.col[p0]) { sum_g += Ak.val[pa]*A0.val[p0]; pa++; p0++; }
                            else if (Ak.col[pa] < A0.col[p0]) pa++;
                            else p0++;
                        }
                        g[c1] = sum_g;

                        for (int c2 = 0; c2 <= c1; ++c2) {
                            int row2 = si[c2];
                            double sum_B = 0.0;
                            int p1 = Ak.ptr[row1], e1 = Ak.ptr[row1+1];
                            int p2 = Ak.ptr[row2], e2 = Ak.ptr[row2+1];
                            while(p1 < e1 && p2 < e2) {
                                if (Ak.col[p1] == Ak.col[p2]) { sum_B += Ak.val[p1]*Ak.val[p2]; p1++; p2++; }
                                else if (Ak.col[p1] < Ak.col[p2]) p1++;
                                else p2++;
                            }
                            B[c1 * n_ls + c2] = sum_B;
                            B[c2 * n_ls + c1] = sum_B;
                        }
                    }

                    if (n_ls <= 3) {
                        solved = solve_micro_kernel(n_ls, B.data(), g.data(), &p_val[s0]);
                    }

                    if (!solved) {
                        char uplo = 'U'; ptrdiff_t pn_chol = n_ls, pnrhs_chol = 1, plda_chol = n_ls, pldb_chol = n_ls, info_chol = 0;
                        dposv_(&uplo, &pn_chol, &pnrhs_chol, B.data(), &plda_chol, g.data(), &pldb_chol, &info_chol);
                        if (info_chol == 0) {
                            for (int ji = 0; ji < n_ls; ++ji) p_val[s0 + ji] = std::isfinite(g[ji]) ? g[ji] : 0.0;
                            solved = true;
                        }
                    }
                }

                if (!solved) {
                    std::fill(Atmp.begin(), Atmp.begin() + (m_ls * n_ls), 0.0);
                    std::fill(f.begin(), f.begin() + ldb, 0.0);

                    for (int ji = 0; ji < nnzi; ++ji) {
                        const int row_k = si[ji], ks = Ak.ptr[row_k], ke = Ak.ptr[row_k + 1];
                        int p_ak = ks, p_ri = 0;
                        while (p_ak < ke && p_ri < nri) {
                            if (Ak.col[p_ak] == ri[p_ri]) { Atmp[ji * m_ls + p_ri] = Ak.val[p_ak]; p_ak++; p_ri++; }
                            else if (Ak.col[p_ak] < ri[p_ri]) p_ak++;
                            else p_ri++;
                        }
                    }

                    {
                        const int as = A0.ptr[i], ae = A0.ptr[i + 1];
                        int p_a0 = as, p_ri = 0;
                        while (p_a0 < ae && p_ri < nri) {
                            if (A0.col[p_a0] == ri[p_ri]) { f[p_ri] = A0.val[p_a0]; p_a0++; p_ri++; }
                            else if (A0.col[p_a0] < ri[p_ri]) p_a0++;
                            else p_ri++;
                        }
                    }

                    if (n_ls <= 64 && m_ls >= n_ls) {
                        char uplo = 'U', trans = 'T';
                        ptrdiff_t pn_blas = n_ls, pk_blas = m_ls;
                        double alpha = 1.0, beta = 0.0;
                        ptrdiff_t plda_blas = m_ls, pldc_blas = n_ls;
                        dsyrk_(&uplo, &trans, &pn_blas, &pk_blas, &alpha, Atmp.data(), &plda_blas, &beta, B.data(), &pldc_blas);

                        char trans_v = 'T'; ptrdiff_t incx = 1, incy = 1;
                        dgemv_(&trans_v, &pk_blas, &pn_blas, &alpha, Atmp.data(), &plda_blas, f.data(), &incx, &beta, g.data(), &incy);

                        char uplo_chol = 'U'; ptrdiff_t pnrhs_chol = 1, info_chol = 0;
                        dposv_(&uplo_chol, &pn_blas, &pnrhs_chol, B.data(), &pldc_blas, g.data(), &pldc_blas, &info_chol);

                        if (info_chol == 0) {
                            for (int ji = 0; ji < n_ls; ++ji) p_val[s0 + ji] = std::isfinite(g[ji]) ? g[ji] : 0.0;
                            solved = true;
                        }
                    }

                    if (!solved) {
                        char trans = 'N'; ptrdiff_t pm = m_ls, pn = n_ls, pnrhs = 1, plda = m_ls, pldb = ldb, plw = static_cast<ptrdiff_t>(work.size()), info = 0;
                        dgels_(&trans, &pm, &pn, &pnrhs, Atmp.data(), &plda, f.data(), &pldb, work.data(), &plw, &info);
                        for (int ji = 0; ji < n_ls; ++ji) p_val[s0 + ji] = (info == 0 && std::isfinite(f[ji])) ? f[ji] : 0.0;
                    }
                }
            }
        }

        outputs[0] = std::move(out_row); outputs[1] = std::move(out_col); outputs[2] = std::move(out_val);
    }
};
