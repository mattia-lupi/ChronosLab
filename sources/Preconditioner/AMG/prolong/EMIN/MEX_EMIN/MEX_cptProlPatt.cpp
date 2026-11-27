#include "mex.h"
#include "matrix.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>

namespace {

struct CandEntry {
    int col;
    double val;
};

struct CSRMatrix {
    int nrows = 0;
    int ncols = 0;
    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    std::vector<double> values;
};

// Convert MATLAB CSC sparse matrix into sorted-column CSR in a single OpenMP parallel region
CSRMatrix cscToCsr(const mxArray *mat, int np)
{
    CSRMatrix out;
    out.nrows = static_cast<int>(mxGetM(mat));
    out.ncols = static_cast<int>(mxGetN(mat));
    out.row_ptr.assign(static_cast<std::size_t>(out.nrows) + 1, 0);

    const mwIndex *jc = mxGetJc(mat);
    const mwIndex *ir = mxGetIr(mat);
    const double  *pr = mxGetDoubles(mat);
    const mwIndex nnz_total = (out.ncols > 0) ? jc[out.ncols] : 0;

    if (out.nrows == 0 || out.ncols == 0 || nnz_total == 0) {
        return out;
    }

    if (np <= 1 || nnz_total < 4096) {
        std::vector<int> counts(out.nrows, 0);
        for (mwIndex k = 0; k < nnz_total; ++k) {
            if (pr[k] != 0.0) counts[ir[k]]++;
        }
        for (int r = 0; r < out.nrows; ++r) {
            out.row_ptr[r + 1] = out.row_ptr[r] + counts[r];
        }
        const int nnz_nz = out.row_ptr[out.nrows];
        out.col_idx.resize(nnz_nz);
        out.values.resize(nnz_nz);
        std::vector<int> next_pos = out.row_ptr;
        for (int c = 0; c < out.ncols; ++c) {
            for (mwIndex k = jc[c]; k < jc[c + 1]; ++k) {
                const double v = pr[k];
                if (v != 0.0) {
                    const int dst = next_pos[ir[k]]++;
                    out.col_idx[dst] = c;
                    out.values[dst] = v;
                }
            }
        }
        return out;
    }

    std::vector<int> thread_counts(static_cast<std::size_t>(np) * out.nrows, 0);
    std::vector<int> thread_offsets(static_cast<std::size_t>(np) * out.nrows, 0);
    std::vector<int> row_totals(out.nrows, 0);

    #pragma omp parallel num_threads(np)
    {
        const int tid = omp_get_thread_num();
        const int nth = omp_get_num_threads();
        const int c0 = static_cast<int>((static_cast<long long>(tid) * out.ncols) / nth);
        const int c1 = static_cast<int>((static_cast<long long>(tid + 1) * out.ncols) / nth);
        int *tc = &thread_counts[static_cast<std::size_t>(tid) * out.nrows];

        for (int c = c0; c < c1; ++c) {
            for (mwIndex k = jc[c]; k < jc[c + 1]; ++k) {
                if (pr[k] != 0.0) tc[ir[k]]++;
            }
        }
        #pragma omp barrier

        #pragma omp for schedule(static)
        for (int r = 0; r < out.nrows; ++r) {
            int sum = 0;
            for (int t = 0; t < nth; ++t) {
                const std::size_t idx = static_cast<std::size_t>(t) * out.nrows + r;
                thread_offsets[idx] = sum;
                sum += thread_counts[idx];
            }
            row_totals[r] = sum;
        }

        #pragma omp single
        {
            for (int r = 0; r < out.nrows; ++r) {
                out.row_ptr[r + 1] = out.row_ptr[r] + row_totals[r];
            }
            const int nnz_nz = out.row_ptr[out.nrows];
            out.col_idx.resize(nnz_nz);
            out.values.resize(nnz_nz);
        }

        int *to = &thread_offsets[static_cast<std::size_t>(tid) * out.nrows];
        for (int r = 0; r < out.nrows; ++r) to[r] += out.row_ptr[r];

        for (int c = c0; c < c1; ++c) {
            for (mwIndex k = jc[c]; k < jc[c + 1]; ++k) {
                const double v = pr[k];
                if (v != 0.0) {
                    const int dst = to[ir[k]]++;
                    out.col_idx[dst] = c;
                    out.values[dst] = v;
                }
            }
        }
    }
    return out;
}

mxArray *arenasToMxSparse(int nrows, int ncols, int nt,
                          const std::vector<int> &row_tid,
                          const std::vector<int> &row_off,
                          const std::vector<int> &row_len,
                          const std::vector<std::vector<int>> &thr_ja,
                          const std::vector<std::vector<double>> &thr_val,
                          int np)
{
    mxArray *res = mxCreateSparse(static_cast<mwSize>(nrows), static_cast<mwSize>(ncols),
                                  static_cast<mwSize>(std::max(1, nt)), mxREAL);
    mwIndex *jc = mxGetJc(res);
    mwIndex *ir = mxGetIr(res);
    double  *pr = mxGetDoubles(res);

    if (nt <= 0 || ncols <= 0 || nrows <= 0) {
        for (int c = 0; c <= ncols; ++c) jc[c] = 0;
        return res;
    }

    if (np <= 1 || nt < 4096) {
        std::vector<int> col_counts(ncols, 0);
        for (int r = 0; r < nrows; ++r) {
            const int len = row_len[r];
            if (len == 0) continue;
            const int *ja = thr_ja[row_tid[r]].data() + row_off[r];
            for (int k = 0; k < len; ++k) col_counts[ja[k]]++;
        }
        jc[0] = 0;
        for (int c = 0; c < ncols; ++c) jc[c + 1] = jc[c] + static_cast<mwIndex>(col_counts[c]);
        std::vector<int> next_pos(ncols);
        for (int c = 0; c < ncols; ++c) next_pos[c] = static_cast<int>(jc[c]);
        for (int r = 0; r < nrows; ++r) {
            const int len = row_len[r];
            if (len == 0) continue;
            const int tid = row_tid[r];
            const int off = row_off[r];
            const int    *ja  = thr_ja[tid].data()  + off;
            const double *val = thr_val[tid].data() + off;
            for (int k = 0; k < len; ++k) {
                const int dst = next_pos[ja[k]]++;
                ir[dst] = static_cast<mwIndex>(r);
                pr[dst] = val[k];
            }
        }
        return res;
    }

    std::vector<int> thread_counts(static_cast<std::size_t>(np) * ncols, 0);
    std::vector<int> thread_offsets(static_cast<std::size_t>(np) * ncols, 0);
    std::vector<int> col_totals(ncols, 0);

    #pragma omp parallel num_threads(np)
    {
        const int tid = omp_get_thread_num();
        const int nth = omp_get_num_threads();
        const int r0 = static_cast<int>((static_cast<long long>(tid) * nrows) / nth);
        const int r1 = static_cast<int>((static_cast<long long>(tid + 1) * nrows) / nth);
        int *tc = &thread_counts[static_cast<std::size_t>(tid) * ncols];

        for (int r = r0; r < r1; ++r) {
            const int len = row_len[r];
            if (len == 0) continue;
            const int *ja = thr_ja[row_tid[r]].data() + row_off[r];
            for (int k = 0; k < len; ++k) tc[ja[k]]++;
        }
        #pragma omp barrier

        #pragma omp for schedule(static)
        for (int c = 0; c < ncols; ++c) {
            int sum = 0;
            for (int t = 0; t < nth; ++t) {
                const std::size_t idx = static_cast<std::size_t>(t) * ncols + c;
                thread_offsets[idx] = sum;
                sum += thread_counts[idx];
            }
            col_totals[c] = sum;
        }

        #pragma omp single
        {
            jc[0] = 0;
            for (int c = 0; c < ncols; ++c) {
                jc[c + 1] = jc[c] + static_cast<mwIndex>(col_totals[c]);
            }
        }

        int *to = &thread_offsets[static_cast<std::size_t>(tid) * ncols];
        for (int c = 0; c < ncols; ++c) to[c] += static_cast<int>(jc[c]);

        // Sequential ascending row traversal inside each stripe guarantees sorted ir per CSC column
        for (int r = r0; r < r1; ++r) {
            const int len = row_len[r];
            if (len == 0) continue;
            const int rt = row_tid[r];
            const int ro = row_off[r];
            const int    *ja  = thr_ja[rt].data()  + ro;
            const double *val = thr_val[rt].data() + ro;
            for (int k = 0; k < len; ++k) {
                const int dst = to[ja[k]]++;
                ir[dst] = static_cast<mwIndex>(r);
                pr[dst] = val[k];
            }
        }
    }

    return res;
}

} // namespace

// P_patt = MEX_cptProlPatt(kpow, S, P, fcnode, nnzr, nnzr_min, np, [verb])
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ((nrhs != 7 && nrhs != 8) || nlhs != 1) {
        mexErrMsgIdAndTxt("MEX_cptProlPatt:badArgCount",
                          "Expected 7 or 8 inputs and 1 output.");
    }
    if (!mxIsSparse(prhs[1]) || !mxIsSparse(prhs[2])) {
        mexErrMsgIdAndTxt("MEX_cptProlPatt:badInputType",
                          "S and P must be sparse double matrices.");
    }

    const int kpow     = static_cast<int>(mxGetScalar(prhs[0]));
    const mxArray *S_mx = prhs[1];
    const mxArray *P_mx = prhs[2];
    const double *fcnode = mxGetDoubles(prhs[3]);
    const int nnzr_raw = static_cast<int>(mxGetScalar(prhs[4]));
    const int nnzr_min = static_cast<int>(mxGetScalar(prhs[5]));
    const int np       = std::max(1, static_cast<int>(mxGetScalar(prhs[6])));
    const int verb     = (nrhs >= 8) ? static_cast<int>(mxGetScalar(prhs[7])) : 0;

    const int nn = static_cast<int>(mxGetM(P_mx));
    const int nc = static_cast<int>(mxGetN(P_mx));
    const int nnzr = std::max(nnzr_raw, nnzr_min);

    if (static_cast<int>(mxGetM(S_mx)) != nn || static_cast<int>(mxGetN(S_mx)) != nn ||
        static_cast<int>(mxGetNumberOfElements(prhs[3])) != nn) {
        mexErrMsgIdAndTxt("MEX_cptProlPatt:dimMismatch",
                          "Dimensions of S, P, and fcnode must match.");
    }

    if (nn == 0 || nc == 0) {
        plhs[0] = mxCreateSparse(static_cast<mwSize>(nn), static_cast<mwSize>(nc), 1, mxREAL);
        mwIndex *jc = mxGetJc(plhs[0]);
        for (int c = 0; c <= nc; ++c) jc[c] = 0;
        return;
    }

    CSRMatrix P0 = cscToCsr(P_mx, np);
    const double nf_d = static_cast<double>(std::max(1, nn - nc));

    if (verb) {
        mexPrintf("Power 1 nnzr: %10.2f\n", static_cast<double>(P0.row_ptr[nn] - nc) / nf_d);
    }

    // Handle kpow <= 1: retain original P0 values on fine rows only
    if (kpow <= 1) {
        std::vector<int> row_tid(nn, 0), row_off(nn, 0), row_len(nn, 0);
        int nt_out = 0;
        for (int i = 0; i < nn; ++i) {
            if (fcnode[i] <= 0.0) {
                row_off[i] = P0.row_ptr[i];
                row_len[i] = P0.row_ptr[i + 1] - P0.row_ptr[i];
                nt_out += row_len[i];
            }
        }
        std::vector<std::vector<int>> thr_ja(1);
        std::vector<std::vector<double>> thr_val(1);
        thr_ja[0] = std::move(P0.col_idx);
        thr_val[0] = std::move(P0.values);
        plhs[0] = arenasToMxSparse(nn, nc, nt_out, row_tid, row_off, row_len, thr_ja, thr_val, np);
        return;
    }

    CSRMatrix S = cscToCsr(S_mx, np);
    const CSRMatrix *P_cur_ptr = &P0;
    CSRMatrix P_curr;

    auto cmp_abs_desc = [](const CandEntry &a, const CandEntry &b) {
        const double aa = std::fabs(a.val);
        const double bb = std::fabs(b.val);
        if (aa != bb) return aa > bb;
        return a.col < b.col;
    };
    auto cmp_col_asc = [](const CandEntry &a, const CandEntry &b) {
        return a.col < b.col;
    };

    // Persistent thread-local output arenas reused across ipow iterations (zero per-chunk allocs)
    std::vector<std::vector<int>> thr_ja(np);
    std::vector<std::vector<double>> thr_val(np);
    std::vector<int> row_tid(nn, 0);
    std::vector<int> row_off(nn, 0);
    std::vector<int> row_len(nn, 0);
    const int est_per_thr = std::max(1024, (nn / np + 64) * std::max(4, nnzr));

    int nt_last = 0;
    for (int ipow = 2; ipow <= kpow; ++ipow) {
        const bool is_last = (ipow == kpow);
        nt_last = 0;

        #pragma omp parallel num_threads(np) reduction(+:nt_last)
        {
            const int tid = omp_get_thread_num();
            // Move into stack-local vectors so push_back never false-shares std::vector headers
            std::vector<int> out_ja;
            std::vector<double> out_val;
            out_ja.swap(thr_ja[tid]);
            out_val.swap(thr_val[tid]);
            out_ja.clear();
            out_val.clear();
            if (static_cast<int>(out_ja.capacity()) < est_per_thr) {
                out_ja.reserve(est_per_thr);
                out_val.reserve(est_per_thr);
            }

            std::vector<int> marker(nc, 0);
            std::vector<double> accum(nc, 0.0);
            std::vector<int> cand_cols;
            std::vector<CandEntry> cand_extra;
            cand_cols.reserve(256);
            cand_extra.reserve(256);
            int row_iter = 0;

            const int *__restrict s_ptr  = S.row_ptr.data();
            const int *__restrict s_col  = S.col_idx.data();
            const double *__restrict s_val = S.values.data();
            const int *__restrict p0_ptr = P0.row_ptr.data();
            const int *__restrict p0_col = P0.col_idx.data();
            const int *__restrict pc_ptr = P_cur_ptr->row_ptr.data();
            const int *__restrict pc_col = P_cur_ptr->col_idx.data();
            const double *__restrict pc_val = P_cur_ptr->values.data();

            #pragma omp for schedule(dynamic, 64)
            for (int i = 0; i < nn; ++i) {
                if (is_last && fcnode[i] > 0.0) {
                    row_len[i] = 0;
                    continue;
                }

                const int p0_s = p0_ptr[i];
                const int p0_e = p0_ptr[i + 1];
                const int n0 = p0_e - p0_s;
                const int budget = std::max(0, nnzr - n0);
                const int off = static_cast<int>(out_ja.size());
                row_tid[i] = tid;
                row_off[i] = off;

                // When P0 already fills nnzr, drop_patt output is identically mask_0 (already sorted)
                if (budget == 0) {
                    for (int k = p0_s; k < p0_e; ++k) {
                        out_ja.push_back(p0_col[k]);
                        out_val.push_back(1.0);
                    }
                    row_len[i] = n0;
                    nt_last += n0;
                    continue;
                }

                if (row_iter > 1000000000) {
                    std::fill(marker.begin(), marker.end(), 0);
                    row_iter = 0;
                }
                const int stamp_p0   = 2 * (++row_iter);
                const int stamp_seen = stamp_p0 + 1;

                for (int k = p0_s; k < p0_e; ++k) {
                    marker[p0_col[k]] = stamp_p0;
                }

                cand_cols.clear();

                const int s_s = s_ptr[i];
                const int s_e = s_ptr[i + 1];
                for (int ks = s_s; ks < s_e; ++ks) {
                    const int j = s_col[ks];
                    const double s_ij = s_val[ks];
                    const int p_s = pc_ptr[j];
                    const int p_e = pc_ptr[j + 1];
                    for (int kp = p_s; kp < p_e; ++kp) {
                        const int col = pc_col[kp];
                        const int m = marker[col];
                        if (m == stamp_seen) {
                            accum[col] += s_ij * pc_val[kp];
                        } else if (m < stamp_p0) {
                            marker[col] = stamp_seen;
                            cand_cols.push_back(col);
                            accum[col] = s_ij * pc_val[kp];
                        }
                    }
                }

                cand_extra.clear();
                for (int col : cand_cols) {
                    const double v = accum[col];
                    if (v != 0.0) {
                        cand_extra.push_back({col, v});
                    }
                }

                if (cand_extra.size() > static_cast<std::size_t>(budget)) {
                    std::nth_element(cand_extra.begin(),
                                     cand_extra.begin() + budget,
                                     cand_extra.end(),
                                     cmp_abs_desc);
                    cand_extra.resize(static_cast<std::size_t>(budget));
                }

                for (int k = p0_s; k < p0_e; ++k) {
                    cand_extra.push_back({p0_col[k], 1.0});
                }
                std::sort(cand_extra.begin(), cand_extra.end(), cmp_col_asc);

                for (const CandEntry &ce : cand_extra) {
                    out_ja.push_back(ce.col);
                    out_val.push_back(ce.val);
                }
                const int rlen = static_cast<int>(cand_extra.size());
                row_len[i] = rlen;
                nt_last += rlen;
            }
            thr_ja[tid].swap(out_ja);
            thr_val[tid].swap(out_val);
        }

        // Only flatten arenas into contiguous CSR if another power iteration follows
        if (!is_last) {
            P_curr.nrows = nn;
            P_curr.ncols = nc;
            P_curr.row_ptr.resize(static_cast<std::size_t>(nn) + 1);
            P_curr.row_ptr[0] = 0;
            for (int i = 0; i < nn; ++i) P_curr.row_ptr[i + 1] = P_curr.row_ptr[i] + row_len[i];
            P_curr.col_idx.resize(nt_last);
            P_curr.values.resize(nt_last);

            #pragma omp parallel for schedule(static) num_threads(np)
            for (int i = 0; i < nn; ++i) {
                const int len = row_len[i];
                if (len == 0) continue;
                const int *src_j = thr_ja[row_tid[i]].data()  + row_off[i];
                const double *src_v = thr_val[row_tid[i]].data() + row_off[i];
                const int dst = P_curr.row_ptr[i];
                for (int k = 0; k < len; ++k) {
                    P_curr.col_idx[dst + k] = src_j[k];
                    P_curr.values[dst + k]  = src_v[k];
                }
            }
            P_cur_ptr = &P_curr;
            if (verb) {
                mexPrintf("Power %d nnzr: %10.2f\n", ipow, static_cast<double>(nt_last) / nf_d);
            }
        }
    }

    if (verb) {
        int nn_poor = 0;
        for (int i = 0; i < nn; ++i) {
            if (fcnode[i] <= 0.0) {
                const int cnt = row_len[i];
                if (cnt > 0 && cnt < nnzr_min) nn_poor++;
            }
        }
        mexPrintf("Final nnzr: %10.2f\n", static_cast<double>(nt_last - nc) / nf_d);
        mexPrintf("There are %d poorly connected nodes\n", nn_poor);
    }

    plhs[0] = arenasToMxSparse(nn, nc, nt_last, row_tid, row_off, row_len, thr_ja, thr_val, np);
}
