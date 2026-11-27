#include "mex.h"
#include "matrix.h"
#include "BAMG.h"
#include "BAMG_params.h"
#include "DebEnv.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstdio>

// Returns [mat_I (sparse), c_mark] via zero-copy inputs, C++ TV transpose, and OpenMP CSR->CSC
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 19 || nlhs != 2) {
        mexErrMsgIdAndTxt("BAMG:badArgCount", "Expected 19 inputs and 2 outputs.");
    }

    const int    level      = static_cast<int>(mxGetScalar(prhs[0]));
    const int    np         = std::max(1, static_cast<int>(mxGetScalar(prhs[1])));
    const int    itmax_vol  = static_cast<int>(mxGetScalar(prhs[2]));
    const int    dist_min   = static_cast<int>(mxGetScalar(prhs[3]));
    const int    dist_max   = static_cast<int>(mxGetScalar(prhs[4]));
    const int    mmax       = static_cast<int>(mxGetScalar(prhs[5]));
    const double maxcond    = mxGetScalar(prhs[6]);
    const double maxrownrm  = mxGetScalar(prhs[7]);
    const double tol_vol    = mxGetScalar(prhs[8]);
    const double eps_prol   = mxGetScalar(prhs[9]);
    const int    nn_S       = static_cast<int>(mxGetScalar(prhs[10]));
    const int32_t *iat_S    = mxGetInt32s(prhs[11]);
    const int32_t *ja_S     = mxGetInt32s(prhs[12]);
    const int    ntv        = static_cast<int>(mxGetScalar(prhs[13]));
    const int32_t *fcnodes  = mxGetInt32s(prhs[14]);
    const double  *TV_colm  = mxGetDoubles(prhs[15]);
    const int    nn_I       = static_cast<int>(mxGetScalar(prhs[16]));
    const int    nc_I       = static_cast<int>(mxGetScalar(prhs[17]));
    const int    verb       = static_cast<int>(mxGetScalar(prhs[18]));

    // Transpose column-major TV [nn_S x ntv] into row-major TV pointers expected by BAMG
    std::vector<double> tv_flat(static_cast<std::size_t>(nn_S) * ntv);
    std::vector<double*> tv_ptr(nn_S);
    for (int i = 0; i < nn_S; ++i) {
        tv_ptr[i] = tv_flat.data() + static_cast<std::size_t>(i) * ntv;
    }
    for (int j = 0; j < ntv; ++j) {
        const double *col = TV_colm + static_cast<std::size_t>(j) * nn_S;
        for (int i = 0; i < nn_S; ++i) {
            tv_ptr[i][j] = col[i];
        }
    }

    if (level == 1) {
        DebEnv.SetDebEnv(np, "w");
    } else {
        DebEnv.OpenDebugLog("a");
    }
    if (DEBUG && BAMG_DEBUG) {
        std::fprintf(DebEnv.r_logfile, "\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n", level);
        std::fflush(DebEnv.r_logfile);
        for (int i = 0; i < np; ++i) {
            std::fprintf(DebEnv.t_logfile[i], "\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n", level);
            std::fflush(DebEnv.t_logfile[i]);
        }
    }

    BAMG_params params;
    params.verbosity  = verb;
    params.itmax_vol  = itmax_vol;
    params.dist_min   = dist_min;
    params.dist_max   = dist_max;
    params.mmax       = mmax;
    params.maxcond    = maxcond;
    params.maxrownrm  = maxrownrm;
    params.tol_vol    = tol_vol;
    params.eps        = eps_prol;

    int nt_I = 0;
    std::vector<int> vec_iat_I, vec_ja_I, vec_c_mark;
    std::vector<double> vec_coef_I;
    int ierr = BAMG(params, np, 0, nn_S, nn_S, iat_S, ja_S, ntv, fcnodes,
                    tv_ptr.empty() ? nullptr : const_cast<const double* const*>(tv_ptr.data()),
                    nt_I, vec_iat_I, vec_ja_I, vec_coef_I, vec_c_mark);
    DebEnv.CloseDebugLog();
    if (ierr != 0) {
        mexErrMsgIdAndTxt("BAMG:computeError", "BAMG failed with error code %d.", ierr);
    }

    // OpenMP CSR -> CSC into MATLAB sparse matrix (ascending row traversal guarantees sorted ir)
    plhs[0] = mxCreateSparse(static_cast<mwSize>(nn_I), static_cast<mwSize>(nc_I),
                             static_cast<mwSize>(std::max(1, nt_I)), mxREAL);
    mwIndex *jc = mxGetJc(plhs[0]);
    mwIndex *ir = mxGetIr(plhs[0]);
    double  *pr = mxGetDoubles(plhs[0]);

    if (nt_I > 0 && nc_I > 0 && nn_I > 0) {
        int nth_used = np;
        std::vector<int> thread_counts(static_cast<std::size_t>(np) * nc_I, 0);
        std::vector<int> thread_offsets(static_cast<std::size_t>(np) * nc_I, 0);
        std::vector<int> col_totals(nc_I, 0);

        #pragma omp parallel num_threads(np)
        {
            int tid = omp_get_thread_num();
            int nth = omp_get_num_threads();
            #pragma omp single
            nth_used = nth;

            int r0 = static_cast<int>((static_cast<long long>(tid) * nn_I) / nth);
            int r1 = static_cast<int>((static_cast<long long>(tid + 1) * nn_I) / nth);
            int *tc = &thread_counts[static_cast<std::size_t>(tid) * nc_I];
            for (int r = r0; r < r1; ++r) {
                for (int k = vec_iat_I[r]; k < vec_iat_I[r + 1]; ++k) {
                    tc[vec_ja_I[k]]++;
                }
            }
            #pragma omp barrier
            #pragma omp for schedule(static)
            for (int c = 0; c < nc_I; ++c) {
                int sum = 0;
                for (int t = 0; t < nth; ++t) {
                    thread_offsets[static_cast<std::size_t>(t) * nc_I + c] = sum;
                    sum += thread_counts[static_cast<std::size_t>(t) * nc_I + c];
                }
                col_totals[c] = sum;
            }
        }
        jc[0] = 0;
        for (int c = 0; c < nc_I; ++c) {
            jc[c + 1] = jc[c] + static_cast<mwIndex>(col_totals[c]);
        }

        #pragma omp parallel num_threads(nth_used)
        {
            int tid = omp_get_thread_num();
            int r0 = static_cast<int>((static_cast<long long>(tid) * nn_I) / nth_used);
            int r1 = static_cast<int>((static_cast<long long>(tid + 1) * nn_I) / nth_used);
            int *to = &thread_offsets[static_cast<std::size_t>(tid) * nc_I];
            for (int c = 0; c < nc_I; ++c) to[c] += static_cast<int>(jc[c]);
            for (int r = r0; r < r1; ++r) {
                for (int k = vec_iat_I[r]; k < vec_iat_I[r + 1]; ++k) {
                    int dst = to[vec_ja_I[k]]++;
                    ir[dst] = static_cast<mwIndex>(r);
                    pr[dst] = vec_coef_I[k];
                }
            }
        }
    } else {
        for (int c = 0; c <= nc_I; ++c) jc[c] = 0;
    }

    plhs[1] = mxCreateDoubleMatrix(static_cast<mwSize>(nn_I), 1, mxREAL);
    double *cm = mxGetDoubles(plhs[1]);
    for (int i = 0; i < nn_I; ++i) cm[i] = static_cast<double>(vec_c_mark[i]);
}
