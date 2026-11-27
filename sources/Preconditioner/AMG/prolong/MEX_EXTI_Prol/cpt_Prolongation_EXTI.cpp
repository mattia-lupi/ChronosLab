#include "mex.h"
#include "matrix.h"
#include "EXTI_prolongation.h"
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cstdint>

// Returns MATLAB sparse matrix directly via C++ vecstart and OpenMP CSR->CSC assembly
// Inputs (11): level, np, nn_A, nt_A, iat_A, ja_A, coef_A, coef_S, fcnodes, nn_I, nc_I
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 11 || nlhs != 1) {
        mexErrMsgIdAndTxt("EXTI_Prol:badArgCount", "Expected 11 inputs and 1 output.");
    }

    const int level    = static_cast<int>(mxGetScalar(prhs[0]));
    const int np       = std::max(1, static_cast<int>(mxGetScalar(prhs[1])));
    const int nn_A     = static_cast<int>(mxGetScalar(prhs[2]));
    const int nt_A     = static_cast<int>(mxGetScalar(prhs[3]));
    const int32_t *iat_A  = mxGetInt32s(prhs[4]);
    const int32_t *ja_A   = mxGetInt32s(prhs[5]);
    const double  *coef_A = mxGetDoubles(prhs[6]);
    const int32_t *coef_S = mxGetInt32s(prhs[7]);
    const int32_t *fcn    = mxGetInt32s(prhs[8]);
    const int nn_I     = static_cast<int>(mxGetScalar(prhs[9]));
    const int nc_I     = static_cast<int>(mxGetScalar(prhs[10]));

    int nt_I = 0;
    int *iat_I = nullptr, *ja_I = nullptr;
    double *coef_I = nullptr;

    int ierr = EXTI_prolongation(level, np, nn_A, nt_A, iat_A, ja_A, coef_A, coef_S,
                                 fcn, nn_I, nc_I, nt_I, iat_I, ja_I, coef_I);
    if (ierr != 0 || !iat_I) {
        delete[] iat_I; delete[] ja_I; delete[] coef_I;
        mexErrMsgIdAndTxt("EXTI_Prol:computeError", "EXTI_prolongation failed with error %d.", ierr);
    }

    // Skip explicit zeros to match MATLAB sparse() and avoid downstream filter corruption
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
                for (int k = iat_I[r]; k < iat_I[r + 1]; ++k) {
                    if (coef_I[k] != 0.0) tc[ja_I[k]]++;
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

        int nt_nz = 0;
        for (int c = 0; c < nc_I; ++c) nt_nz += col_totals[c];

        plhs[0] = mxCreateSparse(static_cast<mwSize>(nn_I), static_cast<mwSize>(nc_I),
                                 static_cast<mwSize>(std::max(1, nt_nz)), mxREAL);
        mwIndex *jc = mxGetJc(plhs[0]);
        mwIndex *ir = mxGetIr(plhs[0]);
        double  *pr = mxGetDoubles(plhs[0]);

        jc[0] = 0;
        for (int c = 0; c < nc_I; ++c) jc[c + 1] = jc[c] + static_cast<mwIndex>(col_totals[c]);

        if (nt_nz > 0) {
            #pragma omp parallel num_threads(nth_used)
            {
                int tid = omp_get_thread_num();
                int r0 = static_cast<int>((static_cast<long long>(tid) * nn_I) / nth_used);
                int r1 = static_cast<int>((static_cast<long long>(tid + 1) * nn_I) / nth_used);
                int *to = &thread_offsets[static_cast<std::size_t>(tid) * nc_I];
                for (int c = 0; c < nc_I; ++c) to[c] += static_cast<int>(jc[c]);
                for (int r = r0; r < r1; ++r) {
                    for (int k = iat_I[r]; k < iat_I[r + 1]; ++k) {
                        double val = coef_I[k];
                        if (val == 0.0) continue;
                        int dst = to[ja_I[k]]++;
                        ir[dst] = static_cast<mwIndex>(r);
                        pr[dst] = val;
                    }
                }
            }
        }
    } else {
        plhs[0] = mxCreateSparse(static_cast<mwSize>(nn_I), static_cast<mwSize>(nc_I), 1, mxREAL);
        mwIndex *jc = mxGetJc(plhs[0]);
        for (int c = 0; c <= nc_I; ++c) jc[c] = 0;
    }

    delete[] iat_I;
    delete[] ja_I;
    delete[] coef_I;
}
