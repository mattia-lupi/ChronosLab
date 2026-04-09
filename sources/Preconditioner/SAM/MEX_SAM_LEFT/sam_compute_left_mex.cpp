// sam_compute_left_mex.cpp  —  LEFT Sparse Approximate Map (MEX, C++17, OpenMP)
//
// Interface (unchanged — MEX_sam_compute_left.m works as-is):
//
//   [row_N, col_N, val_N] = sam_compute_left_mex(iatk, jak, coefk,
//                                                 iat0, ja0, coef0,
//                                                 preproc)
//
// All index arrays (iatk, jak, iat0, ja0, preproc.s_idx, preproc.r_idx)
// are 1-based MATLAB convention; conversion to 0-based is done once,
// up-front, into plain C++ memory before any parallel work begins.
//
// Performance improvements over the original version
// ---------------------------------------------------
// 1. ALL preproc CellArray / StructArray data is extracted into flat
//    C++ CSR-style vectors BEFORE the parallel region.  The MATLAB API
//    is NOT called inside any hot or parallel loop.
//
// 2. LAPACK dgels_ workspace is queried ONCE with worst-case dimensions
//    before the parallel region.  Each thread pre-allocates a single
//    work buffer of that fixed size and reuses it for every row it owns.
//
// 3. Ak and A0 stay in CSR — no CSC conversion.  Left-SAM row problems
//    access Ak by row, which is exactly what CSR provides O(1).
//
// 4. Element lookups inside the LS assembly use std::lower_bound on the
//    sorted column arrays instead of linear scans over pair-vectors.
//
// 5. Thread-local Atmp / f / work vectors are declared once per thread
//    inside the omp parallel region and only grow, never shrink — one
//    heap allocation per thread across the entire sequence of rows.
//
// 6. Output is indexed by pre-computed row offsets (== pp.s_ptr) so
//    each thread writes to a non-overlapping slice with no locks.
//
// 7. Row/col structural indices are filled in one sequential pass before
//    the parallel loop so threads write only the value array.

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
#include <omp.h>

using namespace matlab::data;
using namespace matlab::mex;

// =========================================================================
// LAPACK: dgels_ — overdetermined/underdetermined least squares via QR.
// ptrdiff_t is what MATLAB's mwlapack expects on both ILP64 and LP64.
// =========================================================================
extern "C" {
void dgels_(const char*      trans,
            const ptrdiff_t* m,
            const ptrdiff_t* n,
            const ptrdiff_t* nrhs,
            double* a,       const ptrdiff_t* lda,
            double* b,       const ptrdiff_t* ldb,
            double* work,    const ptrdiff_t* lwork,
            ptrdiff_t* info);
}

// =========================================================================
// Minimal CSR container, all 0-based.
// =========================================================================
struct CSR {
    int n;                   // number of rows
    std::vector<int>    ptr; // size n+1, row pointers
    std::vector<int>    col; // column indices, sorted within each row
    std::vector<double> val;
};

// Build from 1-based MATLAB arrays iat (size n+1), ja (nnz), coef (nnz).
static CSR buildCSR(const TypedArray<double>& iat,
                    const TypedArray<double>& ja,
                    const TypedArray<double>& coef)
{
    CSR csr;
    csr.n   = static_cast<int>(iat.getNumberOfElements()) - 1;
    int nnz = static_cast<int>(ja.getNumberOfElements());

    csr.ptr.resize(csr.n + 1);
    csr.col.resize(nnz);
    csr.val.resize(nnz);

    for (int i = 0; i <= csr.n; ++i)
        csr.ptr[i] = static_cast<int>(iat[i]) - 1;   // 1-based → 0-based

    for (int p = 0; p < nnz; ++p) {
        csr.col[p] = static_cast<int>(ja[p]) - 1;
        csr.val[p] = coef[p];
    }
    return csr;
}

// =========================================================================
// Flat (cache-friendly) representation of the SAM preproc struct.
//
//   s_data[ s_ptr[i] .. s_ptr[i+1]-1 ]  =  s_i  (0-based column indices)
//   r_data[ r_ptr[i] .. r_ptr[i+1]-1 ]  =  r_i  (0-based, sorted)
// =========================================================================
struct PreprocFlat {
    int n, total_nnz;
    std::vector<int> s_ptr, s_data;
    std::vector<int> r_ptr, r_data;
};

// -------------------------------------------------------------------------
// Read preproc struct produced by MEX_sam_preprocess_left.
//
// NEW (fast) path — fields s_ptr / s_data / r_ptr / r_data are flat double
// vectors written by the updated sam_preprocess_left_mex.cpp.
// Extraction cost: 4 linear casts (memcpy-equivalent), O(nnz) total,
// with zero MATLAB CellArray API calls.
//
// LEGACY path — fields s_idx / r_idx are CellArrays (old MATLAB or old MEX).
// Kept for backward compatibility; costs O(n) MATLAB API calls.
// -------------------------------------------------------------------------
static PreprocFlat extractPreproc(const StructArray& preproc)
{
    PreprocFlat pp;

    TypedArray<double> n_arr   = preproc[0]["n"];
    TypedArray<double> nnz_arr = preproc[0]["nnz_total"];
    pp.n         = static_cast<int>(n_arr[0]);
    pp.total_nnz = static_cast<int>(nnz_arr[0]);

    // Detect which format the struct carries by checking for "s_ptr" field
    bool has_flat = false;
    for (const auto& f : preproc.getFieldNames())
       if (std::string(f) == "s_ptr") { has_flat = true; break; }

    if (has_flat) {
        // ---- NEW path: flat 0-based double arrays ----------------------
        TypedArray<double> sp = preproc[0]["s_ptr"];
        TypedArray<double> sd = preproc[0]["s_data"];
        TypedArray<double> rp = preproc[0]["r_ptr"];
        TypedArray<double> rd = preproc[0]["r_data"];

        int np1  = static_cast<int>(sp.getNumberOfElements());  // n+1
        int sns  = static_cast<int>(sd.getNumberOfElements());
        int rns  = static_cast<int>(rd.getNumberOfElements());

        pp.s_ptr.resize(np1);
        pp.s_data.resize(sns);
        pp.r_ptr.resize(np1);
        pp.r_data.resize(rns);

        // Arrays are already 0-based; one cast loop per array
        for (int i = 0; i < np1; ++i) pp.s_ptr[i]  = static_cast<int>(sp[i]);
        for (int i = 0; i < sns; ++i) pp.s_data[i] = static_cast<int>(sd[i]);
        for (int i = 0; i < np1; ++i) pp.r_ptr[i]  = static_cast<int>(rp[i]);
        for (int i = 0; i < rns; ++i) pp.r_data[i] = static_cast<int>(rd[i]);

    } else {
        // ---- LEGACY path: CellArray s_idx / r_idx (1-based) -----------
        CellArray s_cell = preproc[0]["s_idx"];
        CellArray r_cell = preproc[0]["r_idx"];

        pp.s_ptr.resize(pp.n + 1, 0);
        pp.r_ptr.resize(pp.n + 1, 0);
        pp.s_data.reserve(pp.total_nnz);
        pp.r_data.reserve(pp.total_nnz * 3);

        for (int i = 0; i < pp.n; ++i) {
            TypedArray<double> si = s_cell[i];
            TypedArray<double> ri = r_cell[i];
            int nsi = static_cast<int>(si.getNumberOfElements());
            int nri = static_cast<int>(ri.getNumberOfElements());

            pp.s_ptr[i + 1] = pp.s_ptr[i] + nsi;
            pp.r_ptr[i + 1] = pp.r_ptr[i] + nri;

            for (int k = 0; k < nsi; ++k)
                pp.s_data.push_back(static_cast<int>(si[k]) - 1);  // 1→0
            for (int k = 0; k < nri; ++k)
                pp.r_data.push_back(static_cast<int>(ri[k]) - 1);  // 1→0
        }
    }
    return pp;
}

// =========================================================================
// Query optimal dgels_ LWORK with worst-case dimensions.
// Called once, single-threaded, before entering the parallel region.
// =========================================================================
static int queryLwork(int max_m, int max_n)
{
    if (max_m <= 0 || max_n <= 0) return 128;

    // Build a small diagonal matrix so the query is non-degenerate
    int mn = std::min(max_m, max_n);
    std::vector<double> Atmp(max_m * max_n, 0.0);
    for (int k = 0; k < mn; ++k) Atmp[k * max_m + k] = 1.0;

    char      trans  = 'N';
    ptrdiff_t pm = max_m, pn = max_n, pnrhs = 1;
    ptrdiff_t plda = pm, pldb = std::max(pm, pn);
    ptrdiff_t plw  = -1;
    ptrdiff_t info = 0;
    double wq = 0.0;
    std::vector<double> ftmp(std::max(max_m, max_n), 0.0);

    dgels_(&trans, &pm, &pn, &pnrhs,
           Atmp.data(), &plda,
           ftmp.data(), &pldb,
           &wq, &plw, &info);

    return std::max(128, static_cast<int>(wq));
}

// =========================================================================
// MexFunction
// =========================================================================
class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        ArrayFactory factory;

        if (inputs.size() != 7) {
            getEngine()->feval(u"error", 0,
                std::vector<Array>({factory.createScalar(
                    "sam_compute_left_mex: requires "
                    "iatk, jak, coefk, iat0, ja0, coef0, preproc.")}));
            return;
        }

        // ---- Phase 1: extract all data from MATLAB into C++ memory -----
        // This entire phase is single-threaded and touches the MATLAB API.
        // Nothing after this point calls the MATLAB API inside a hot loop.

        TypedArray<double> iatk  = std::move(inputs[0]);
        TypedArray<double> jak   = std::move(inputs[1]);
        TypedArray<double> coefk = std::move(inputs[2]);
        TypedArray<double> iat0  = std::move(inputs[3]);
        TypedArray<double> ja0   = std::move(inputs[4]);
        TypedArray<double> coef0 = std::move(inputs[5]);
        StructArray preproc      = std::move(inputs[6]);

        const CSR        Ak = buildCSR(iatk, jak, coefk);
        const CSR        A0 = buildCSR(iat0, ja0, coef0);
        const PreprocFlat pp = extractPreproc(preproc);

        const int n         = pp.n;
        const int total_nnz = pp.total_nnz;

        // ---- Phase 2: compute worst-case LS dimensions, query LWORK ----
        int max_si = 0, max_ri = 0;
        for (int i = 0; i < n; ++i) {
            max_si = std::max(max_si, pp.s_ptr[i+1] - pp.s_ptr[i]);
            max_ri = std::max(max_ri, pp.r_ptr[i+1] - pp.r_ptr[i]);
        }
        const int LWORK = queryLwork(max_ri, max_si);

        // ---- Phase 3: allocate output COO arrays (plain C++ vectors) ---
        std::vector<double> row_N(total_nnz), col_N(total_nnz), val_N(total_nnz);

        // Fill structural indices sequentially (cheap, no data dependency)
        for (int i = 0; i < n; ++i) {
            for (int p = pp.s_ptr[i]; p < pp.s_ptr[i+1]; ++p) {
                row_N[p] = static_cast<double>(i + 1);              // 1-based
                col_N[p] = static_cast<double>(pp.s_data[p] + 1);  // 1-based
            }
        }

        // ---- Phase 4: parallel row loop --------------------------------
        // Dynamic scheduling with a chunk of 16 rows balances the variable
        // subproblem sizes without excessive scheduling overhead.
        #pragma omp parallel
        {
            // Thread-local persistent buffers.
            // Declared here so they survive across omp for iterations
            // and are only reallocated when a larger size is needed.
            std::vector<double> Atmp, f, work(LWORK);

            #pragma omp for schedule(dynamic, 16) nowait
            for (int i = 0; i < n; ++i)
            {
                const int s0   = pp.s_ptr[i];
                const int s1   = pp.s_ptr[i + 1];
                const int nnzi = s1 - s0;
                if (nnzi == 0) continue;

                const int r0  = pp.r_ptr[i];
                const int r1  = pp.r_ptr[i + 1];
                const int nri = r1 - r0;

                const int* si = pp.s_data.data() + s0;  // 0-based col indices
                const int* ri = pp.r_data.data() + r0;  // 0-based, sorted

                // Zero out the output values for this row up front
                for (int p = s0; p < s1; ++p) val_N[p] = 0.0;
                if (nri == 0) continue;

                // ---- Assemble Atmp = Ak(s_i, r_i)^T  -------------------
                // Size: nri (rows) × nnzi (cols), column-major for LAPACK.
                // Column ji of Atmp corresponds to row si[ji] of Ak.
                // We need Atmp[ji * nri + ri_idx] = Ak[ si[ji], ri[ri_idx] ].
                //
                // Ak is CSR → row si[ji] is a sorted list of (col, val) pairs.
                // For each entry in that row, binary-search its column in ri.
                const int m_ls = nri, n_ls = nnzi;
                Atmp.assign(m_ls * n_ls, 0.0);

                for (int ji = 0; ji < nnzi; ++ji) {
                    const int row_k = si[ji];
                    const int ks    = Ak.ptr[row_k];
                    const int ke    = Ak.ptr[row_k + 1];
                    for (int p = ks; p < ke; ++p) {
                        const int c = Ak.col[p];
                        // Binary search: ri is sorted, typically short
                        const int* pos = std::lower_bound(ri, ri + nri, c);
                        if (pos != ri + nri && *pos == c) {
                            Atmp[ji * m_ls + static_cast<int>(pos - ri)] = Ak.val[p];
                        }
                    }
                }

                // ---- Assemble f = A0(i, r_i) ----------------------------
                // Row i of A0 in CSR; binary-search each needed column.
                const int ldb = std::max(m_ls, n_ls);
                f.assign(ldb, 0.0);
                {
                    const int as = A0.ptr[i];
                    const int ae = A0.ptr[i + 1];
                    for (int p = as; p < ae; ++p) {
                        const int c = A0.col[p];
                        const int* pos = std::lower_bound(ri, ri + nri, c);
                        if (pos != ri + nri && *pos == c) {
                            f[static_cast<int>(pos - ri)] = A0.val[p];
                        }
                    }
                }

                // ---- Solve min_z || Atmp z - f ||_2 via dgels_ ----------
                char      trans  = 'N';
                ptrdiff_t pm = m_ls, pn = n_ls, pnrhs = 1;
                ptrdiff_t plda = pm, pldb = ldb;
                ptrdiff_t plw  = static_cast<ptrdiff_t>(work.size());
                ptrdiff_t info = 0;

                dgels_(&trans, &pm, &pn, &pnrhs,
                       Atmp.data(), &plda,
                       f.data(),    &pldb,
                       work.data(), &plw, &info);

                // Solution is in f[0 .. nnzi-1]; write to val_N
                for (int ji = 0; ji < nnzi; ++ji) {
                    double v  = f[ji];
                    val_N[s0 + ji] = (info == 0 && std::isfinite(v)) ? v : 0.0;
                }
            }
        } // end omp parallel

        // ---- Phase 5: copy C++ COO into MATLAB output arrays -----------
        // One linear copy; MATLAB array creation cannot be parallelised.
        const size_t sz = static_cast<size_t>(total_nnz);
        TypedArray<double> out_row = factory.createArray<double>({sz, 1});
        TypedArray<double> out_col = factory.createArray<double>({sz, 1});
        TypedArray<double> out_val = factory.createArray<double>({sz, 1});

        for (int p = 0; p < total_nnz; ++p) {
            out_row[p] = row_N[p];
            out_col[p] = col_N[p];
            out_val[p] = val_N[p];
        }

        outputs[0] = std::move(out_row);
        outputs[1] = std::move(out_col);
        outputs[2] = std::move(out_val);
    }
};
