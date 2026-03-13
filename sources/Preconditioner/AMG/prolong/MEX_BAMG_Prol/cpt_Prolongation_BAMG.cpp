//----------------------------------------------------------------------------------------
// cpt_Prolongation_BAMG.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [nt_I, iat_I, ja_I, coef_I, c_mark] =
//       cpt_Prolongation_BAMG(level, np, itmax_vol, dist_min, dist_max, mmax,
//                             maxcond, maxrownrm, tol_vol, eps_prol,
//                             nn_S, nt_S, iat_S, ja_S, coef_S,
//                             ntv, fcnodes, TV, nn_I, nc_I)
//
// Build command:
//   See compile.m — ensure -R2018a is on its own line
//
// ALL FIXES APPLIED (consistent with NSY_rFSAI_compute and EMIN_Prolong_compute)
// -----------------------------------------------------------------------
// [FIX-A] mexPrintf() / mexErrMsgIdAndTxt() not declared in the pure C++
//         MEX API. Replaced with mprint() / throwError() helpers that route
//         through getEngine()->feval() with factory.createCharArray().
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         Input arrays copied into std::vector<T>; .data() passed to kernel.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() requires std::is_arithmetic<T>.
//         All string arguments use factory.createCharArray() instead.
//
// [NEW-1] Removed #include <iostream> / using namespace std / cout.
//
// [NEW-2] TV_2D uses a DEEP copy (each row individually malloc'd) — as in
//         the original. All rows plus the outer pointer array are managed
//         by TV2DGuard, a custom RAII class that frees every row then the
//         outer array, guaranteed even if the kernel throws.
//
// [NEW-3] mxGetPr() is deprecated and was incorrectly used for int arrays
//         (iat_S, ja_S, coef_S, fcnodes). Replaced with typed copies via
//         TypedArray<int32_t> (int arrays) and TypedArray<double> (TV).
//         coef_S is kept as int32 since the kernel signature is int*const.
//
// [NEW-4] Output arrays iat_I, ja_I, c_mark were stored as double in the
//         original (comment: "Typed Data Access NOT working on RUSSEL").
//         Preserved as double outputs for MATLAB compatibility — the int→
//         double cast is kept explicit via static_cast<double>.
//
// [NEW-5] ierr check routes through throwError() — no bare return that
//         would leave plhs[] unset.
//
// [NEW-6] std::size_t / mwSize used for all array sizes and loop bounds.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "BAMG.h"
#include "BAMG_params.h"
#include "DebEnv.h"

#include <cstdlib>    // malloc(), free()
#include <string>
#include <vector>
#include <algorithm>  // std::copy

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// [NEW-2] RAII owner for the deep-copy TV_2D (double**) structure.
//         Each row is independently malloc'd — all rows plus the outer
//         pointer array are freed on destruction, in any exit path.
//----------------------------------------------------------------------------------------
class TV2DGuard {
public:
    double **ptr  = nullptr;
    int      rows = 0;

    TV2DGuard() = default;

    // Allocate outer array and deep-copy rows from a flat column-major buffer.
    // Layout convention (same as original): flat[i*ntv + j] → ptr[i][j]
    bool allocate(int nn_S, int ntv, const double *flat)
    {
        rows = nn_S;
        ptr  = static_cast<double**>(malloc(static_cast<std::size_t>(nn_S) *
                                            sizeof(double*)));
        if (!ptr) return false;

        for (int i = 0; i < nn_S; ++i) ptr[i] = nullptr;   // safe for partial-alloc free

        for (int i = 0; i < nn_S; ++i) {
            ptr[i] = static_cast<double*>(malloc(static_cast<std::size_t>(ntv) *
                                                 sizeof(double)));
            if (!ptr[i]) return false;
            for (int j = 0; j < ntv; ++j)
                ptr[i][j] = flat[static_cast<std::size_t>(i) * ntv + j];
        }
        return true;
    }

    ~TV2DGuard()
    {
        if (ptr) {
            for (int i = 0; i < rows; ++i) free(ptr[i]);
            free(ptr);
        }
    }

    TV2DGuard(const TV2DGuard&)            = delete;
    TV2DGuard& operator=(const TV2DGuard&) = delete;
};

//----------------------------------------------------------------------------------------
// Computational wrapper — unchanged logic, lives in the MEX file as before
//----------------------------------------------------------------------------------------
int cpt_Prolongation_BAMG_wrapper(const int level, const BAMG_params params,
                                  const int nthreads,
                                  const int nn_S, const int nt_S,
                                  const int *const iat_S, const int *const ja_S,
                                  const int *const coef_S, const int ntv,
                                  const int *const fcnodes,
                                  const double *const *const TV,
                                  const int nc_I, int &nt_I,
                                  std::vector<int> &vec_iat_I,
                                  std::vector<int> &vec_ja_I,
                                  std::vector<double> &vec_coef_I,
                                  std::vector<int> &vec_c_mark)
{
    // Init debug
    if (level == 1) {
        DebEnv.SetDebEnv(nthreads, "w");
    } else {
        DebEnv.OpenDebugLog("a");
    }
    if (DEBUG && BAMG_DEBUG) {
        fprintf(DebEnv.r_logfile,
                "\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n", level);
        fflush(DebEnv.r_logfile);
        for (int i = 0; i < nthreads; ++i) {
            fprintf(DebEnv.t_logfile[i],
                    "\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n", level);
            fflush(DebEnv.t_logfile[i]);
        }
    }

    // nn_L = 0 (useless without MPI), nn_C = nn_S
    iReg nn_L = 0;
    iReg nn_C = nn_S;

    int ierr = BAMG(params, nthreads, nn_L, nn_C, nn_S,
                    iat_S, ja_S, ntv, fcnodes, TV,
                    nt_I, vec_iat_I, vec_ja_I, vec_coef_I, vec_c_mark);

    DebEnv.CloseDebugLog();
    return ierr;
}

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    // [FIX-A][FIX-E]
    void mprint(const std::string& msg)
    {
        getEngine()->feval(u"fprintf", 0,
            std::vector<Array>{ factory.createCharArray(msg) });
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Read input scalars (inputs 0–11 and 18–19)
        // -----------------------------------------------------------------------
        const int    level      = static_cast<int>   (TypedArray<double>(inputs[ 0])[0]);
        const int    np         = static_cast<int>   (TypedArray<double>(inputs[ 1])[0]);
        const int    itmax_vol  = static_cast<int>   (TypedArray<double>(inputs[ 2])[0]);
        const int    dist_min   = static_cast<int>   (TypedArray<double>(inputs[ 3])[0]);
        const int    dist_max   = static_cast<int>   (TypedArray<double>(inputs[ 4])[0]);
        const int    mmax       = static_cast<int>   (TypedArray<double>(inputs[ 5])[0]);
        const double maxcond    = static_cast<double>(TypedArray<double>(inputs[ 6])[0]);
        const double maxrownrm  = static_cast<double>(TypedArray<double>(inputs[ 7])[0]);
        const double tol_vol    = static_cast<double>(TypedArray<double>(inputs[ 8])[0]);
        const double eps_prol   = static_cast<double>(TypedArray<double>(inputs[ 9])[0]);
        const int    nn_S       = static_cast<int>   (TypedArray<double>(inputs[10])[0]);
        const int    nt_S       = static_cast<int>   (TypedArray<double>(inputs[11])[0]);
        // inputs 12–17: arrays (see below)
        const int    ntv        = static_cast<int>   (TypedArray<double>(inputs[15])[0]);
        const int    nn_I       = static_cast<int>   (TypedArray<double>(inputs[18])[0]);
        const int    nc_I       = static_cast<int>   (TypedArray<double>(inputs[19])[0]);

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-3] Copy input arrays into std::vector for raw pointer access.
        //   iat_S, ja_S, coef_S, fcnodes → int32 (kernel takes int*const)
        //   TV (flat buffer)             → double
        // -----------------------------------------------------------------------
        const TypedArray<int32_t> iat_S_arr   = inputs[12];
        const TypedArray<int32_t> ja_S_arr    = inputs[13];
        const TypedArray<int32_t> coef_S_arr  = inputs[14];   // [NEW-3] int, not double
        const TypedArray<int32_t> fcnodes_arr = inputs[16];
        const TypedArray<double>  TV_arr      = inputs[17];

        std::vector<int32_t> iat_S_vec  (iat_S_arr.begin(),   iat_S_arr.end());
        std::vector<int32_t> ja_S_vec   (ja_S_arr.begin(),    ja_S_arr.end());
        std::vector<int32_t> coef_S_vec (coef_S_arr.begin(),  coef_S_arr.end());
        std::vector<int32_t> fcnodes_vec(fcnodes_arr.begin(), fcnodes_arr.end());
        std::vector<double>  TV_flat    (TV_arr.begin(),       TV_arr.end());

        // -----------------------------------------------------------------------
        // [NEW-2] Build deep-copy TV_2D — each row independently malloc'd.
        //         TV2DGuard frees everything on any exit path.
        // -----------------------------------------------------------------------
        TV2DGuard tv2d;
        if (!tv2d.allocate(nn_S, ntv, TV_flat.data()))
            throwError("BAMG_Prol:allocError",
                       "Failed to allocate TV_2D — out of memory.");

        // -----------------------------------------------------------------------
        // Fill BAMG_params struct
        // -----------------------------------------------------------------------
        BAMG_params params;
        params.verbosity  = VERB_LEV;
        params.itmax_vol  = itmax_vol;
        params.dist_min   = dist_min;
        params.dist_max   = dist_max;
        params.mmax       = mmax;
        params.maxcond    = maxcond;
        params.maxrownrm  = maxrownrm;
        params.tol_vol    = tol_vol;
        params.eps        = eps_prol;

        // -----------------------------------------------------------------------
        // Call the computational wrapper
        // Outputs come back as std::vector — already RAII-managed, no free() needed.
        // -----------------------------------------------------------------------
        int              nt_I = 0;
        std::vector<int>    vec_iat_I, vec_ja_I, vec_c_mark;
        std::vector<double> vec_coef_I;

        int ierr = cpt_Prolongation_BAMG_wrapper(
                       level, params, np,
                       nn_S, nt_S,
                       iat_S_vec.data(), ja_S_vec.data(), coef_S_vec.data(),
                       ntv, fcnodes_vec.data(),
                       const_cast<const double* const*>(tv2d.ptr),
                       nc_I, nt_I,
                       vec_iat_I, vec_ja_I, vec_coef_I, vec_c_mark);

        // tv2d destructs here — all TV_2D rows and outer array freed

        // [NEW-5] Route error through MATLAB exception machinery
        if (ierr != 0)
            throwError("BAMG_Prol:computeError",
                       "cpt_Prolongation_BAMG returned error code: " +
                       std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB output arrays.
        // [NEW-4] iat_I, ja_I, c_mark are stored as double — preserved for
        //         MATLAB-side compatibility (original comment: "Typed Data
        //         Access NOT working on RUSSEL").  int→double cast is explicit.
        // [NEW-6] std::size_t / mwSize for all sizes
        // -----------------------------------------------------------------------
        const std::size_t sz_iat   = static_cast<std::size_t>(nn_I + 1);
        const std::size_t sz_nt    = static_cast<std::size_t>(nt_I);
        const std::size_t sz_cmark = static_cast<std::size_t>(nn_I);

        // --- nt_I : scalar double ---------------------------------------------
        TypedArray<double> out_nt_I = factory.createScalar<double>(
                                          static_cast<double>(nt_I));

        // --- iat_I : double, length nn_I+1, 0-based → 1-based ----------------
        TypedArray<double> out_iat_I = factory.createArray<double>({sz_iat, 1});
        {
            auto it = out_iat_I.begin();
            for (int i = 0; i <= nn_I; ++i, ++it)
                *it = static_cast<double>(vec_iat_I[static_cast<std::size_t>(i)] + 1);
        }

        // --- ja_I : double, length nt_I, 0-based → 1-based -------------------
        TypedArray<double> out_ja_I = factory.createArray<double>({sz_nt, 1});
        {
            auto it = out_ja_I.begin();
            for (std::size_t i = 0; i < sz_nt; ++i, ++it)
                *it = static_cast<double>(vec_ja_I[i] + 1);
        }

        // --- coef_I : double, length nt_I -------------------------------------
        TypedArray<double> out_coef_I = factory.createArray<double>({sz_nt, 1});
        std::copy(vec_coef_I.begin(), vec_coef_I.end(), out_coef_I.begin());

        // --- c_mark : double, length nn_I -------------------------------------
        TypedArray<double> out_c_mark = factory.createArray<double>({sz_cmark, 1});
        {
            auto it = out_c_mark.begin();
            for (std::size_t i = 0; i < sz_cmark; ++i, ++it)
                *it = static_cast<double>(vec_c_mark[i]);
        }

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(out_nt_I);
        outputs[1] = std::move(out_iat_I);
        outputs[2] = std::move(out_ja_I);
        outputs[3] = std::move(out_coef_I);
        outputs[4] = std::move(out_c_mark);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 20)
            throwError("BAMG_Prol:badInputCount",
                       "Expected 20 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 5)
            throwError("BAMG_Prol:badOutputCount",
                       "Expected 5 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Scalar double inputs: 0–11, 15, 18, 19
        for (std::size_t i : {0u,1u,2u,3u,4u,5u,6u,7u,8u,9u,10u,11u,15u,18u,19u})
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("BAMG_Prol:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // int32 array inputs: iat_S(12), ja_S(13), coef_S(14), fcnodes(16)
        for (std::size_t i : {12u, 13u, 14u, 16u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("BAMG_Prol:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // double array input: TV(17)
        if (inputs[17].getType() != ArrayType::DOUBLE)
            throwError("BAMG_Prol:badArray",
                       "Input argument 18 (TV) must be a double array.");
    }

    // [FIX-E] createCharArray for strings; routes through MATLAB error()
    void throwError(const std::string& id, const std::string& msg)
    {
        getEngine()->feval(u"error", 0,
            std::vector<Array>{
                factory.createCharArray(id),
                factory.createCharArray(msg)
            });
    }
};

//----------------------------------------------------------------------------------------
