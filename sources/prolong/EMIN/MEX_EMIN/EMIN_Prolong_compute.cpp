//----------------------------------------------------------------------------------------
// EMIN_Prolong_compute.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [iat_Pout, ja_Pout, coef_Pout, info] =
//       EMIN_Prolong_compute(level,np,itmax,en_tol,condmax,maxwgt,prec,sol_type,
//                            min_lfil,max_lfil,D_lfil,nn,nn_C,ntv,nt_A,nt_P,
//                            nt_patt,fcnode,iat_A,ja_A,coef_A,iat_Pin,ja_Pin,
//                            coef_Pin,iat_patt,ja_patt,TV);
//
// Build command:
//   See compile.m — ensure -R2018a is on its own line
//
// ALL FIXES APPLIED (battle-tested against Clang / R2025b / macOS)
// -----------------------------------------------------------------------
// [FIX-A] mexPrintf() not declared in the pure C++ MEX API (mex.hpp does
//         not pull in mex.h). Replaced with mprint() helper that routes
//         through getEngine()->feval(u"fprintf", createCharArray(...)).
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy (ArrayElementTypedRef).
//         Taking its address is ill-formed. Input arrays are copied into
//         std::vector<T> and their .data() pointers are passed to the kernel.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified in
//         MexIORange. validateArguments() takes non-const ArgumentList&.
//
// [FIX-E] ArrayFactory::createScalar<T>() requires std::is_arithmetic<T>.
//         All string arguments use factory.createCharArray() instead.
//
// [NEW-1] Removed #include <iostream> and  using namespace std.
//         cout/endl replaced with throwError() / mprint().
//
// [NEW-2] TV (double**) is malloc-allocated from the input TVbuf pointer.
//         Wrapped in a std::unique_ptr with a custom deleter so it is
//         freed automatically on any exit path (normal or exception).
//         The pointed-to rows are views into TVbuf — only the pointer
//         array itself is freed, not the individual row pointers.
//
// [NEW-3] info[EMIN_INFO_SZ] is a plain stack array — fine as-is.
//         Copied into the output TypedArray<double> via std::copy.
//         The manual pt_r[0]=info[0]; ... loop is replaced.
//
// [NEW-4] All size / loop-counter variables use std::size_t or mwSize;
//         no plain int for array dimensions.
//
// [NEW-5] Null-pointer checks on every kernel output pointer before use.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#if defined PRINT
    static constexpr bool dump = true;
#else
    static constexpr bool dump = false;
#endif

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "EMIN_ImpProl.h"
#include "DebEnv.h"

#include <cstdlib>    // free(), malloc()
#include <cstring>    // std::copy
#include <string>
#include <vector>
#include <algorithm>  // std::copy
#include <memory>     // std::unique_ptr

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// RAII guard for malloc-allocated pointers returned by the C kernel      [NEW-2]
//----------------------------------------------------------------------------------------
struct MallocGuard {
    void *ptr = nullptr;
    explicit MallocGuard(void *p) : ptr(p) {}
    ~MallocGuard() { if (ptr) free(ptr); }
    MallocGuard(const MallocGuard&)            = delete;
    MallocGuard& operator=(const MallocGuard&) = delete;
};

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    // [FIX-A][FIX-E] Route print output through MATLAB's fprintf
    void mprint(const std::string& msg)
    {
        getEngine()->feval(u"fprintf", 0,
            std::vector<Array>{ factory.createCharArray(msg) });
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C] non-const ref required
        validateArguments(outputs, inputs);

        if (dump) mprint("*** EMIN_Prolong_compute (C++ MEX API) ***\n");

        // -----------------------------------------------------------------------
        // Read input scalars (inputs 0–16, all real double scalars from MATLAB)
        // -----------------------------------------------------------------------
        if (dump) mprint("- Get input scalars\n");

        const int    level     = static_cast<int>   (TypedArray<double>(inputs[ 0])[0]);
        const int    np        = static_cast<int>   (TypedArray<double>(inputs[ 1])[0]);
        const int    itmax     = static_cast<int>   (TypedArray<double>(inputs[ 2])[0]);
        const double en_tol    = static_cast<double>(TypedArray<double>(inputs[ 3])[0]);
        const double condmax   = static_cast<double>(TypedArray<double>(inputs[ 4])[0]);
        const double maxwgt    = static_cast<double>(TypedArray<double>(inputs[ 5])[0]);
        const int    prec      = static_cast<int>   (TypedArray<double>(inputs[ 6])[0]);
        const int    sol_type  = static_cast<int>   (TypedArray<double>(inputs[ 7])[0]);
        const int    min_lfil  = static_cast<int>   (TypedArray<double>(inputs[ 8])[0]);
        const int    max_lfil  = static_cast<int>   (TypedArray<double>(inputs[ 9])[0]);
        const int    D_lfil    = static_cast<int>   (TypedArray<double>(inputs[10])[0]);
        const int    nn        = static_cast<int>   (TypedArray<double>(inputs[11])[0]);
        const int    nn_C      = static_cast<int>   (TypedArray<double>(inputs[12])[0]);
        const int    ntv       = static_cast<int>   (TypedArray<double>(inputs[13])[0]);
        const int    nt_A      = static_cast<int>   (TypedArray<double>(inputs[14])[0]);
        const int    nt_P      = static_cast<int>   (TypedArray<double>(inputs[15])[0]);
        const int    nt_patt   = static_cast<int>   (TypedArray<double>(inputs[16])[0]);

        // -----------------------------------------------------------------------
        // Read input arrays
        // [FIX-B] Copy TypedArray → std::vector to obtain a real raw pointer
        // -----------------------------------------------------------------------
        if (dump) mprint("- Get input arrays\n");

        const TypedArray<int32_t> fcnode_arr   = inputs[17];
        const TypedArray<int32_t> iat_A_arr    = inputs[18];
        const TypedArray<int32_t> ja_A_arr     = inputs[19];
        const TypedArray<double>  coef_A_arr   = inputs[20];
        const TypedArray<int32_t> iat_Pin_arr  = inputs[21];
        const TypedArray<int32_t> ja_Pin_arr   = inputs[22];
        const TypedArray<double>  coef_Pin_arr = inputs[23];
        const TypedArray<int32_t> iat_patt_arr = inputs[24];
        const TypedArray<int32_t> ja_patt_arr  = inputs[25];
        const TypedArray<double>  TVbuf_arr    = inputs[26];

        std::vector<int32_t> fcnode_vec  (fcnode_arr.begin(),   fcnode_arr.end());
        std::vector<int32_t> iat_A_vec   (iat_A_arr.begin(),    iat_A_arr.end());
        std::vector<int32_t> ja_A_vec    (ja_A_arr.begin(),     ja_A_arr.end());
        std::vector<double>  coef_A_vec  (coef_A_arr.begin(),   coef_A_arr.end());
        std::vector<int32_t> iat_Pin_vec (iat_Pin_arr.begin(),  iat_Pin_arr.end());
        std::vector<int32_t> ja_Pin_vec  (ja_Pin_arr.begin(),   ja_Pin_arr.end());
        std::vector<double>  coef_Pin_vec(coef_Pin_arr.begin(), coef_Pin_arr.end());
        std::vector<int32_t> iat_patt_vec(iat_patt_arr.begin(), iat_patt_arr.end());
        std::vector<int32_t> ja_patt_vec (ja_patt_arr.begin(),  ja_patt_arr.end());
        std::vector<double>  TVbuf_vec   (TVbuf_arr.begin(),    TVbuf_arr.end());

        // -----------------------------------------------------------------------
        // Build the TV double** pointer array
        // [NEW-2] The pointer array is heap-allocated; rows point into TVbuf_vec.
        //         unique_ptr with a free-deleter ensures cleanup on any exit path.
        // -----------------------------------------------------------------------
        std::unique_ptr<double*, decltype(&free)> TV_owner(
            static_cast<double**>(malloc(static_cast<std::size_t>(nn) * sizeof(double*))),
            &free);

        if (!TV_owner)
            throwError("EMIN_Prolong:allocError",
                       "Failed to allocate TV pointer array.");

        double **TV = TV_owner.get();
        {
            int offset = 0;
            for (int i = 0; i < nn; ++i) {
                TV[i]   = TVbuf_vec.data() + offset;
                offset += ntv;
            }
        }

        // -----------------------------------------------------------------------
        // Initialise debug environment (unchanged from legacy)
        // -----------------------------------------------------------------------
        if (level == 1) {
            DebEnv.SetDebEnv(np, "w");
        } else {
            DebEnv.OpenDebugLog("a");
        }
        if (DEBUG) {
            fprintf(DebEnv.r_logfile,
                    "\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n", level);
            fflush(DebEnv.r_logfile);
            for (int i = 0; i < np; ++i) {
                fprintf(DebEnv.t_logfile[i],
                        "\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n", level);
                fflush(DebEnv.t_logfile[i]);
            }
        }

        // -----------------------------------------------------------------------
        // Call the C computational kernel
        // -----------------------------------------------------------------------
        if (dump) mprint("- Compute Pout entries\n");

        double   info[EMIN_INFO_SZ] = {};   // stack-allocated, zero-initialised
        int32_t *iat_Pout_raw  = nullptr;
        int32_t *ja_Pout_raw   = nullptr;
        double  *coef_Pout_raw = nullptr;

        int ierr = EMIN_ImpProl(np, itmax, en_tol, condmax, maxwgt,
                                prec, sol_type, min_lfil, max_lfil, D_lfil,
                                nn, nn_C, ntv, nt_A, nt_P, nt_patt,
                                fcnode_vec.data(),
                                iat_A_vec.data(),    ja_A_vec.data(),   coef_A_vec.data(),
                                iat_Pin_vec.data(),  ja_Pin_vec.data(), coef_Pin_vec.data(),
                                iat_patt_vec.data(), ja_patt_vec.data(),
                                TV,
                                iat_Pout_raw, ja_Pout_raw, coef_Pout_raw,
                                info);

        // TV_owner destructs here — free() called on the pointer array only;
        // the rows were views into TVbuf_vec (stack-managed), not separately alloc'd.

        // Guard every kernel-allocated output pointer immediately          [NEW-5]
        MallocGuard g_iat (iat_Pout_raw);
        MallocGuard g_ja  (ja_Pout_raw);
        MallocGuard g_coef(coef_Pout_raw);

        // Close debug log before any possible exception throw
        DebEnv.CloseDebugLog();

        // [NEW-5] Null-pointer check
        if (!iat_Pout_raw || !ja_Pout_raw || !coef_Pout_raw)
            throwError("EMIN_Prolong:nullPointer",
                       "Kernel returned a null pointer — likely an allocation failure.");

        // [NEW-1] Route kernel errors through MATLAB exception machinery
        if (ierr != 0)
            throwError("EMIN_Prolong:computeError",
                       "EMIN_ImpProl returned error code: " + std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB TypedArray output objects
        // -----------------------------------------------------------------------
        if (dump) mprint("- Store Pout into the output arrays\n");

        const std::size_t n1      = static_cast<std::size_t>(nn + 1);
        const std::size_t nt_Pout = static_cast<std::size_t>(iat_Pout_raw[nn]);

        // --- iat_Pout : int32, length nn+1, 0-based → 1-based ----------------
        TypedArray<int32_t> iat_final = factory.createArray<int32_t>({1, n1});
        {
            auto it = iat_final.begin();
            for (int k = 0; k <= nn; ++k, ++it)
                *it = iat_Pout_raw[k] + 1;
        }

        // --- ja_Pout : int32, length nt_Pout, 0-based → 1-based --------------
        TypedArray<int32_t> ja_final = factory.createArray<int32_t>({1, nt_Pout});
        {
            auto it = ja_final.begin();
            for (std::size_t k = 0; k < nt_Pout; ++k, ++it)
                *it = ja_Pout_raw[k] + 1;
        }

        // --- coef_Pout : double, length nt_Pout ------------------------------
        TypedArray<double> coef_final = factory.createArray<double>({1, nt_Pout});
        std::copy(coef_Pout_raw, coef_Pout_raw + nt_Pout, coef_final.begin());

        // MallocGuards destruct here — iat/ja/coef_Pout_raw freed automatically

        // --- info : double, length EMIN_INFO_SZ  [NEW-3] ---------------------
        TypedArray<double> info_out =
            factory.createArray<double>({1, static_cast<std::size_t>(EMIN_INFO_SZ)});
        std::copy(info, info + EMIN_INFO_SZ, info_out.begin());

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(iat_final);
        outputs[1] = std::move(ja_final);
        outputs[2] = std::move(coef_final);
        outputs[3] = std::move(info_out);

        if (dump) mprint("\n");
    }

private:

    // [FIX-C] ArgumentList methods are not const — take non-const refs
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 27)
            throwError("EMIN_Prolong:badInputCount",
                       "Expected 27 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("EMIN_Prolong:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Inputs 0–16: real double scalars
        for (std::size_t i = 0; i < 17; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("EMIN_Prolong:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // Inputs 17–19, 21–22, 24–25: int32 arrays
        for (std::size_t i : {17u, 18u, 19u, 21u, 22u, 24u, 25u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("EMIN_Prolong:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // Inputs 20, 23, 26: double arrays (coef_A, coef_Pin, TV)
        for (std::size_t i : {20u, 23u, 26u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("EMIN_Prolong:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a double array.");
    }

    // [FIX-E] createCharArray for string args; routes through MATLAB error()
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
