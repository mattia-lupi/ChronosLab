//----------------------------------------------------------------------------------------
// cpt_Prolongation_EXTI.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [nt_I, iat_I, ja_I, coef_I] =
//       cpt_Prolongation_EXTI(level, np, vecstart, nn_A, nt_A,
//                             iat_A, ja_A, coef_A, coef_S,
//                             iat_C, ja_C, coef_C, fcnodes,
//                             nn_I, nc_I)
//
// Build command:
//   See compile.m — ensure -R2018a is on its own line
//
// ALL FIXES APPLIED (consistent with previous modernized MEX files)
// -----------------------------------------------------------------------
// [FIX-A] mexPrintf() / mexErrMsgIdAndTxt() not declared in the pure C++
//         MEX API. Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
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
// [NEW-1] mxGetPr() was used for integer arrays (vecstart, iat_A, ja_A,
//         coef_S, iat_C, ja_C, fcnodes) — this is a latent bug: mxGetPr
//         returns double*, casting to int* is undefined behaviour.
//         Replaced with TypedArray<int32_t> copies via std::vector<int32_t>.
//
// [NEW-2] Kernel output pointers (iat_I, ja_I, coef_I) are malloc-allocated
//         inside the kernel. Wrapped in MallocGuard immediately after the
//         kernel returns — guarantees free() on any exit path including
//         exceptions, replacing the manual free() calls at the end.
//
// [NEW-3] Null-pointer check on all kernel output pointers before use.
//
// [NEW-4] Output arrays iat_I, ja_I stored as double (same as original —
//         the "Typed Data Access NOT working on RUSSEL" workaround is
//         preserved). int→double cast made explicit via static_cast<double>.
//
// [NEW-5] std::size_t / mwSize used for all array sizes and loop bounds.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "EXTI_prolongation.h"

#include <cstdlib>    // free()
#include <string>
#include <vector>
#include <algorithm>  // std::copy

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// [NEW-2] RAII guard for kernel-allocated raw pointers
//----------------------------------------------------------------------------------------
struct MallocGuard {
    void *ptr = nullptr;
    explicit MallocGuard(void *p) : ptr(p) {}
    ~MallocGuard() { if (ptr) free(ptr); }
    MallocGuard(const MallocGuard&)            = delete;
    MallocGuard& operator=(const MallocGuard&) = delete;
};

//----------------------------------------------------------------------------------------
// Computational wrapper — logic unchanged from original
//----------------------------------------------------------------------------------------
static int cpt_Prolongation_EXTI(const int level, const int nthreads,
                                 const int *const vecstart,
                                 const int nn_A, const int nt_A,
                                 const int *const iat_A, const int *const ja_A,
                                 const double *const coef_A,
                                 const int *const coef_S,
                                 const int *const iat_C, const int *const ja_C,
                                 const double *const coef_C,
                                 const int *const fcnodes,
                                 const int nr_I, const int nc_I,
                                 int &nt_I, int *&iat_I,
                                 int *&ja_I, double *&coef_I)
{
    return EXTI_prolongation(level, nthreads, vecstart,
                             nn_A, nt_A, iat_A, ja_A, coef_A, coef_S,
                             iat_C, ja_C, coef_C, fcnodes,
                             nr_I, nc_I, nt_I, iat_I, ja_I, coef_I);
}

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Read input scalars
        // -----------------------------------------------------------------------
        const int level  = static_cast<int>(TypedArray<double>(inputs[ 0])[0]);
        const int np     = static_cast<int>(TypedArray<double>(inputs[ 1])[0]);
        const int nn_A   = static_cast<int>(TypedArray<double>(inputs[ 3])[0]);
        const int nt_A   = static_cast<int>(TypedArray<double>(inputs[ 4])[0]);
        const int nn_I   = static_cast<int>(TypedArray<double>(inputs[13])[0]);
        const int nc_I   = static_cast<int>(TypedArray<double>(inputs[14])[0]);

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-1] Copy input arrays into std::vector for raw pointer access.
        //   Integer arrays use TypedArray<int32_t> — NOT mxGetPr (which is double*).
        //   coef_A and coef_C are genuine double arrays.
        // -----------------------------------------------------------------------
        const TypedArray<int32_t> vecstart_arr = inputs[ 2];
        const TypedArray<int32_t> iat_A_arr    = inputs[ 5];
        const TypedArray<int32_t> ja_A_arr     = inputs[ 6];
        const TypedArray<double>  coef_A_arr   = inputs[ 7];
        const TypedArray<int32_t> coef_S_arr   = inputs[ 8];   // [NEW-1] int, not double
        const TypedArray<int32_t> iat_C_arr    = inputs[ 9];
        const TypedArray<int32_t> ja_C_arr     = inputs[10];
        const TypedArray<double>  coef_C_arr   = inputs[11];
        const TypedArray<int32_t> fcnodes_arr  = inputs[12];

        std::vector<int32_t> vecstart_vec(vecstart_arr.begin(), vecstart_arr.end());
        std::vector<int32_t> iat_A_vec   (iat_A_arr.begin(),    iat_A_arr.end());
        std::vector<int32_t> ja_A_vec    (ja_A_arr.begin(),     ja_A_arr.end());
        std::vector<double>  coef_A_vec  (coef_A_arr.begin(),   coef_A_arr.end());
        std::vector<int32_t> coef_S_vec  (coef_S_arr.begin(),   coef_S_arr.end());
        std::vector<int32_t> iat_C_vec   (iat_C_arr.begin(),    iat_C_arr.end());
        std::vector<int32_t> ja_C_vec    (ja_C_arr.begin(),     ja_C_arr.end());
        std::vector<double>  coef_C_vec  (coef_C_arr.begin(),   coef_C_arr.end());
        std::vector<int32_t> fcnodes_vec (fcnodes_arr.begin(),  fcnodes_arr.end());

        // -----------------------------------------------------------------------
        // Call the computational wrapper
        // -----------------------------------------------------------------------
        int     nt_I     = 0;
        int    *iat_I_raw  = nullptr;
        int    *ja_I_raw   = nullptr;
        double *coef_I_raw = nullptr;

        int ierr = cpt_Prolongation_EXTI(
                       level, np,
                       vecstart_vec.data(),
                       nn_A, nt_A,
                       iat_A_vec.data(),  ja_A_vec.data(),  coef_A_vec.data(),
                       coef_S_vec.data(),
                       iat_C_vec.data(),  ja_C_vec.data(),  coef_C_vec.data(),
                       fcnodes_vec.data(),
                       nn_I, nc_I,
                       nt_I, iat_I_raw, ja_I_raw, coef_I_raw);

        // [NEW-2] Guard kernel-allocated output pointers immediately
        MallocGuard g_iat (iat_I_raw);
        MallocGuard g_ja  (ja_I_raw);
        MallocGuard g_coef(coef_I_raw);

        // [NEW-3] Null-pointer check before any dereference
        if (!iat_I_raw || !ja_I_raw || !coef_I_raw)
            throwError("EXTI_Prol:nullPointer",
                       "Kernel returned a null pointer — likely an allocation failure.");

        // [FIX-A] Route kernel error through MATLAB exception machinery
        if (ierr != 0)
            throwError("EXTI_Prol:computeError",
                       "cpt_Prolongation_EXTI returned error code: " +
                       std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB output arrays.
        // [NEW-4] iat_I and ja_I stored as double — preserved from original
        //         ("Typed Data Access NOT working on RUSSEL"); int→double explicit.
        // [NEW-5] std::size_t for all sizes and loop bounds
        // -----------------------------------------------------------------------
        const std::size_t sz_iat = static_cast<std::size_t>(nn_I + 1);
        const std::size_t sz_nt  = static_cast<std::size_t>(nt_I);

        // --- nt_I : scalar double ---------------------------------------------
        TypedArray<double> out_nt_I =
            factory.createScalar<double>(static_cast<double>(nt_I));

        // --- iat_I : double, length nn_I+1, 0-based → 1-based ----------------
        TypedArray<double> out_iat_I = factory.createArray<double>({sz_iat, 1});
        {
            auto it = out_iat_I.begin();
            for (int i = 0; i <= nn_I; ++i, ++it)
                *it = static_cast<double>(iat_I_raw[i] + 1);
        }

        // --- ja_I : double, length nt_I, 0-based → 1-based -------------------
        TypedArray<double> out_ja_I = factory.createArray<double>({sz_nt, 1});
        {
            auto it = out_ja_I.begin();
            for (std::size_t i = 0; i < sz_nt; ++i, ++it)
                *it = static_cast<double>(ja_I_raw[i] + 1);
        }

        // --- coef_I : double, length nt_I -------------------------------------
        TypedArray<double> out_coef_I = factory.createArray<double>({sz_nt, 1});
        std::copy(coef_I_raw, coef_I_raw + sz_nt, out_coef_I.begin());

        // MallocGuards destruct here — iat_I, ja_I, coef_I freed automatically [NEW-2]

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(out_nt_I);
        outputs[1] = std::move(out_iat_I);
        outputs[2] = std::move(out_ja_I);
        outputs[3] = std::move(out_coef_I);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 15)
            throwError("EXTI_Prol:badInputCount",
                       "Expected 15 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("EXTI_Prol:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Scalar double inputs: level(0), np(1), nn_A(3), nt_A(4), nn_I(13), nc_I(14)
        for (std::size_t i : {0u, 1u, 3u, 4u, 13u, 14u})
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("EXTI_Prol:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // int32 array inputs: vecstart(2), iat_A(5), ja_A(6), coef_S(8),
        //                     iat_C(9), ja_C(10), fcnodes(12)
        for (std::size_t i : {2u, 5u, 6u, 8u, 9u, 10u, 12u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("EXTI_Prol:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // double array inputs: coef_A(7), coef_C(11)
        for (std::size_t i : {7u, 11u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("EXTI_Prol:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a double array.");
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
