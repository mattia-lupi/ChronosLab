//----------------------------------------------------------------------------------------
// NSY_rFSAI_compute.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [iat_FL, ja_FL, coef_FL, iat_FU, ja_FU, coef_FU] = ...
//       NSY_rFSAI_compute(nstep, step_size, epsilon, nn_A, iat_A, ja_A, coef_A)
//
// Build command (compile.m):
//   See compile.m — ensure -R2018a is on its OWN line (see compile.m fix notes)
//
// FIXES APPLIED (vs previous version — all found by Clang on R2025b/macOS)
// -----------------------------------------------------------------------
// [FIX-A] mexPrintf() is NOT declared when using the pure C++ MEX API
//         (mex.hpp does not pull in mex.h).  Replaced with a helper
//         mprint() that routes through getEngine()->feval(u"fprintf").
//
// [FIX-B] TypedArray<T>::operator[] returns a PROXY object
//         (ArrayElementTypedRef<T>), NOT a real reference.  Taking its
//         address (&arr[0]) is ill-formed.  Fixed by copying input arrays
//         into std::vector<T> and passing vec.data() to the kernel.
//         This is the only correct way to get a raw pointer from a
//         TypedArray in the C++ MEX API.
//
// [FIX-C] ArgumentList (MexIORange) has size() and operator[] defined as
//         NON-const member functions, so validateArguments() cannot take
//         its arguments as  const ArgumentList&.  Changed to non-const
//         references  (ArgumentList&).
//
// [FIX-D] compile.m had '-R2018a' placed AFTER '...' on the same line:
//              '-lmwlapack', ... '-R2018a',...
//         In MATLAB, everything after '...' on the same physical line is
//         silently ignored — so -R2018a was never passed to mex, meaning
//         the interleaved complex ABI was never enabled.  Fixed in
//         compile.m (see companion file).
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#if defined PRINT
    static constexpr bool dump = true;
#else
    static constexpr bool dump = false;
#endif

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "Compute_nsy_rfsai.h"

#include <cstdlib>   // free()
#include <string>
#include <vector>
#include <algorithm> // std::copy

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// RAII guard for malloc-allocated pointers returned by the C kernel
//----------------------------------------------------------------------------------------
struct MallocGuard {
    void *ptr = nullptr;
    explicit MallocGuard(void *p) : ptr(p) {}
    ~MallocGuard() { if (ptr) free(ptr); }
    MallocGuard(const MallocGuard&)            = delete;
    MallocGuard& operator=(const MallocGuard&) = delete;
};

//----------------------------------------------------------------------------------------
// MexFunction — mandatory class name, discovered by mexAdapter.hpp
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    //------------------------------------------------------------------------------------
    // [FIX-A] mexPrintf is not declared in the pure C++ MEX API.
    //         Route diagnostic output through MATLAB's fprintf instead.
    // [FIX-E] createScalar<T> only accepts arithmetic types.
    //         Use createCharArray() for std::string arguments.
    //------------------------------------------------------------------------------------
    void mprint(const std::string& msg)
    {
        getEngine()->feval(u"fprintf", 0,
            std::vector<Array>{ factory.createCharArray(msg) });
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C] validateArguments takes non-const refs — ArgumentList
        //         methods (size, operator[]) are not const-qualified
        validateArguments(outputs, inputs);

        if (dump) mprint("*** NSY_rFSAI_compute (C++ MEX API) ***\n");

        // -----------------------------------------------------------------------
        // Read input scalars
        // -----------------------------------------------------------------------
        if (dump) mprint("- Get input scalars\n");

        const TypedArray<double> s0 = inputs[0];
        const TypedArray<double> s1 = inputs[1];
        const TypedArray<double> s2 = inputs[2];
        const TypedArray<double> s3 = inputs[3];

        const int    nstep     = static_cast<int>(s0[0]);
        const int    step_size = static_cast<int>(s1[0]);
        const double epsilon   = static_cast<double>(s2[0]);
        const int    nn_A      = static_cast<int>(s3[0]);

        // -----------------------------------------------------------------------
        // [FIX-B] TypedArray<T>::operator[] returns a PROXY, not a real ref.
        //         You cannot take its address.  The correct pattern is to copy
        //         into a std::vector<T> and hand vec.data() to the C kernel.
        //         Cost: one allocation + memcpy per input array — unavoidable
        //         when bridging the MATLAB Data API to a raw-pointer C kernel.
        // -----------------------------------------------------------------------
        if (dump) mprint("- Get input arrays\n");

        const TypedArray<int32_t> iat_A_arr  = inputs[4];
        const TypedArray<int32_t> ja_A_arr   = inputs[5];
        const TypedArray<double>  coef_A_arr = inputs[6];

        std::vector<int32_t> iat_A_vec (iat_A_arr.begin(),  iat_A_arr.end());
        std::vector<int32_t> ja_A_vec  (ja_A_arr.begin(),   ja_A_arr.end());
        std::vector<double>  coef_A_vec(coef_A_arr.begin(), coef_A_arr.end());

        int32_t *iat_A  = iat_A_vec.data();
        int32_t *ja_A   = ja_A_vec.data();
        double  *coef_A = coef_A_vec.data();

        // -----------------------------------------------------------------------
        // Call the C computational kernel
        // -----------------------------------------------------------------------
        if (dump) mprint("- Compute FL and FU entries\n");

        int32_t *iat_FL_raw  = nullptr, *ja_FL_raw   = nullptr;
        int32_t *iat_FU_raw  = nullptr, *ja_FU_raw   = nullptr;
        double  *coef_FL_raw = nullptr, *coef_FU_raw = nullptr;

        int ierr = Compute_nsy_rfsai(
                       nstep, step_size, epsilon, nn_A,
                       iat_A, ja_A, coef_A,
                       iat_FL_raw, ja_FL_raw,  coef_FL_raw,
                       iat_FU_raw, ja_FU_raw,  coef_FU_raw);

        // Guard every kernel-allocated pointer immediately
        MallocGuard g_iat_FL (iat_FL_raw);
        MallocGuard g_ja_FL  (ja_FL_raw);
        MallocGuard g_cFL    (coef_FL_raw);
        MallocGuard g_iat_FU (iat_FU_raw);
        MallocGuard g_ja_FU  (ja_FU_raw);
        MallocGuard g_cFU    (coef_FU_raw);

        // Null-pointer check (catches kernel malloc failures before deref)
        if (!iat_FL_raw || !ja_FL_raw  || !coef_FL_raw ||
            !iat_FU_raw || !ja_FU_raw  || !coef_FU_raw)
            throwError("NSY_rFSAI:nullPointer",
                       "Kernel returned a null pointer — likely an allocation failure.");

        if (ierr != 0)
            throwError("NSY_rFSAI:computeError",
                       "Compute_nsy_rfsai returned error code: " +
                       std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB TypedArray output objects
        // -----------------------------------------------------------------------
        if (dump) mprint("- Store FL and FU into output arrays\n");

        const std::size_t n1    = static_cast<std::size_t>(nn_A + 1);
        const std::size_t nt_FL = static_cast<std::size_t>(iat_FL_raw[nn_A]);
        const std::size_t nt_FU = static_cast<std::size_t>(iat_FU_raw[nn_A]);

        // --- iat_FL : int32, length nn_A+1, 0-based → 1-based ----------------
        TypedArray<int32_t> iat_FL_out = factory.createArray<int32_t>({1, n1});
        {
            auto it = iat_FL_out.begin();
            for (int k = 0; k <= nn_A; ++k, ++it)
                *it = iat_FL_raw[k] + 1;
        }

        // --- ja_FL : int32, length nt_FL, 0-based → 1-based ------------------
        TypedArray<int32_t> ja_FL_out = factory.createArray<int32_t>({1, nt_FL});
        {
            auto it = ja_FL_out.begin();
            for (std::size_t k = 0; k < nt_FL; ++k, ++it)
                *it = ja_FL_raw[k] + 1;
        }

        // --- coef_FL : double, length nt_FL -----------------------------------
        TypedArray<double> coef_FL_out = factory.createArray<double>({1, nt_FL});
        std::copy(coef_FL_raw, coef_FL_raw + nt_FL, coef_FL_out.begin());

        // --- iat_FU : int32, length nn_A+1, 0-based → 1-based ----------------
        TypedArray<int32_t> iat_FU_out = factory.createArray<int32_t>({1, n1});
        {
            auto it = iat_FU_out.begin();
            for (int k = 0; k <= nn_A; ++k, ++it)
                *it = iat_FU_raw[k] + 1;
        }

        // --- ja_FU : int32, length nt_FU, 0-based → 1-based ------------------
        TypedArray<int32_t> ja_FU_out = factory.createArray<int32_t>({1, nt_FU});
        {
            auto it = ja_FU_out.begin();
            for (std::size_t k = 0; k < nt_FU; ++k, ++it)
                *it = ja_FU_raw[k] + 1;
        }

        // --- coef_FU : double, length nt_FU -----------------------------------
        TypedArray<double> coef_FU_out = factory.createArray<double>({1, nt_FU});
        std::copy(coef_FU_raw, coef_FU_raw + nt_FU, coef_FU_out.begin());

        // MallocGuards destruct here — all six free() calls fire automatically

        // -----------------------------------------------------------------------
        // Hand ownership of output arrays to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(iat_FL_out);
        outputs[1] = std::move(ja_FL_out);
        outputs[2] = std::move(coef_FL_out);
        outputs[3] = std::move(iat_FU_out);
        outputs[4] = std::move(ja_FU_out);
        outputs[5] = std::move(coef_FU_out);

        if (dump) mprint("\n");
    }

private:

    //------------------------------------------------------------------------------------
    // [FIX-C] ArgumentList::size() and operator[] are NOT const-qualified in
    //         MexIORange — parameters must be non-const references
    //------------------------------------------------------------------------------------
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 7)
            throwError("NSY_rFSAI:badInputCount",
                       "Expected 7 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 6)
            throwError("NSY_rFSAI:badOutputCount",
                       "Expected 6 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Inputs 0–3: real double scalars
        for (std::size_t i = 0; i < 4; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("NSY_rFSAI:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // Inputs 4–5: int32 arrays
        for (std::size_t i = 4; i <= 5; ++i)
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("NSY_rFSAI:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // Input 6: double array
        if (inputs[6].getType() != ArrayType::DOUBLE)
            throwError("NSY_rFSAI:badArray",
                       "Input argument 7 must be a double array.");
    }

    // [FIX-E] createScalar<T> is restricted to arithmetic types.
    //         createCharArray() is the correct factory method for strings.
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
