//----------------------------------------------------------------------------------------
// EMIN_Prolong_compute.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [iat_Pout, ja_Pout, coef_Pout, info] =
//       EMIN_Prolong_compute(level,np,itmax,en_tol,condmax,prec,sol_type,
//                            min_lfil,max_lfil,D_lfil,nn,nn_C,ntv,
//                            nt_patt,fcnode,iat_A,ja_A,coef_A,iat_Pin,ja_Pin,
//                            coef_Pin,iat_patt,ja_patt,TV);
//
//----------------------------------------------------------------------------------------

#if defined PRINT
    static constexpr bool dump = true;
#else
    static constexpr bool dump = false;
#endif

#include <cstdint>
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

struct MallocGuard {
    void *ptr = nullptr;
    explicit MallocGuard(void *p) : ptr(p) {}
    ~MallocGuard() { if (ptr) free(ptr); }
    void* release() { void* p = ptr; ptr = nullptr; return p; }
    MallocGuard(const MallocGuard&)            = delete;
    MallocGuard& operator=(const MallocGuard&) = delete;
};

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    void mprint(const std::string& msg)
    {
        getEngine()->feval(u"fprintf", 0,
            std::vector<Array>{ factory.createCharArray(msg) });
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
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
        const int    prec      = static_cast<int>   (TypedArray<double>(inputs[ 5])[0]);
        const int    sol_type  = static_cast<int>   (TypedArray<double>(inputs[ 6])[0]);
        const int    nn        = static_cast<int>   (TypedArray<double>(inputs[ 7])[0]);
        const int    nn_C      = static_cast<int>   (TypedArray<double>(inputs[ 8])[0]);
        const int    ntv       = static_cast<int>   (TypedArray<double>(inputs[ 9])[0]);
        const int    nt_patt   = static_cast<int>   (TypedArray<double>(inputs[10])[0]);
        bool         verb      = static_cast<bool>  (TypedArray<double>(inputs[21])[0]);

        // -----------------------------------------------------------------------
        // Read input arrays
        // -----------------------------------------------------------------------
        if (dump) mprint("- Get input arrays\n");

        const TypedArray<int32_t> fcnode_arr   = inputs[11];
        const TypedArray<int32_t> iat_A_arr    = inputs[12];
        const TypedArray<int32_t> ja_A_arr     = inputs[13];
        const TypedArray<double>  coef_A_arr   = inputs[14];
        const TypedArray<int32_t> iat_Pin_arr  = inputs[15];
        const TypedArray<int32_t> ja_Pin_arr   = inputs[16];
        const TypedArray<double>  coef_Pin_arr = inputs[17];
        const TypedArray<int32_t> iat_patt_arr = inputs[18];
        const TypedArray<int32_t> ja_patt_arr  = inputs[19];
        const TypedArray<double>  TVbuf_arr    = inputs[20];

        // Direct raw pointers to MATLAB data buffers — ZERO memory allocation and copying
        const int32_t* fcnode_ptr   = (fcnode_arr.getNumberOfElements() > 0)   ? (&(*fcnode_arr.begin()))   : nullptr;
        const int32_t* iat_A_ptr    = (iat_A_arr.getNumberOfElements() > 0)    ? (&(*iat_A_arr.begin()))    : nullptr;
        const int32_t* ja_A_ptr     = (ja_A_arr.getNumberOfElements() > 0)     ? (&(*ja_A_arr.begin()))     : nullptr;
        const double*  coef_A_ptr   = (coef_A_arr.getNumberOfElements() > 0)   ? (&(*coef_A_arr.begin()))   : nullptr;
        const int32_t* iat_Pin_ptr  = (iat_Pin_arr.getNumberOfElements() > 0)  ? (&(*iat_Pin_arr.begin()))  : nullptr;
        const int32_t* ja_Pin_ptr   = (ja_Pin_arr.getNumberOfElements() > 0)   ? (&(*ja_Pin_arr.begin()))   : nullptr;
        const double*  coef_Pin_ptr = (coef_Pin_arr.getNumberOfElements() > 0) ? (&(*coef_Pin_arr.begin())) : nullptr;
        const int32_t* iat_patt_ptr = (iat_patt_arr.getNumberOfElements() > 0) ? (&(*iat_patt_arr.begin())) : nullptr;
        const int32_t* ja_patt_ptr  = (ja_patt_arr.getNumberOfElements() > 0)  ? (&(*ja_patt_arr.begin()))  : nullptr;
        const double*  TVbuf_ptr    = (TVbuf_arr.getNumberOfElements() > 0)    ? (&(*TVbuf_arr.begin()))    : nullptr;

        // -----------------------------------------------------------------------
        // Build the TV const double** pointer array view into TVbuf
        // -----------------------------------------------------------------------
        std::unique_ptr<const double*[], decltype(&free)> TV_owner(
            static_cast<const double**>(malloc(static_cast<std::size_t>(nn) * sizeof(const double*))),
            &free);

        if (!TV_owner)
            throwError("EMIN_Prolong:allocError",
                       "Failed to allocate TV pointer array.");

        const double **TV = TV_owner.get();
        {
            int offset = 0;
            for (int i = 0; i < nn; ++i) {
                TV[i]   = TVbuf_ptr + offset;
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

        int ierr = EMIN_ImpProl(np, itmax, en_tol, condmax,
                                prec, sol_type, nn, nn_C, ntv,
                                nt_patt, fcnode_ptr,
                                iat_A_ptr,    ja_A_ptr,   coef_A_ptr,
                                iat_Pin_ptr,  ja_Pin_ptr, coef_Pin_ptr,
                                iat_patt_ptr, ja_patt_ptr,
                                TV,
                                iat_Pout_raw, ja_Pout_raw, coef_Pout_raw,
                                info,verb);

        // TV_owner destructs automatically

        MallocGuard g_iat (iat_Pout_raw);
        MallocGuard g_ja  (ja_Pout_raw);
        MallocGuard g_coef(coef_Pout_raw);

        // Close debug log before any possible exception throw
        DebEnv.CloseDebugLog();

        if (!iat_Pout_raw || !ja_Pout_raw || !coef_Pout_raw)
            throwError("EMIN_Prolong:nullPointer",
                       "Kernel returned a null pointer — likely an allocation failure.");

        if (ierr != 0)
            throwError("EMIN_Prolong:computeError",
                       "EMIN_ImpProl returned error code: " + std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB TypedArray output objects (zero-copy buffer transfer)
        // -----------------------------------------------------------------------
        if (dump) mprint("- Store Pout into the output arrays\n");

        const std::size_t n1      = static_cast<std::size_t>(nn + 1);
        const std::size_t nt_Pout = static_cast<std::size_t>(iat_Pout_raw[nn]);

        // Convert 0-based to 1-based indexing in-place
        #pragma omp parallel for num_threads(np)
        for (int k = 0; k <= nn; ++k)
            iat_Pout_raw[k] += 1;

        #pragma omp parallel for num_threads(np)
        for (std::size_t k = 0; k < nt_Pout; ++k)
            ja_Pout_raw[k] += 1;

        // Release pointers from guards and transfer ownership directly to MATLAB
        g_iat.release();
        g_ja.release();
        g_coef.release();

        buffer_ptr_t<int32_t> iat_buf(iat_Pout_raw, [](int32_t* p) { if (p) free(p); });
        TypedArray<int32_t> iat_final = factory.createArrayFromBuffer<int32_t>({1, n1}, std::move(iat_buf));

        buffer_ptr_t<int32_t> ja_buf(ja_Pout_raw, [](int32_t* p) { if (p) free(p); });
        TypedArray<int32_t> ja_final = factory.createArrayFromBuffer<int32_t>({1, nt_Pout}, std::move(ja_buf));

        buffer_ptr_t<double> coef_buf(coef_Pout_raw, [](double* p) { if (p) free(p); });
        TypedArray<double> coef_final = factory.createArrayFromBuffer<double>({1, nt_Pout}, std::move(coef_buf));

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

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 22)
            throwError("EMIN_Prolong:badInputCount",
                       "Expected 22 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("EMIN_Prolong:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Inputs 0–10: real double scalars
        for (std::size_t i = 0; i < 11; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("EMIN_Prolong:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // Inputs 22: real double scalar
        if (inputs[21].getType() != ArrayType::DOUBLE ||
                inputs[21].getNumberOfElements() != 1)
                throwError("EMIN_Prolong:badScalar",
                           "Input argument 22 must be a real double scalar.");

        // Inputs 11–13, 15–16, 18–19: int32 arrays
        for (std::size_t i : {11u, 12u, 13u, 15u, 16u, 18u, 19u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("EMIN_Prolong:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // Inputs 14, 17, 20: double arrays (coef_A, coef_Pin, TV)
        for (std::size_t i : {14u, 17u, 20u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("EMIN_Prolong:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a double array.");
    }

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

