//----------------------------------------------------------------------------------------
// EMIN_Prolong_compute.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [iat_Pout, ja_Pout, coef_Pout, info] =
//       EMIN_Prolong_compute(level,np,itmax,en_tol,condmax,prec,sol_type,
//                            min_lfil,max_lfil,D_lfil,nn,nn_C,ntv,
//                            nt_patt,SPDflag,fcnode,iat_A,ja_A,coef_A,
//                            iat_Pin,ja_Pin,coef_Pin,iat_patt,ja_patt,TV);
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
        // Read input scalars (inputs 0–11, all real double scalars from MATLAB)
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
        const bool   SPDflag   = static_cast<bool>  (TypedArray<double>(inputs[11])[0]);
        bool         verb      = static_cast<bool>  (TypedArray<double>(inputs[22])[0]);

        // -----------------------------------------------------------------------
        // Read input arrays
        // -----------------------------------------------------------------------
        if (dump) mprint("- Get input arrays\n");

        TypedArray<int32_t> fcnode_arr   = std::move(inputs[12]);
        TypedArray<int32_t> iat_A_arr    = std::move(inputs[13]);
        TypedArray<int32_t> ja_A_arr     = std::move(inputs[14]);
        TypedArray<double>  coef_A_arr   = std::move(inputs[15]);
        TypedArray<int32_t> iat_Pin_arr  = std::move(inputs[16]);
        TypedArray<int32_t> ja_Pin_arr   = std::move(inputs[17]);
        TypedArray<double>  coef_Pin_arr = std::move(inputs[18]);
        TypedArray<int32_t> iat_patt_arr = std::move(inputs[19]);
        TypedArray<int32_t> ja_patt_arr  = std::move(inputs[20]);
        TypedArray<double>  TVbuf_arr    = std::move(inputs[21]);

        auto fcnode_buf   = fcnode_arr.release();
        auto iat_A_buf    = iat_A_arr.release();
        auto ja_A_buf     = ja_A_arr.release();
        auto coef_A_buf   = coef_A_arr.release();
        auto iat_Pin_buf  = iat_Pin_arr.release();
        auto ja_Pin_buf   = ja_Pin_arr.release();
        auto coef_Pin_buf = coef_Pin_arr.release();
        auto iat_patt_buf = iat_patt_arr.release();
        auto ja_patt_buf  = ja_patt_arr.release();
        auto TVbuf_buf    = TVbuf_arr.release();

        // -----------------------------------------------------------------------
        // Build the TV double** pointer array
        // -----------------------------------------------------------------------
        std::unique_ptr<double*, decltype(&free)> TV_owner(
            static_cast<double**>(malloc(static_cast<std::size_t>(nn) * sizeof(double*))),
            &free);

        if (!TV_owner)
            throwError("EMIN_Prolong:allocError",
                       "Failed to allocate TV pointer array.");

        double **TV = TV_owner.get();
        double *TVbuf_ptr = TVbuf_buf.get();
        {
            std::ptrdiff_t offset = 0;
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
                                nt_patt, fcnode_buf.get(), SPDflag,
                                iat_A_buf.get(),    ja_A_buf.get(),   coef_A_buf.get(),
                                iat_Pin_buf.get(),  ja_Pin_buf.get(), coef_Pin_buf.get(),
                                iat_patt_buf.get(), ja_patt_buf.get(),
                                TV,
                                iat_Pout_raw, ja_Pout_raw, coef_Pout_raw,
                                info,verb);

        // TV_owner destructs here — free() called on the pointer array only;
        // the rows were views into TVbuf_buf, not separately alloc'd.

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
        // Pack results into MATLAB TypedArray output objects
        // -----------------------------------------------------------------------
        if (dump) mprint("- Store Pout into the output arrays\n");

        const std::size_t n1      = static_cast<std::size_t>(nn + 1);
        const std::size_t nt_Pout = static_cast<std::size_t>(iat_Pout_raw[nn]);

        auto iat_buf  = factory.createBuffer<int32_t>(n1);
        auto ja_buf   = factory.createBuffer<int32_t>(nt_Pout);
        auto coef_buf = factory.createBuffer<double>(nt_Pout);

        int32_t *iat_dst  = iat_buf.get();
        int32_t *ja_dst   = ja_buf.get();
        double  *coef_dst = coef_buf.get();

        for (int k = 0; k <= nn; ++k)
            iat_dst[k] = iat_Pout_raw[k] + 1;

        for (std::size_t k = 0; k < nt_Pout; ++k) {
            ja_dst[k]   = ja_Pout_raw[k] + 1;
            coef_dst[k] = coef_Pout_raw[k];
        }

        TypedArray<int32_t> iat_final  = factory.createArrayFromBuffer<int32_t>({1, n1}, std::move(iat_buf));
        TypedArray<int32_t> ja_final   = factory.createArrayFromBuffer<int32_t>({1, nt_Pout}, std::move(ja_buf));
        TypedArray<double>  coef_final = factory.createArrayFromBuffer<double>({1, nt_Pout}, std::move(coef_buf));

        // MallocGuards destruct here — iat/ja/coef_Pout_raw freed automatically

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
        if (inputs.size() != 23)
            throwError("EMIN_Prolong:badInputCount",
                       "Expected 23 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("EMIN_Prolong:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Inputs 0–11: real double scalars
        for (std::size_t i = 0; i < 12; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("EMIN_Prolong:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // Inputs 23: real double scalar
        if (inputs[22].getType() != ArrayType::DOUBLE ||
                inputs[22].getNumberOfElements() != 1)
                throwError("EMIN_Prolong:badScalar",
                           "Input argument 22 must be a real double scalar.");

        // Inputs 12–14, 16–17, 19–20: int32 arrays
        for (std::size_t i : {12u, 13u, 14u, 16u, 17u, 19u, 20u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("EMIN_Prolong:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // Inputs 15, 18, 21: double arrays (coef_A, coef_Pin, TV)
        for (std::size_t i : {15u, 18u, 21u})
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

