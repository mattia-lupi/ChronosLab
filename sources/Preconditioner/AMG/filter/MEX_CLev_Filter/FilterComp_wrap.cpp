//----------------------------------------------------------------------------------------
// FilterComp_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [nt_AC, iat_AC, ja_AC, coef_AC] =
//       FilterComp_wrap(np, tau, nn_A, iat_A, ja_A, coef_A,
//                       nt_patt, iat_patt, ja_patt, ntv, TV)
//
//----------------------------------------------------------------------------------------

#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "FilterComp.h"

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
// RAII guard for malloc-allocated pointers
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
static int FilterComp_wrapper(const int nthreads, const double tau,
                              const int nn_A,
                              int *iat_A, int *ja_A, double *coef_A,
                              const int nt_patt,
                              const int *iat_patt, const int *ja_patt,
                              const int ntv, const double *const *TV,
                              int &nt_AC, int *&iat_AC,
                              int *&ja_AC, double *&coef_AC)
{
    return FilterComp(nthreads, tau, nn_A, iat_A, ja_A, coef_A,
                      nt_patt, iat_patt, ja_patt,
                      ntv, TV,
                      nt_AC, iat_AC, ja_AC, coef_AC);
}

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Read input
        // -----------------------------------------------------------------------
        const int    np      = static_cast<int>   (TypedArray<double>(inputs[ 0])[0]);
        const double tau     = static_cast<double>(TypedArray<double>(inputs[ 1])[0]);
        const int    nn_A    = static_cast<int>   (TypedArray<double>(inputs[ 2])[0]);
        const int    nt_patt = static_cast<int>   (TypedArray<double>(inputs[ 6])[0]);
        const int    ntv     = static_cast<int>   (TypedArray<double>(inputs[ 9])[0]);

        const TypedArray<int32_t> iat_A_arr    = inputs[3];
        const TypedArray<int32_t> ja_A_arr     = inputs[4];
        const TypedArray<double>  coef_A_arr   = inputs[5];
        const TypedArray<int32_t> iat_patt_arr = inputs[7];
        const TypedArray<int32_t> ja_patt_arr  = inputs[8];
        const TypedArray<double>  TV_arr       = inputs[10];

        std::vector<int32_t> iat_A_vec   (iat_A_arr.begin(),    iat_A_arr.end());
        std::vector<int32_t> ja_A_vec    (ja_A_arr.begin(),     ja_A_arr.end());
        std::vector<double>  coef_A_vec  (coef_A_arr.begin(),   coef_A_arr.end());
        std::vector<int32_t> iat_patt_vec(iat_patt_arr.begin(), iat_patt_arr.end());
        std::vector<int32_t> ja_patt_vec (ja_patt_arr.begin(),  ja_patt_arr.end());

        std::vector<double> tv_buffer_vec(TV_arr.begin(), TV_arr.end());

        double **TV_2D_raw = static_cast<double**>(
                                 malloc(static_cast<std::size_t>(nn_A) *
                                        sizeof(double*)));
        if (!TV_2D_raw)
            throwError("FilterComp:allocError",
                       "Failed to allocate TV_2D pointer array — out of memory.");

        MallocGuard g_tv2d(TV_2D_raw);

        for (int i = 0; i < nn_A; ++i)
            TV_2D_raw[i] = tv_buffer_vec.data() +
                           static_cast<std::size_t>(i) *
                           static_cast<std::size_t>(ntv);

        // -----------------------------------------------------------------------
        // Call the computational wrapper
        // -----------------------------------------------------------------------
        int     nt_AC       = 0;
        int    *iat_AC_raw  = nullptr;
        int    *ja_AC_raw   = nullptr;
        double *coef_AC_raw = nullptr;

        int ierr = FilterComp_wrapper(np, tau, nn_A,
                                      iat_A_vec.data(),    ja_A_vec.data(),
                                      coef_A_vec.data(),
                                      nt_patt,
                                      iat_patt_vec.data(), ja_patt_vec.data(),
                                      ntv, TV_2D_raw,
                                      nt_AC, iat_AC_raw, ja_AC_raw, coef_AC_raw);

        // g_tv2d destructs here — TV_2D_raw freed; tv_buffer_vec auto-freed by vector

        MallocGuard g_iat (iat_AC_raw);
        MallocGuard g_ja  (ja_AC_raw);
        MallocGuard g_coef(coef_AC_raw);

        if (!iat_AC_raw || !ja_AC_raw || !coef_AC_raw)
            throwError("FilterComp:nullPointer",
                       "Kernel returned a null pointer — likely an allocation failure.");

        if (ierr != 0)
            throwError("FilterComp:computeError",
                       "FilterComp returned error code: " + std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB output arrays.
        // -----------------------------------------------------------------------
        const std::size_t sz_iat = static_cast<std::size_t>(nn_A + 1);
        const std::size_t sz_nt  = static_cast<std::size_t>(nt_AC);

        // --- nt_AC : scalar double --------------------------------------------
        TypedArray<double> out_nt_AC =
            factory.createScalar<double>(static_cast<double>(nt_AC));

        // --- iat_AC : double, length nn_A+1, 0-based → 1-based ---------------
        TypedArray<double> out_iat_AC = factory.createArray<double>({sz_iat, 1});
        {
            auto it = out_iat_AC.begin();
            for (int i = 0; i <= nn_A; ++i, ++it)
                *it = static_cast<double>(iat_AC_raw[i] + 1);
        }

        // --- ja_AC : double, length nt_AC, 0-based → 1-based -----------------
        TypedArray<double> out_ja_AC = factory.createArray<double>({sz_nt, 1});
        {
            auto it = out_ja_AC.begin();
            for (std::size_t i = 0; i < sz_nt; ++i, ++it)
                *it = static_cast<double>(ja_AC_raw[i] + 1);
        }

        // --- coef_AC : double, length nt_AC -----------------------------------
        TypedArray<double> out_coef_AC = factory.createArray<double>({sz_nt, 1});
        std::copy(coef_AC_raw, coef_AC_raw + sz_nt, out_coef_AC.begin());

        // MallocGuards destruct here — iat/ja/coef_AC_raw freed automatically

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(out_nt_AC);
        outputs[1] = std::move(out_iat_AC);
        outputs[2] = std::move(out_ja_AC);
        outputs[3] = std::move(out_coef_AC);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 11)
            throwError("FilterComp:badInputCount",
                       "Expected 11 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("FilterComp:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Scalar double inputs: np(0), tau(1), nn_A(2), nt_patt(6), ntv(9)
        for (std::size_t i : {0u, 1u, 2u, 6u, 9u})
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("FilterComp:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // int32 array inputs: iat_A(3), ja_A(4), iat_patt(7), ja_patt(8)
        for (std::size_t i : {3u, 4u, 7u, 8u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("FilterComp:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // double array inputs: coef_A(5), TV(10)
        for (std::size_t i : {5u, 10u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("FilterComp:badArray",
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
