//----------------------------------------------------------------------------------------
// FilterProl_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [nt_PF, iat_PF, ja_PF, coef_PF] =
//       FilterProl_wrap(np, perc, tol, nn_P, iat_P, ja_P, coef_P,
//                       nr_TV, ntv, TV)
//
//----------------------------------------------------------------------------------------

#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "FilterProl.h"

#include <cstdlib>    // malloc(), free()
#include <cstring>    // memcpy()
#include <string>
#include <vector>
#include <algorithm>  // std::copy

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

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Read input
        // -----------------------------------------------------------------------
        const int    np    = static_cast<int>   (TypedArray<double>(inputs[0])[0]);
        const double perc  = static_cast<double>(TypedArray<double>(inputs[1])[0]);
        const double tol   = static_cast<double>(TypedArray<double>(inputs[2])[0]);
        const int    nn_P  = static_cast<int>   (TypedArray<double>(inputs[3])[0]);
        const int    nr_TV = static_cast<int>   (TypedArray<double>(inputs[7])[0]);
        const int    ntv   = static_cast<int>   (TypedArray<double>(inputs[8])[0]);

        const TypedArray<int32_t> iat_P_arr = inputs[4];
        const TypedArray<int32_t> ja_P_arr  = inputs[5];
        const TypedArray<double>  coef_P_arr= inputs[6];
        const TypedArray<double>  TV_arr    = inputs[9];

        std::vector<int32_t> iat_P_vec (iat_P_arr.begin(),  iat_P_arr.end());
        std::vector<int32_t> ja_P_vec  (ja_P_arr.begin(),   ja_P_arr.end());
        std::vector<double>  coef_P_vec(coef_P_arr.begin(), coef_P_arr.end());

        // -----------------------------------------------------------------------
        //   Original approach:
        //     buffer  = malloc(nr_TV * ntv * sizeof(double))   ← flat data copy
        //     TV_2D   = malloc(nr_TV * sizeof(double*))         ← pointer array
        //     TV_2D[i] = buffer + i*ntv
        //     memcpy(buffer, TV, ...)
        //
        //   Modern approach:
        //     tv_buffer_vec  — std::vector owns the flat copy (no malloc/free)
        //     TV_2D_raw      — still malloc'd (kernel takes double**);
        //                      freed via MallocGuard on any exit path.
        // -----------------------------------------------------------------------
        const std::size_t tv_total = static_cast<std::size_t>(nr_TV) *
                                     static_cast<std::size_t>(ntv);

        // Flat contiguous copy of the TV input — owned by vector, auto-freed
        std::vector<double> tv_buffer_vec(TV_arr.begin(), TV_arr.end());

        // Outer pointer array — malloc'd because the kernel takes double**
        double **TV_2D_raw = static_cast<double**>(
                                 malloc(static_cast<std::size_t>(nr_TV) *
                                        sizeof(double*)));
        if (!TV_2D_raw)
            throwError("FilterProl:allocError",
                       "Failed to allocate TV_2D pointer array — out of memory.");

        MallocGuard g_tv2d(TV_2D_raw);   // freed on scope exit

        for (int i = 0; i < nr_TV; ++i)
            TV_2D_raw[i] = tv_buffer_vec.data() +
                           static_cast<std::size_t>(i) * static_cast<std::size_t>(ntv);

        // -----------------------------------------------------------------------
        // Call the computational kernel
        // -----------------------------------------------------------------------
        int     nt_PF       = 0;
        int    *iat_PF_raw  = nullptr;
        int    *ja_PF_raw   = nullptr;
        double *coef_PF_raw = nullptr;

        int ierr = FilterProl(np, perc, tol,
                              nn_P,
                              iat_P_vec.data(), ja_P_vec.data(), coef_P_vec.data(),
                              ntv, TV_2D_raw,
                              nt_PF, iat_PF_raw, ja_PF_raw, coef_PF_raw);

        // g_tv2d destructs here — TV_2D_raw freed; tv_buffer_vec auto-freed by vector

        MallocGuard g_iat (iat_PF_raw);
        MallocGuard g_ja  (ja_PF_raw);
        MallocGuard g_coef(coef_PF_raw);

        if (!iat_PF_raw || !ja_PF_raw || !coef_PF_raw)
            throwError("FilterProl:nullPointer",
                       "Kernel returned a null pointer — likely an allocation failure.");

        if (ierr != 0)
            throwError("FilterProl:computeError",
                       "FilterProl returned error code: " + std::to_string(ierr));

        // -----------------------------------------------------------------------
        // Pack results into MATLAB output arrays.
        // -----------------------------------------------------------------------
        const std::size_t sz_iat = static_cast<std::size_t>(nn_P + 1);
        const std::size_t sz_nt  = static_cast<std::size_t>(nt_PF);

        // --- nt_PF : scalar double --------------------------------------------
        TypedArray<double> out_nt_PF =
            factory.createScalar<double>(static_cast<double>(nt_PF));

        // --- iat_PF : double, length nn_P+1, 0-based → 1-based ---------------
        TypedArray<double> out_iat_PF = factory.createArray<double>({sz_iat, 1});
        {
            auto it = out_iat_PF.begin();
            for (int i = 0; i <= nn_P; ++i, ++it)
                *it = static_cast<double>(iat_PF_raw[i] + 1);
        }

        // --- ja_PF : double, length nt_PF, 0-based → 1-based -----------------
        TypedArray<double> out_ja_PF = factory.createArray<double>({sz_nt, 1});
        {
            auto it = out_ja_PF.begin();
            for (std::size_t i = 0; i < sz_nt; ++i, ++it)
                *it = static_cast<double>(ja_PF_raw[i] + 1);
        }

        // --- coef_PF : double, length nt_PF -----------------------------------
        TypedArray<double> out_coef_PF = factory.createArray<double>({sz_nt, 1});
        std::copy(coef_PF_raw, coef_PF_raw + sz_nt, out_coef_PF.begin());

        // MallocGuards destruct here — iat/ja/coef_PF_raw freed automatically

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(out_nt_PF);
        outputs[1] = std::move(out_iat_PF);
        outputs[2] = std::move(out_ja_PF);
        outputs[3] = std::move(out_coef_PF);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 10)
            throwError("FilterProl:badInputCount",
                       "Expected 10 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("FilterProl:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Scalar double inputs: np(0), perc(1), tol(2), nn_P(3), nr_TV(7), ntv(8)
        for (std::size_t i : {0u, 1u, 2u, 3u, 7u, 8u})
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("FilterProl:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // int32 array inputs: iat_P(4), ja_P(5)
        for (std::size_t i : {4u, 5u})
            if (inputs[i].getType() != ArrayType::INT32)
                throwError("FilterProl:badArray",
                           "Input argument " + std::to_string(i + 1) +
                           " must be an int32 array.");

        // double array inputs: coef_P(6), TV(9)
        for (std::size_t i : {6u, 9u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("FilterProl:badArray",
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
