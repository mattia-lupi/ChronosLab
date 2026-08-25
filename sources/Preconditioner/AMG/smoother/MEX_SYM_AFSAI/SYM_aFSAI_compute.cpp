#include <cstdint>   // int64_t, int32_t
#include "mex.hpp"
#include "mexAdapter.hpp"

#include <string>
#include <cassert>
#include <algorithm>

#include "compute_local_fsai.h"

//----------------------------------------------------------------------------------------
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [nterm_G, iat_G, ja_G, coef_G] =
//       compute_local_fsai_wrap(nthread, n_step, step_size, tau, eps,
//                               nrows, nrows_M, nterm_M,
//                               iat_M, ja_M, coef_M)


static_assert(sizeof(iExt) == sizeof(int64_t),
              "iExt (long int) must be 64 bits on this platform. "
              "Update the typedef or the TypedArray<int64_t> mapping.");

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    // [FIX-A] Route diagnostic output through MATLAB's fprintf
    void mprint(const std::string& msg)
    {
        getEngine()->feval(u"fprintf", 0,
            std::vector<Array>{ factory.createCharArray(msg) });
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Read input
        // -----------------------------------------------------------------------

        const iReg nthread   = static_cast<iReg>(TypedArray<double>(inputs[0])[0]);
        const iReg n_step    = static_cast<iReg>(TypedArray<double>(inputs[1])[0]);
        const iReg step_size = static_cast<iReg>(TypedArray<double>(inputs[2])[0]);
        const rExt tau       = static_cast<rExt>(TypedArray<double>(inputs[3])[0]);
        const rExt eps       = static_cast<rExt>(TypedArray<double>(inputs[4])[0]);
        const iReg nrows     = static_cast<iReg>(TypedArray<double>(inputs[5])[0]);
        const iReg nrows_M   = static_cast<iReg>(TypedArray<double>(inputs[6])[0]);
        const iExt nterm_M   = static_cast<iExt>(TypedArray<double>(inputs[7])[0]);

        const TypedArray<int64_t> iat_M_arr  = inputs[8];
        const TypedArray<int32_t> ja_M_arr   = inputs[9];
        const TypedArray<double>  coef_M_arr = inputs[10];

        std::vector<iExt> iat_M_vec (iat_M_arr.begin(),  iat_M_arr.end());
        std::vector<iReg> ja_M_vec  (ja_M_arr.begin(),   ja_M_arr.end());
        std::vector<rExt> coef_M_vec(coef_M_arr.begin(), coef_M_arr.end());

        const iExt kmax    = 1 + static_cast<iExt>(n_step) * static_cast<iExt>(step_size);
        const iExt nzmax_G = static_cast<iExt>(nrows) * kmax;

        iExt              nterm_G_val = 0;                      // scalar output
        std::vector<iExt> iat_G_vec(static_cast<std::size_t>(nrows) + 1, 0);
        std::vector<iReg> ja_G_vec (static_cast<std::size_t>(nzmax_G),   0);
        std::vector<rExt> coef_G_vec(static_cast<std::size_t>(nzmax_G),  0.0);

        // -----------------------------------------------------------------------
        // Call the computational kernel
        // -----------------------------------------------------------------------

        compute_local_fsai(nthread, n_step, step_size, tau, eps,
                           nrows, nrows_M, nterm_M,
                           iat_M_vec.data(),
                           ja_M_vec.data(),
                           coef_M_vec.data(),
                           &nterm_G_val,
                           iat_G_vec.data(),
                           ja_G_vec.data(),
                           coef_G_vec.data());

        const std::size_t sz_rows1 = static_cast<std::size_t>(nrows) + 1;
        const std::size_t sz_nt    = static_cast<std::size_t>(nterm_G_val);

        // --- nterm_G : int64 scalar -------------------------------------------
        TypedArray<int64_t> out_nterm_G =
            factory.createArray<int64_t>({1, 1});
        out_nterm_G[0] = static_cast<int64_t>(nterm_G_val);

        // --- iat_G : int64, length nrows+1 ------------------------------------
        TypedArray<int64_t> out_iat_G =
            factory.createArray<int64_t>({1, sz_rows1});
        {
            auto it = out_iat_G.begin();
            for (std::size_t i = 0; i < sz_rows1; ++i, ++it)
                *it = static_cast<int64_t>(iat_G_vec[i]);
        }

        // --- ja_G : int32, length nterm_G_val (trimmed from oversized work) ---
        TypedArray<int32_t> out_ja_G =
            factory.createArray<int32_t>({1, sz_nt});
        {
            auto it = out_ja_G.begin();
            for (std::size_t i = 0; i < sz_nt; ++i, ++it)
                *it = static_cast<int32_t>(ja_G_vec[i]);
        }

        // --- coef_G : double, length nterm_G_val ------------------------------
        TypedArray<double> out_coef_G =
            factory.createArray<double>({1, sz_nt});
        std::copy(coef_G_vec.begin(),
                  coef_G_vec.begin() + static_cast<std::ptrdiff_t>(sz_nt),
                  out_coef_G.begin());

        // Work vectors (ja_G_vec, coef_G_vec) destruct here automatically —
        // equivalent to mxDestroyArray(ja_W) / mxDestroyArray(coef_W)

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(out_nterm_G);
        outputs[1] = std::move(out_iat_G);
        outputs[2] = std::move(out_ja_G);
        outputs[3] = std::move(out_coef_G);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 11)
            throwError("LocalFSAI:badInputCount",
                       "Expected 11 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("LocalFSAI:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Scalar double inputs: nthread(0)..nterm_M(7)
        for (std::size_t i = 0; i < 8; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("LocalFSAI:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // int64 array input: iat_M(8)
        if (inputs[8].getType() != ArrayType::INT64)
            throwError("LocalFSAI:badArray",
                       "Input argument 9 (iat_M) must be an int64 array.");

        // int32 array input: ja_M(9)
        if (inputs[9].getType() != ArrayType::INT32)
            throwError("LocalFSAI:badArray",
                       "Input argument 10 (ja_M) must be an int32 array.");

        // double array input: coef_M(10)
        if (inputs[10].getType() != ArrayType::DOUBLE)
            throwError("LocalFSAI:badArray",
                       "Input argument 11 (coef_M) must be a double array.");
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

