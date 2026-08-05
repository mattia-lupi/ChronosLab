// sam_compute_left_mex.cpp
#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include "sam_compute_left.h"

#include <vector>
#include <cstdint>

using namespace matlab::data;
using namespace matlab::mex;

static_assert(sizeof(int) == sizeof(int32_t), "int size must match int32_t");

class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        ArrayFactory factory;

        if (inputs.size() < 11 || inputs.size() > 12) {
            getEngine()->feval(u"error", 0,
                std::vector<Array>({factory.createScalar(
                    "sam_compute_left_mex: requires 11 or 12 inputs (iatk, jak, coefk, iat0, ja0, coef0, s_ptr, s_data, r_ptr, r_data, nnz_total, [num_threads]).")}));
            return;
        }

        TypedArray<double> iatk  = std::move(inputs[0]);
        TypedArray<double> jak   = std::move(inputs[1]);
        TypedArray<double> coefk = std::move(inputs[2]);

        TypedArray<double> iat0  = std::move(inputs[3]);
        TypedArray<double> ja0   = std::move(inputs[4]);
        TypedArray<double> coef0 = std::move(inputs[5]);

        TypedArray<int32_t> s_ptr_arr = std::move(inputs[6]);
        TypedArray<int32_t> s_dat_arr = std::move(inputs[7]);
        TypedArray<int32_t> r_ptr_arr = std::move(inputs[8]);
        TypedArray<int32_t> r_dat_arr = std::move(inputs[9]);

        int total_nnz;
        if (inputs[10].getType() == ArrayType::INT32) {
            total_nnz = static_cast<TypedArray<int32_t>>(inputs[10])[0];
        } else {
            total_nnz = static_cast<int>(static_cast<TypedArray<double>>(inputs[10])[0]);
        }

        int num_threads = 0;
        if (inputs.size() == 12) {
            if (inputs[11].getType() == ArrayType::INT32) {
                num_threads = static_cast<TypedArray<int32_t>>(inputs[11])[0];
            } else {
                num_threads = static_cast<int>(static_cast<TypedArray<double>>(inputs[11])[0]);
            }
        }

        const int n = static_cast<int>(iatk.getNumberOfElements()) - 1;
        const int nnz_ak = static_cast<int>(jak.getNumberOfElements());
        const int nnz_a0 = static_cast<int>(ja0.getNumberOfElements());

        // Convert 1-based MATLAB CSR arrays into 0-based C++ int vectors
        std::vector<int> ak_ptr(n + 1);
        std::vector<int> ak_col(nnz_ak);
        std::vector<int> a0_ptr(n + 1);
        std::vector<int> a0_col(nnz_a0);

        for (int i = 0; i <= n; ++i) ak_ptr[i] = static_cast<int>(iatk[i]) - 1;
        for (int p = 0; p < nnz_ak; ++p) ak_col[p] = static_cast<int>(jak[p]) - 1;

        for (int i = 0; i <= n; ++i) a0_ptr[i] = static_cast<int>(iat0[i]) - 1;
        for (int p = 0; p < nnz_a0; ++p) a0_col[p] = static_cast<int>(ja0[p]) - 1;

        // Zero-copy pointers for values and pre-converted 0-based arrays
        const double* ak_val = &(*coefk.begin());
        const double* a0_val = &(*coef0.begin());

        const int32_t* s_ptr  = &(*s_ptr_arr.begin());
        const int32_t* s_data = &(*s_dat_arr.begin());
        const int32_t* r_ptr  = &(*r_ptr_arr.begin());
        const int32_t* r_data = &(*r_dat_arr.begin());

        const size_t sz = static_cast<size_t>(total_nnz);
        TypedArray<double> out_row = factory.createArray<double>({sz, 1});
        TypedArray<double> out_col = factory.createArray<double>({sz, 1});
        TypedArray<double> out_val = factory.createArray<double>({sz, 1});

        double* p_row = &(*out_row.begin());
        double* p_col = &(*out_col.begin());
        double* p_val = &(*out_val.begin());

        // Fill 1-based row/col outputs for MATLAB sparse matrix construction
        for (int i = 0; i < n; ++i) {
            for (int p = s_ptr[i]; p < s_ptr[i + 1]; ++p) {
                p_row[p] = static_cast<double>(i + 1);
                p_col[p] = static_cast<double>(s_data[p] + 1);
            }
        }

        // Call Standalone C++ Compute Function
        sam_compute_left(n,
                         ak_ptr.data(), ak_col.data(), ak_val,
                         a0_ptr.data(), a0_col.data(), a0_val,
                         reinterpret_cast<const int*>(s_ptr),
                         reinterpret_cast<const int*>(s_data),
                         reinterpret_cast<const int*>(r_ptr),
                         reinterpret_cast<const int*>(r_data),
                         p_val,
                         num_threads);

        outputs[0] = std::move(out_row);
        outputs[1] = std::move(out_col);
        outputs[2] = std::move(out_val);
    }
};
