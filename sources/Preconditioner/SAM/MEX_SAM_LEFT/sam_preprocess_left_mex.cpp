// sam_preprocess_left_mex.cpp
#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include "sam_preprocess_left.h"

#include <vector>
#include <cstdint>

using namespace matlab::data;
using namespace matlab::mex;

class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        ArrayFactory factory;

        if (inputs.size() < 4 || inputs.size() > 5) {
            getEngine()->feval(u"error", 0,
                std::vector<Array>({factory.createScalar(
                    "sam_preprocess_left_mex: requires 4 or 5 inputs (iatk, jak, iats, jas, [num_threads]).")}));
            return;
        }

        TypedArray<double> iatk = std::move(inputs[0]);
        TypedArray<double> jak  = std::move(inputs[1]);
        TypedArray<double> iats = std::move(inputs[2]);
        TypedArray<double> jas  = std::move(inputs[3]);

        int num_threads = 0;
        if (inputs.size() == 5) {
            num_threads = static_cast<int>(inputs[4][0]);
        }

        const int n = static_cast<int>(iatk.getNumberOfElements()) - 1;
        const int nnz_a = static_cast<int>(jak.getNumberOfElements());
        const int nnz_s = static_cast<int>(jas.getNumberOfElements());

        // Convert 1-based MATLAB double arrays into 0-based C++ int vectors
        std::vector<int> a_ptr(n + 1);
        std::vector<int> a_col(nnz_a);
        std::vector<int> s_ptr(n + 1);
        std::vector<int> s_data(nnz_s);

        for (int i = 0; i <= n; ++i) a_ptr[i] = static_cast<int>(iatk[i]) - 1;
        for (int p = 0; p < nnz_a; ++p) a_col[p] = static_cast<int>(jak[p]) - 1;

        for (int i = 0; i <= n; ++i) s_ptr[i] = static_cast<int>(iats[i]) - 1;
        for (int p = 0; p < nnz_s; ++p) s_data[p] = static_cast<int>(jas[p]) - 1;

        std::vector<int> r_ptr;
        std::vector<int> r_data;

        // Call Standalone C++ Function
        sam_preprocess_left(n, a_ptr.data(), a_col.data(),
                            s_ptr.data(), s_data.data(),
                            r_ptr, r_data, num_threads);

        // Package outputs into MATLAB TypedArrays (int32_t)
        const size_t sz_sp = static_cast<size_t>(n + 1);
        const size_t sz_sd = static_cast<size_t>(nnz_s);
        const size_t sz_rp = static_cast<size_t>(n + 1);
        const size_t sz_rd = r_data.size();

        TypedArray<int32_t> out_sptr = factory.createArray<int32_t>({sz_sp, 1});
        TypedArray<int32_t> out_sdat = factory.createArray<int32_t>({sz_sd, 1});
        TypedArray<int32_t> out_rptr = factory.createArray<int32_t>({sz_rp, 1});
        TypedArray<int32_t> out_rdat = factory.createArray<int32_t>({sz_rd, 1});

        for (size_t i = 0; i < sz_sp; ++i) out_sptr[i] = static_cast<int32_t>(s_ptr[i]);
        for (size_t i = 0; i < sz_sd; ++i) out_sdat[i] = static_cast<int32_t>(s_data[i]);
        for (size_t i = 0; i < sz_rp; ++i) out_rptr[i] = static_cast<int32_t>(r_ptr[i]);
        for (size_t i = 0; i < sz_rd; ++i) out_rdat[i] = static_cast<int32_t>(r_data[i]);

        outputs[0] = std::move(out_sptr);
        outputs[1] = std::move(out_sdat);
        outputs[2] = std::move(out_rptr);
        outputs[3] = std::move(out_rdat);
        outputs[4] = factory.createScalar(static_cast<int32_t>(sz_sd));
    }
};
