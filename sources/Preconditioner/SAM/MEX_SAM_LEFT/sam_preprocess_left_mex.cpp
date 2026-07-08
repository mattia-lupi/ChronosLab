// sam_preprocess_left_mex.cpp
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <vector>
#include <algorithm>
#include <omp.h>

using namespace matlab::data;
using namespace matlab::mex;

static void buildCSR(const TypedArray<double>& iat,
                     const TypedArray<double>& ja,
                     int n,
                     std::vector<int>& ptr,
                     std::vector<int>& col)
{
    int nnz = static_cast<int>(ja.getNumberOfElements());
    ptr.resize(n + 1);
    col.resize(nnz);
    for (int i = 0; i <= n; ++i) ptr[i] = static_cast<int>(iat[i]) - 1;
    for (int p = 0; p < nnz; ++p) col[p] = static_cast<int>(ja[p]) - 1;
}

class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        ArrayFactory factory;

        if (inputs.size() != 4) {
            getEngine()->feval(u"error", 0,
                std::vector<Array>({factory.createScalar(
                    "sam_preprocess_left_mex: requires iatk, jak, iats, jas.")}));
            return;
        }

        TypedArray<double> iatk = std::move(inputs[0]);
        TypedArray<double> jak  = std::move(inputs[1]);
        TypedArray<double> iats = std::move(inputs[2]);
        TypedArray<double> jas  = std::move(inputs[3]);

        const int n = static_cast<int>(iatk.getNumberOfElements()) - 1;

        std::vector<int> a_ptr, a_col, s_ptr, s_data;
        buildCSR(iatk, jak, n, a_ptr, a_col);
        buildCSR(iats, jas, n, s_ptr, s_data);

        std::vector<int> r_ptr(n + 1, 0);

        #pragma omp parallel
        {
            std::vector<int> marker(n, -1);
            #pragma omp for schedule(guided)
            for (int i = 0; i < n; ++i) {
                int count = 0;
                int s0 = s_ptr[i], s1 = s_ptr[i + 1];

                for (int ji = s0; ji < s1; ++ji) {
                    int j = s_data[ji];
                    int a0 = a_ptr[j], a1 = a_ptr[j + 1];
                    for (int p = a0; p < a1; ++p) {
                        int c = a_col[p];
                        if (marker[c] != i) {
                            marker[c] = i;
                            count++;
                        }
                    }
                }
                r_ptr[i + 1] = count;
            }
        }

        for (int i = 0; i < n; ++i) {
            r_ptr[i + 1] += r_ptr[i];
        }
        
        int total_r_nnz = r_ptr[n];
        std::vector<int> r_data(total_r_nnz);

        #pragma omp parallel
        {
            std::vector<int> marker(n, -1);
            #pragma omp for schedule(guided)
            for (int i = 0; i < n; ++i) {
                int s0 = s_ptr[i], s1 = s_ptr[i + 1];
                int offset = r_ptr[i];
                int count = 0;

                for (int ji = s0; ji < s1; ++ji) {
                    int j = s_data[ji];
                    int a0 = a_ptr[j], a1 = a_ptr[j + 1];
                    for (int p = a0; p < a1; ++p) {
                        int c = a_col[p];
                        if (marker[c] != i) {
                            marker[c] = i;
                            r_data[offset + count] = c;
                            count++;
                        }
                    }
                }
                std::sort(r_data.begin() + offset, r_data.begin() + offset + count);
            }
        }

        const size_t sz_sp = static_cast<size_t>(n + 1);
        const size_t sz_sd = s_data.size();
        const size_t sz_rp = static_cast<size_t>(n + 1);
        const size_t sz_rd = static_cast<size_t>(total_r_nnz);

        // Type downcasting to int32_t to halve memory bandwidth
        TypedArray<int32_t> out_sptr = factory.createArray<int32_t>({sz_sp, 1});
        TypedArray<int32_t> out_sdat = factory.createArray<int32_t>({sz_sd, 1});
        TypedArray<int32_t> out_rptr = factory.createArray<int32_t>({sz_rp, 1});
        TypedArray<int32_t> out_rdat = factory.createArray<int32_t>({sz_rd, 1});

        for (int i = 0; i <= n; ++i) out_sptr[i] = static_cast<int32_t>(s_ptr[i]);
        for (size_t i = 0; i < sz_sd; ++i) out_sdat[i] = static_cast<int32_t>(s_data[i]);
        for (int i = 0; i <= n; ++i) out_rptr[i] = static_cast<int32_t>(r_ptr[i]);
        for (size_t i = 0; i < sz_rd; ++i) out_rdat[i] = static_cast<int32_t>(r_data[i]);

        outputs[0] = std::move(out_sptr);
        outputs[1] = std::move(out_sdat);
        outputs[2] = std::move(out_rptr);
        outputs[3] = std::move(out_rdat);
        outputs[4] = factory.createScalar(static_cast<int32_t>(sz_sd));
    }
};