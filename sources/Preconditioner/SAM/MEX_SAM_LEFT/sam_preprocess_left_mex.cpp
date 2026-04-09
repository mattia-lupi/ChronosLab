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

        int max_threads = omp_get_max_threads();
        std::vector<std::vector<int>> thread_r_data(max_threads);
        std::vector<std::vector<int>> thread_r_ptr(max_threads);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();

            int chunk = (n + nthreads - 1) / nthreads;
            int start = std::min(tid * chunk, n);
            int end   = std::min(start + chunk, n);
            int local_n = end - start;

            thread_r_ptr[tid].resize(local_n + 1, 0);
            thread_r_data[tid].reserve(local_n * 32); 

            std::vector<int> marker(n, -1);
            std::vector<int> local_r;
            local_r.reserve(512);

            for (int i = start; i < end; ++i) {
                local_r.clear();
                int s0 = s_ptr[i], s1 = s_ptr[i + 1];

                for (int ji = s0; ji < s1; ++ji) {
                    int j = s_data[ji];
                    int a0 = a_ptr[j], a1 = a_ptr[j + 1];
                    for (int p = a0; p < a1; ++p) {
                        int c = a_col[p];
                        if (marker[c] != i) {
                            marker[c] = i; 
                            local_r.push_back(c);
                        }
                    }
                }
                
                std::sort(local_r.begin(), local_r.end());
                thread_r_data[tid].insert(thread_r_data[tid].end(), local_r.begin(), local_r.end());
                thread_r_ptr[tid][i - start + 1] = static_cast<int>(thread_r_data[tid].size());
            }
        }

        std::vector<int> global_r_ptr(n + 1, 0);
        int total_r_nnz = 0;
        int active_threads = omp_get_max_threads();

        for (int t = 0; t < active_threads; ++t) {
            int chunk = (n + active_threads - 1) / active_threads;
            int start = std::min(t * chunk, n);
            int end   = std::min(start + chunk, n);
            int local_n = end - start;

            for (int i = 0; i < local_n; ++i) {
                global_r_ptr[start + i + 1] = total_r_nnz + thread_r_ptr[t][i + 1];
            }
            total_r_nnz += static_cast<int>(thread_r_data[t].size());
        }

        std::vector<int> global_r_data(total_r_nnz);

        #pragma omp parallel for schedule(static, 1)
        for (int t = 0; t < active_threads; ++t) {
            int chunk = (n + active_threads - 1) / active_threads;
            int start = std::min(t * chunk, n);
            if (start < n) {
                std::copy(thread_r_data[t].begin(), thread_r_data[t].end(), 
                          global_r_data.begin() + global_r_ptr[start]);
            }
        }

        const size_t sz_sp = static_cast<size_t>(n + 1);
        const size_t sz_sd = s_data.size();
        const size_t sz_rp = static_cast<size_t>(n + 1);
        const size_t sz_rd = static_cast<size_t>(total_r_nnz);

        TypedArray<double> out_sptr = factory.createArray<double>({sz_sp, 1});
        TypedArray<double> out_sdat = factory.createArray<double>({sz_sd, 1});
        TypedArray<double> out_rptr = factory.createArray<double>({sz_rp, 1});
        TypedArray<double> out_rdat = factory.createArray<double>({sz_rd, 1});

        for (int i = 0; i <= n; ++i) out_sptr[i] = static_cast<double>(s_ptr[i]);
        for (size_t i = 0; i < sz_sd; ++i) out_sdat[i] = static_cast<double>(s_data[i]);
        for (int i = 0; i <= n; ++i) out_rptr[i] = static_cast<double>(global_r_ptr[i]);
        for (size_t i = 0; i < sz_rd; ++i) out_rdat[i] = static_cast<double>(global_r_data[i]);

        outputs[0] = std::move(out_sptr);
        outputs[1] = std::move(out_sdat);
        outputs[2] = std::move(out_rptr);
        outputs[3] = std::move(out_rdat);
        outputs[4] = factory.createScalar(static_cast<double>(sz_sd));
    }
};
