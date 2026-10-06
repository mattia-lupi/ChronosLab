#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <stdexcept>
#include "cpt_sam_adaptive_left.h"

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        checkArguments(outputs, inputs);

        // Extract scalar inputs directly
        const int nthread   = static_cast<int>(inputs[7][0]);
        const int n_step    = static_cast<int>(inputs[8][0]);
        const int step_size = static_cast<int>(inputs[9][0]);
        const double eps    = static_cast<double>(inputs[10][0]);

        // Extract typed array views (Zero Copy)
        const matlab::data::TypedArray<double> jatk_in   = inputs[0];
        const matlab::data::TypedArray<double> iak_in    = inputs[1];
        const matlab::data::TypedArray<double> coefk_in  = inputs[2];
        const matlab::data::TypedArray<double> coefkT_in = inputs[3];
        const matlab::data::TypedArray<double> jat0_in   = inputs[4];
        const matlab::data::TypedArray<double> ia0_in    = inputs[5];
        const matlab::data::TypedArray<double> coef0_in  = inputs[6];

        const int nn_A = static_cast<int>(jatk_in.getNumberOfElements()) - 1;

        // --- 1. ZERO-COPY RAW POINTER ACCESS FOR COEFFICIENT INPUTS ---
        const double* ptr_coefk  = &(*coefk_in.begin());
        const double* ptr_coefkT = &(*coefkT_in.begin());
        const double* ptr_coef0  = &(*coef0_in.begin());

        // --- 2. FAST SINGLE-PASS INDEX CONVERSION (double -> 0-based int) ---
        const double* ptr_jatk = &(*jatk_in.begin());
        std::vector<int> iatk(jatk_in.getNumberOfElements());
        for (size_t i = 0; i < iatk.size(); ++i) {
            iatk[i] = static_cast<int>(ptr_jatk[i]) - 1;
        }

        const double* ptr_iak = &(*iak_in.begin());
        std::vector<int> jak(iak_in.getNumberOfElements());
        for (size_t i = 0; i < jak.size(); ++i) {
            jak[i] = static_cast<int>(ptr_iak[i]) - 1;
        }

        const double* ptr_jat0 = &(*jat0_in.begin());
        std::vector<int> iat0(jat0_in.getNumberOfElements());
        for (size_t i = 0; i < iat0.size(); ++i) {
            iat0[i] = static_cast<int>(ptr_jat0[i]) - 1;
        }

        const double* ptr_ia0 = &(*ia0_in.begin());
        std::vector<int> ja0(ia0_in.getNumberOfElements());
        for (size_t i = 0; i < ja0.size(); ++i) {
            ja0[i] = static_cast<int>(ptr_ia0[i]) - 1;
        }

        // Output handles to be populated by the core function
        int* iatN = nullptr;
        int* jaN = nullptr;
        double* coefN = nullptr;
        double avg_resRelNorm = 0.0;

        // Execute core algorithm
        cpt_sam_adaptive_left(iatk.data(), jak.data(),
                               const_cast<double*>(ptr_coefk),
                               const_cast<double*>(ptr_coefkT),
                               iat0.data(), ja0.data(),
                               const_cast<double*>(ptr_coef0),
                               nthread, n_step, step_size, eps, nn_A,
                               iatN, jaN, coefN, avg_resRelNorm);

        const size_t nnz_N = static_cast<size_t>(iatN[nn_A]);

        // --- 3. ZERO-COPY OUTPUT ALLOCATION & RAW POINTER WRITING ---
        matlab::data::ArrayFactory factory;
        matlab::data::TypedArray<double> row_N = factory.createArray<double>({nnz_N, 1});
        matlab::data::TypedArray<double> col_N = factory.createArray<double>({nnz_N, 1});
        matlab::data::TypedArray<double> val_N = factory.createArray<double>({nnz_N, 1});

        // Extract raw C++ double pointers directly from MATLAB output arrays
        double* ptr_row = (nnz_N > 0) ? &(*row_N.begin()) : nullptr;
        double* ptr_col = (nnz_N > 0) ? &(*col_N.begin()) : nullptr;
        double* ptr_val = (nnz_N > 0) ? &(*val_N.begin()) : nullptr;

        size_t idx = 0;
        for (int col = 0; col < nn_A; ++col) {
            const double col_val = static_cast<double>(col + 1);
            const int p_end = iatN[col + 1];
            for (int p = iatN[col]; p < p_end; ++p) {
                ptr_row[idx] = static_cast<double>(jaN[p] + 1);
                ptr_col[idx] = col_val;
                ptr_val[idx] = coefN[p];
                idx++;
            }
        }

        // Return outputs
        outputs[0] = std::move(row_N);
        outputs[1] = std::move(col_N);
        outputs[2] = std::move(val_N);
        outputs[3] = factory.createScalar(avg_resRelNorm);

        // Free dynamically allocated memory from the C function
        delete[] iatN;
        delete[] jaN;
        delete[] coefN;
    }

private:
    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        if (inputs.size() != 11) {
            throw std::invalid_argument("Eleven inputs required.");
        }
        if (outputs.size() > 4) {
            throw std::invalid_argument("Too many output arguments.");
        }
    }
};
