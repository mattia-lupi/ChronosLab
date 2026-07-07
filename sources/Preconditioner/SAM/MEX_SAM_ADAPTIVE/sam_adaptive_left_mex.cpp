#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <cstdint>
#include "cpt_sam_adaptive_left.h"


class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        checkArguments(outputs, inputs);

        // Extract scalar inputs
        const int nthread   = static_cast<int>(inputs[6][0]);
        const int n_step    = static_cast<int>(inputs[7][0]);
        const int step_size = static_cast<int>(inputs[8][0]);
        const double eps    = static_cast<double>(inputs[9][0]);

        // Extract array inputs
        matlab::data::TypedArray<double> jatk_in  = std::move(inputs[0]);
        matlab::data::TypedArray<double> iak_in   = std::move(inputs[1]);
        matlab::data::TypedArray<double> coefk_in = std::move(inputs[2]);
        
        matlab::data::TypedArray<double> jat0_in  = std::move(inputs[3]);
        matlab::data::TypedArray<double> ia0_in   = std::move(inputs[4]);
        matlab::data::TypedArray<double> coef0_in = std::move(inputs[5]);

        const int nn_A = static_cast<int>(jatk_in.getNumberOfElements()) - 1;

        // Convert 1-based MATLAB indexing to 0-based C++ indexing for CSC structures
        std::vector<int> iatk(jatk_in.begin(), jatk_in.end());
        std::vector<int> jak(iak_in.begin(), iak_in.end());
        for (auto& val : iatk) val--;
        for (auto& val : jak)  val--;

        std::vector<int> iat0(jat0_in.begin(), jat0_in.end());
        std::vector<int> ja0(ia0_in.begin(), ia0_in.end());
        for (auto& val : iat0) val--;
        for (auto& val : ja0)  val--;

        // Copy data to vector to obtain valid contiguous raw pointers, avoiding ArrayElementTypedRef errors
        std::vector<double> coefk(coefk_in.begin(), coefk_in.end());
        std::vector<double> coef0(coef0_in.begin(), coef0_in.end());

        // Output references to be populated by the function
        int* iatN = nullptr;
        int* jaN = nullptr;
        double* coefN = nullptr;
        double avg_resRelNorm = 0.0;

        // Execute core algorithm
        cpt_sam_adaptive_left(iatk.data(), jak.data(), coefk.data(),
                               iat0.data(), ja0.data(), coef0.data(),
                               nthread, n_step, step_size, eps, nn_A,
                               iatN, jaN, coefN, avg_resRelNorm);

        // Determine number of non-zero elements from the generated CSC column pointer array
        const int nnz_N = iatN[nn_A];

        // Allocate MATLAB output arrays using correct API factory method
        matlab::data::ArrayFactory factory;
        matlab::data::TypedArray<double> row_N = factory.createArray<double>({static_cast<size_t>(nnz_N), 1});
        matlab::data::TypedArray<double> col_N = factory.createArray<double>({static_cast<size_t>(nnz_N), 1});
        matlab::data::TypedArray<double> val_N = factory.createArray<double>({static_cast<size_t>(nnz_N), 1});

        // Uncompress CSC into COO format and convert to 1-based indexing
        int idx = 0;
        for (int col = 0; col < nn_A; ++col) {
            for (int p = iatN[col]; p < iatN[col + 1]; ++p) {
                row_N[idx] = static_cast<double>(jaN[p] + 1);
                col_N[idx] = static_cast<double>(col + 1);
                val_N[idx] = coefN[p];
                idx++;
            }
        }

        // Assign outputs
        outputs[0] = row_N;
        outputs[1] = col_N;
        outputs[2] = val_N;
        outputs[3] = factory.createScalar(avg_resRelNorm);

        // Free dynamically allocated memory from the C++ function
        delete[] iatN;
        delete[] jaN;
        delete[] coefN;
    }

private:
    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        matlab::data::ArrayFactory factory;
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        if (inputs.size() != 10) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({
                factory.createScalar("Ten inputs required.") }));
        }
        if (outputs.size() > 4) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({
                factory.createScalar("Too many output arguments.") }));
        }
    }
};