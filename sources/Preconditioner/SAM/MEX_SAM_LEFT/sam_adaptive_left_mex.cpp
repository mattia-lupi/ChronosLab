// sam_adaptive_left_mex.cpp
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include "cpt_sam_adaptive_left.h"
#include <vector>


class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        ArrayFactory factory;

        if (inputs.size() != 11) {
            getEngine()->feval(u"error", 0, std::vector<Array>({factory.createScalar("sam_compute_left_mex: requires iatk, jak, coefk, iat0, ja0, coef0, s_ptr, s_data, r_ptr, r_data, nnz_total.")}));
            return;
        }

        TypedArray<double> iatk  = std::move(inputs[0]); TypedArray<double> jak   = std::move(inputs[1]); TypedArray<double> coefk = std::move(inputs[2]);
        TypedArray<double> iat0  = std::move(inputs[3]); TypedArray<double> ja0   = std::move(inputs[4]); TypedArray<double> coef0 = std::move(inputs[5]);
        
        const int nthread   = static_cast<int>(TypedArray<double>(inputs[6])[0]);
        const int n_step    = static_cast<int>(TypedArray<double>(inputs[7])[0]);
        const int step_size = static_cast<int>(TypedArray<double>(inputs[8])[0]);
        const double eps    = static_cast<double>(TypedArray<double>(inputs[9])[0]);// Tolerance of SAM
        const long int nn_A = static_cast<int>(TypedArray<double>(inputs[10])[0]);// size of A
        
        long int max_nnz = n_step * step_size * nn_A;

        // Allocate the sparse matrix for the result (max size possible)
        std::vector<int> iatN_vec(nn_A+1);
        std::vector<long int> jaN_vec(max_nnz);
        std::vector<double> coefN_vec(max_nnz);

        // Get handles
        int *iatN = iatN_vec.data();
        long int *jaN = jaN_vec.data();
        double *coefN = coefN_vec.data();

        cpt_sam_adaptive_left(iatk,jak,coefk,iat0,ja0,coef0,nthread,n_step,step_size,eps,nn_A,iatN,jaN,coefN);



        
    }
};


