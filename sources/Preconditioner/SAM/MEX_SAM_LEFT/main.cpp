#include <iostream>
#include <vector>
#include <iomanip>
#include "cpt_sam_adaptive_left.h"

int main() {
    // 1. Target Matrix (4x4, 6 Non-zero elements)
    // Row 0: 1.625619, 0.000000, 0.081126, 0.000000 (2 nnz)
    // Row 1: 0.000000, 1.929386, 0.000000, 0.000000 (1 nnz)
    // Row 2: 0.775713, 0.000000, 1.435859, 0.000000 (2 nnz)
    // Row 3: 0.000000, 0.000000, 0.000000, 1.000000 (1 nnz)
    
    std::vector<double> mat_vals = {1.625618560729690,0.081125768865785,1.929385970968730,
                                    0.775712678608402,1.435858588580919,1.000000000000000};

    std::vector<int> mat_cols = {0,2,1,0,2,3};

    std::vector<int> mat_row_ptr = {0,2,3,5,6};

    // 2. Identity Matrix (4x4, 4 Non-zero elements)
    std::vector<double> eye_vals = {1.0, 1.0, 1.0, 1.0};
    std::vector<int> eye_cols = {0, 1, 2, 3};
    std::vector<int> eye_row_ptr = {0, 1, 2, 3, 4};

    cpt_sam_adaptive_left(mat_row_ptr.data(),mat_cols.data(),mat_vals.data(),
                     eye_row_ptr.data(),eye_cols.data(),eye_vals.data(),
                     1,1,1,1e-5,4,
                     nullptr,nullptr,nullptr);

    return 0;
}