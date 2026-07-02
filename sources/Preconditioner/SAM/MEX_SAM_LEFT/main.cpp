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
   
   std::vector<double> mat_vals = {1.625619,0.081126,0.156789,
                                   1.929386,0.14568, 
                                   0.775713,1.435859,
                                   -0.45245,1.0000,
                                   0.145680,1.0000};

   std::vector<int> mat_cols = {0, 2, 3, 1, 4, 0, 2, 0, 3, 1, 4};

   std::vector<int> mat_row_ptr = {0, 3, 5, 7, 9, 11};

   // 2. Identity Matrix (4x4, 4 Non-zero elements)
   std::vector<double> eye_vals = {1.0, 1.0, 1.0, 1.0,1,0};
   std::vector<int> eye_cols = {0, 1, 2, 3, 4};
   std::vector<int> eye_row_ptr = {0, 1, 2, 3, 4, 5};

   int nthread = 5;
   int n_step = 5;
   int step_size = 1;
   double eps = 1e-5;
   int nn_A = mat_row_ptr.size() - 1;
   int *iat_N = nullptr;
   int *ja_N = nullptr;
   double *coef_N = nullptr;
   double resNorm = 0;

   cpt_sam_adaptive_left(mat_row_ptr.data(),mat_cols.data(),mat_vals.data(),
                    eye_row_ptr.data(),eye_cols.data(),eye_vals.data(),
                    nthread,n_step,step_size,eps,nn_A,
                    iat_N,ja_N,coef_N,resNorm);

   return 0;
}