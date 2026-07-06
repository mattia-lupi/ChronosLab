#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "cpt_sam_adaptive_left.h"

// Helper function to load data from a file into a vector
template <typename T>
std::vector<T> load_vector_from_file(const std::string& filename) {
    std::vector<T> data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    T value;
    while (file >> value) {
        data.push_back(value);
    }
    
    return data;
}

int main() {
    // 1. Load the matrix from .dat files
    std::vector<double> mat_vals    = load_vector_from_file<double>("coef.dat");
    std::vector<int> mat_cols       = load_vector_from_file<int>("ja.dat");
    std::vector<int> mat_row_ptr    = load_vector_from_file<int>("iat.dat");

    // Check if loading was successful
    if (mat_row_ptr.empty()) {
        std::cerr << "Failed to load matrix data properly." << std::endl;
        return 1;
    }

    // 2. Dynamically create the CSR Identity Matrix
    // The number of rows is determined by (row_ptr.size() - 1)
    size_t num_rows = mat_row_ptr.size() - 1;

    std::vector<double> eye_vals(num_rows, 1.0);
    std::vector<int> eye_cols(num_rows);
    std::vector<int> eye_row_ptr(num_rows + 1);

    for (size_t i = 0; i < num_rows; ++i) {
        eye_cols[i] = i;         // Diagonal element column index matches row index
        eye_row_ptr[i] = i;      // Each row has exactly 1 non-zero element before it
    }
    eye_row_ptr[num_rows] = num_rows; // Final element is total NNZ

    // --- Verification Output ---
    // std::cout << "Loaded Matrix Rows: " << num_rows << "\n";
    // std::cout << "Identity Matrix row_ptr sizes match: " << eye_row_ptr.size() << std::endl;

   int nthread = 20;
   int n_step = 25;
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
