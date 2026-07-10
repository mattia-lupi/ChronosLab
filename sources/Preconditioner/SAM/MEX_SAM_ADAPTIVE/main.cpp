#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>

#include "cpt_sam_adaptive_left.h"

// Helper function to load data from a file into a vector
template <typename T>
std::vector<T> load_vector_from_file(const std::string& filename) {
    std::vector<T> data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        // Soft warning instead of hard error, since fallback is handled
        std::cout << "Notice: Could not open file " << filename << std::endl;
        return data;
    }

    T value;
    while (file >> value) {
        data.push_back(value);
    }
    
    return data;
}

int main() {
    // 1. Load the FIRST matrix from .dat files
    std::vector<double> mat_vals    = load_vector_from_file<double>("coef.dat");
    std::vector<int> mat_cols       = load_vector_from_file<int>("ja.dat");
    std::vector<int> mat_row_ptr    = load_vector_from_file<int>("iat.dat");

    // Check if loading the first matrix was successful
    if (mat_row_ptr.empty() || mat_vals.empty() || mat_cols.empty()) {
        std::cerr << "Failed to load the primary matrix data properly." << std::endl;
        return 1;
    }

    size_t num_rows = mat_row_ptr.size() - 1;

    // 2. Attempt to load the SECOND matrix from .dat files
    std::vector<double> mat2_vals   = load_vector_from_file<double>("coef2.dat");
    std::vector<int> mat2_cols      = load_vector_from_file<int>("ja2.dat");
    std::vector<int> mat2_row_ptr   = load_vector_from_file<int>("iat2.dat");

    // 3. Fallback to Identity Matrix if loading the second matrix failed/is incomplete
    if (mat2_row_ptr.empty() || mat2_vals.empty() || mat2_cols.empty()) {
        std::cout << "--> Second matrix files not found or invalid. Falling back to Identity Matrix.\n";
        
        mat2_vals.assign(num_rows, 1.0);
        mat2_cols.resize(num_rows);
        mat2_row_ptr.resize(num_rows + 1);

        for (size_t i = 0; i < num_rows; ++i) {
            mat2_cols[i] = i;        // Diagonal element column index matches row index
            mat2_row_ptr[i] = i;     // Each row has exactly 1 non-zero element before it
        }
        mat2_row_ptr[num_rows] = num_rows; // Final element is total NNZ
    } else {
        std::cout << "--> Successfully loaded second matrix from files.\n";
    }

    // Variables for the function
    int nthread = 8;
    int n_step = 25;
    int step_size = 1;
    double eps = 1e-5;
    int nn_A = mat_row_ptr.size() - 1;
    int *iat_N = nullptr;
    int *ja_N = nullptr;
    double *coef_N = nullptr;
    double resNorm = 0;

    auto start = std::chrono::steady_clock::now();
    
    // Pass mat2_* which now holds either your loaded data or the fallback identity data
    cpt_sam_adaptive_left(mat_row_ptr.data(), mat_cols.data(), mat_vals.data(),
                          mat2_row_ptr.data(), mat2_cols.data(), mat2_vals.data(),
                          nthread, n_step, step_size, eps, nn_A,
                          iat_N, ja_N, coef_N, resNorm);

    auto end = std::chrono::steady_clock::now();

    // 4. Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
    std::cout << "Time taken: " << duration.count() / 1000.0 << " milliseconds" << std::endl;

    delete[] iat_N;
    delete[] ja_N;
    delete[] coef_N;

    return 0;
}
