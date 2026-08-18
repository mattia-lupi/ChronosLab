#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>
#include <algorithm>

#include "cpt_sam_adaptive_left.h"

// Helper function to load data from a file into a vector
template <typename T>
std::vector<T> load_vector_from_file(const std::string& filename) {
    std::vector<T> data;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cout << "Notice: Could not open file " << filename << std::endl;
        return data;
    }

    T value;
    while (file >> value) {
        data.push_back(value);
    }
    
    return data;
}

// Helper function to save vector data to a file
template <typename T>
bool save_vector_to_file(const std::string& filename, const T* data, size_t size) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return false;
    }
    for (size_t i = 0; i < size; ++i) {
        file << std::setprecision(16) << data[i] << "\n";
    }
    return true;
}

// Helper to print a preview of a vector
template <typename T>
void print_preview(const std::string& name, const T* data, size_t size, size_t max_items = 10) {
    std::cout << name << " [size = " << size << "]: [";
    if (data != nullptr && size > 0) {
        size_t preview_len = std::min(size, max_items);
        for (size_t i = 0; i < preview_len; ++i) {
            std::cout << data[i] << (i + 1 < preview_len ? ", " : "");
        }
        if (size > max_items) {
            std::cout << ", ... (" << size - max_items << " more elements)";
        }
    }
    std::cout << "]\n";
}

// Verification function with separate index and numerical tolerances
bool verify_correctness(const int* iat_N, size_t iat_size,
                        const int* ja_N, size_t ja_size,
                        const double* coef_N, size_t coef_size,
                        int tol_idx = 0,
                        double tol_coef = 1e-6) {
    
    const std::string ref_iat_file  = "iat_N_ref.dat";
    const std::string ref_ja_file   = "ja_N_ref.dat";
    const std::string ref_coef_file = "coef_N_ref.dat";

    std::vector<int> ref_iat     = load_vector_from_file<int>(ref_iat_file);
    std::vector<int> ref_ja      = load_vector_from_file<int>(ref_ja_file);
    std::vector<double> ref_coef = load_vector_from_file<double>(ref_coef_file);

    // Save baseline if reference files do not exist
    if (ref_iat.empty() || ref_ja.empty() || ref_coef.empty()) {
        std::cout << "\n[!] Baseline not found. Saving ground truth files...\n";
        save_vector_to_file(ref_iat_file, iat_N, iat_size);
        save_vector_to_file(ref_ja_file, ja_N, ja_size);
        save_vector_to_file(ref_coef_file, coef_N, coef_size);
        std::cout << "Baseline saved successfully (" 
                  << ref_iat_file << ", " << ref_ja_file << ", " << ref_coef_file << ").\n";
        return true;
    }

    std::cout << "\n=== Verification (tol_idx = " << tol_idx 
              << ", tol_coef = " << tol_coef << ") ===\n";

    // 1. Check iat_N (Row Pointers)
    if (ref_iat.size() != iat_size) {
        std::cerr << "FAILED: iat_N size mismatch. Got " << iat_size 
                  << ", Expected " << ref_iat.size() << "\n";
        return false;
    }
    size_t iat_errors = 0;
    for (size_t i = 0; i < iat_size; ++i) {
        int diff = std::abs(iat_N[i] - ref_iat[i]);
        if (diff > tol_idx) {
            if (iat_errors < 5) {
                std::cerr << "iat_N mismatch at [" << i << "]: Got " 
                          << iat_N[i] << ", Ref " << ref_iat[i] << " (diff: " << diff << ")\n";
            }
            iat_errors++;
        }
    }
    if (iat_errors > 0) {
        std::cerr << "FAILED: " << iat_errors << " iat_N elements exceeded index tolerance.\n";
        return false;
    }
    std::cout << "-> iat_N:  PASSED\n";

    // 2. Check ja_N (Column Indices)
    if (ref_ja.size() != ja_size) {
        std::cerr << "FAILED: ja_N size mismatch. Got " << ja_size 
                  << ", Expected " << ref_ja.size() << "\n";
        return false;
    }
    size_t ja_errors = 0;
    for (size_t i = 0; i < ja_size; ++i) {
        int diff = std::abs(ja_N[i] - ref_ja[i]);
        if (diff > tol_idx) {
            if (ja_errors < 5) {
                std::cerr << "ja_N mismatch at [" << i << "]: Got " 
                          << ja_N[i] << ", Ref " << ref_ja[i] << " (diff: " << diff << ")\n";
            }
            ja_errors++;
        }
    }
    if (ja_errors > 0) {
        std::cerr << "FAILED: " << ja_errors << " ja_N elements exceeded index tolerance.\n";
        return false;
    }
    std::cout << "-> ja_N:   PASSED\n";

    // 3. Check coef_N (Values up to floating-point tolerance)
    if (ref_coef.size() != coef_size) {
        std::cerr << "FAILED: coef_N size mismatch. Got " << coef_size 
                  << ", Expected " << ref_coef.size() << "\n";
        return false;
    }

    double max_abs_diff = 0.0;
    double max_rel_diff = 0.0;
    size_t coef_errors = 0;

    for (size_t i = 0; i < coef_size; ++i) {
        double diff = std::fabs(coef_N[i] - ref_coef[i]);
        double denom = std::max(1.0, std::fabs(ref_coef[i]));
        double rel_diff = diff / denom;

        if (diff > max_abs_diff) max_abs_diff = diff;
        if (rel_diff > max_rel_diff) max_rel_diff = rel_diff;

        if (diff > tol_coef && rel_diff > tol_coef) {
            if (coef_errors < 5) {
                std::cerr << "coef_N mismatch at [" << i << "]: Computed " 
                          << coef_N[i] << ", Ref " << ref_coef[i] 
                          << " (Abs Diff: " << diff << ")\n";
            }
            coef_errors++;
        }
    }

    std::cout << "-> coef_N: PASSED (Max Abs Diff: " << max_abs_diff 
              << ", Max Rel Diff: " << max_rel_diff << ")\n";

    if (coef_errors > 0) {
        std::cerr << "FAILED: " << coef_errors << " coef_N elements exceeded tolerance (" 
                  << tol_coef << ")\n";
        return false;
    }

    std::cout << "=== Overall Status: PASSED ===\n\n";
    return true;
}

int main() {
    // 1. Load the FIRST matrix from .dat files
    std::vector<double> mat_vals    = load_vector_from_file<double>("coef.dat");
    std::vector<double> mat_valsT   = load_vector_from_file<double>("coefT.dat");
    std::vector<int> mat_cols       = load_vector_from_file<int>("ja.dat");
    std::vector<int> mat_row_ptr    = load_vector_from_file<int>("iat.dat");

    if (mat_row_ptr.empty() || mat_vals.empty() || mat_cols.empty()) {
        std::cerr << "Failed to load the primary matrix data properly." << std::endl;
        return 1;
    }

    size_t num_rows = mat_row_ptr.size() - 1;

    // 2. Load the SECOND matrix from .dat files
    std::vector<double> mat2_vals   = load_vector_from_file<double>("coef2.dat");
    std::vector<int> mat2_cols      = load_vector_from_file<int>("ja2.dat");
    std::vector<int> mat2_row_ptr   = load_vector_from_file<int>("iat2.dat");

    // 3. Fallback to Identity Matrix if needed
    if (mat2_row_ptr.empty() || mat2_vals.empty() || mat2_cols.empty()) {
        std::cout << "--> Second matrix files not found. Falling back to Identity Matrix.\n";
        
        mat2_vals.assign(num_rows, 1.0);
        mat2_cols.resize(num_rows);
        mat2_row_ptr.resize(num_rows + 1);

        for (size_t i = 0; i < num_rows; ++i) {
            mat2_cols[i] = i;
            mat2_row_ptr[i] = i;
        }
        mat2_row_ptr[num_rows] = num_rows;
    } else {
        std::cout << "--> Successfully loaded second matrix from files.\n";
    }

    // Function variables
    int nthread = 1;
    int n_step = 5;
    int step_size = 1;
    double eps = 1e-4;
    int nn_A = static_cast<int>(num_rows);
    int *iat_N = nullptr;
    int *ja_N = nullptr;
    double *coef_N = nullptr;
    double resNorm = 0;

    auto start = std::chrono::steady_clock::now();
    
    cpt_sam_adaptive_left(mat_row_ptr.data(), mat_cols.data(), mat_vals.data(), mat_valsT.data(),
                          mat2_row_ptr.data(), mat2_cols.data(), mat2_vals.data(),
                          nthread, n_step, step_size, eps, nn_A,
                          iat_N, ja_N, coef_N, resNorm);

    auto end = std::chrono::steady_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "\nComputation finished.\n";
    std::cout << "Time taken: " << duration.count() / 1000.0 << " ms\n";
    std::cout << "Residual Norm (resNorm): " << resNorm << "\n\n";

    // 4. Print Matrix N details
    size_t iat_size = nn_A + 1;
    size_t nnz_N = (iat_N != nullptr) ? static_cast<size_t>(iat_N[nn_A] - iat_N[0]) : 0;


    // 5. Verification Check
    int tol_idx = 0;           // Tolerance for integer index offsets (iat, ja)
    double tol_coef = 1e-15;    // Tolerance for floating-point matrix values

    bool is_correct = verify_correctness(iat_N, iat_size, 
                                         ja_N, nnz_N, 
                                         coef_N, nnz_N, 
                                         tol_idx, tol_coef);

    // Cleanup dynamic memory
    delete[] iat_N;
    delete[] ja_N;
    delete[] coef_N;

    return is_correct ? 0 : 2;
}
