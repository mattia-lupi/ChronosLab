#ifndef SAM_COMPUTE_LEFT_H
#define SAM_COMPUTE_LEFT_H

#include <vector>

// Solves local micro least-squares systems to compute SAM output values.
void sam_compute_left(int n,
                      const int* ak_ptr, const int* ak_col, const double* ak_val,
                      const int* a0_ptr, const int* a0_col, const double* a0_val,
                      const int* s_ptr,  const int* s_data,
                      const int* r_ptr,  const int* r_data,
                      double* out_val,
                      int num_threads = 0);

#endif // SAM_COMPUTE_LEFT_H
