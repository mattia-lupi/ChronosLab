/**
 * @file FilterComp_stripe.h
 * @brief This header is used to filter a stripe of the operator. 
 * @date April 2020
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "lapacke.h"
#include <algorithm> // tu use: max
using namespace std;

#include "abs_ri_heapsort.h"
#include "bin_search.h"
#include "ir_heapsort.h"

/**
 * @brief Inner part of FilterComp that is used to distributed work among threads. 
 */
int FilterComp_stripe(const double tau, const int shift, const int nrows, 
                      const int nrows_A, int* iat_A, int* ja_A, double* coef_A, 
                      const int nt_patt, const int* iat_patt, const int* ja_patt, 
                      const int ntv, const double *const *TV, int &nt_AC_loc, 
                      int *&iat_out, int *&ja_out, double *&coef_out);
