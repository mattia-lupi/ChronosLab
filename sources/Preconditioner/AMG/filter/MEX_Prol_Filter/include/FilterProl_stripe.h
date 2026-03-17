/**
 * @file FilterProl_stripe.h
 * @brief This header is used to filter a stripe of the prolongation operator. 
 * @date May 2020
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "lapacke.h"
#include <algorithm> // tu use: max,min

#include "abs_norm.h"
#include "ir_heapsort.h"
#include "abs_ri_heapsort.h"

int FilterProl_stripe(const double perc, const double tol, const int firstrow,
                      const int nrows, const int nn_P, int *iat_P, int *ja_P,
                      double *coef_P, const int ntv, const double *const *TV,
                      int &nt_PF_loc, int *&iat_PF, int *&ja_PF, double *&coef_PF);
