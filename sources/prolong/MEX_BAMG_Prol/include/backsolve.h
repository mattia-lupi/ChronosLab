/**
 * @file backsolve.h
 * @brief This header is used to manage a back substitution with a upper triangular matrix.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "precision.h" // to use: iReg, rExt

/**
 * @brief Performs a back substitution with the upper triangular (non unitary) matrix R
 *        solving Z_in = Z_out*R
 * @param [in] n number of rows/columns of R and columns of Z
 * @param [in] m number of rows of Z
 * @param [in] R (non unitary) upper triangular matrix
 * @param [inout] Z matrix that is backsubstituted in place
 */
void backsolve(const iReg n,const iReg m, const rExt *const *const R, rExt **Z);
