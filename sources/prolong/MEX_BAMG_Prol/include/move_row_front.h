/**
 * @file move_row_front.h
 * @brief This header is used to manage the movement of the j-th row of A.
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
#include "inl_blas1.h" // to use: inl_dnrm2

/**
 * @brief Moves to the first position the j-th row of A whose norm ( A(j,icol:end) )
 *        is maximal. It also returns the value of the norm.
 */
void move_row_front(const iReg jcol, const iReg n, const iReg m, rExt **A,
                    iReg &irow, rExt &maxnorm);
