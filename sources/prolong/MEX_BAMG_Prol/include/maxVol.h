/**
 * @file maxVol.h
 * @brief This header is used to manage the application of the max volume algorithm.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include <cmath>     // to use: fabs
#include <algorithm> // to use: min

#include "precision.h"    // to use: iReg, rExt
#include "linsol_error.h" // to throw errors

#include "move_row_front.h"
#include "mk_HouHolVec.h"
#include "Apply_HouHol_Rot.h"
#include "backsolve.h"
#include "SWAP.h"
#include "maxVol_inner.h"

/**
 * @brief Applies the max volume algorithm to retrieve the best possible basis among the
 *        rows of the input matrix A[n,m].
 * @details The routine returns the list of the best rows along with the rank of the matrix.
 *          There are two parameters controlling the algorithm itmax and delta
 *          For reference see the paper "How to find a good submatrix" by S. A. Goreinov,
 *          I. V. Oseledets, D. V. Savostyanov, E. E. Tyrtyshnikov, N. L. Zamarashkin
 * @param [in] mmax max number of members in the basis
 * @param [in] condmax used to determine the rank of A (the rank of A must give
 *                     a conditioning below condmax)
 * @param [in] itmax max number of iterations of the max volume algoritgm
 * @param [in] delta the algorithm stops when zij < 1+delta
 * @param [in] n number of rows of A
 * @param [in] m number of columns of A
 * @param [inout] A the matrix whose rows are analyzed to find the best basis
 * @param [out] rank rank of A
 * @param [out] list list of the rows representing the best basis
 */
void maxVol(const iReg mmax, const rExt condmax, const iReg itmax, const rExt delta,
            const iReg n, const iReg m, rExt **A, iReg &rank, iReg *list);


