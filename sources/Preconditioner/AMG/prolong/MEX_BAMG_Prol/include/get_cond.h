/**
 * @file get_cond.h
 * @brief This header is used to manage the conditioning of an input matrix A.
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

#include "precision.h" // to use: iReg, rExt

#include "move_row_front.h"
#include "mk_HouHolVec.h"
#include "Apply_HouHol_Rot.h"

/**
 * @brief Retrieves rank and conditioning of an input matrix A of size n x m using a rank
 *        revealing QR based on Householder rotations.
 *        The conditioning must be below condmax, in other words condmax is used as an
 *        indicator for numerical singularity.
 */
void get_cond(const rExt condmax, const iReg n, const iReg m, rExt **A,
              iReg &rank, rExt &cond);
