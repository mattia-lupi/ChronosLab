/**
 * @file maxVol_inner.h
 * @brief This header is used to manage the inner function of maxVol.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include <cmath> // to use: fabs

#include "precision.h"    // to use: iReg, rExt
#include "linsol_error.h" // to throw errors

#include "SWAP.h"        // to swap two integers

/**
 * @brief Inner function of maxVol, list(n+r) list of best columns whose first r are the best.
 */
void maxVol_inner(const iReg itmax, const rExt delta, const iReg n, const iReg r, iReg *list, rExt** Z);
