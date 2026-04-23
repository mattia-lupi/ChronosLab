/**
 * @file ir_heapsort.h
 * @brief This header is used to sort an integer array (and a real one in the same way).
 * @date August 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

typedef int iReg;
typedef double rExt;
#include "swapi.h"
#include "swapr.h"

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

/**
 * @brief Sorts an integer array x1 in such a way that x1(i) <= x1(i+1). 
 *        and a real array x2 in the same way.
 * @param [inout] x1 array of integers to be sorted.
 * @param [inout] x2 array of reals to be sorted.
 * @param [in] n number of components of x1.
 */
void ir_heapsort(iReg* RESTRICT x1, rExt* RESTRICT x2, const iReg n);

