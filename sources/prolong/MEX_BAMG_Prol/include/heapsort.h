/**
 * @file heapsort.h
 * @brief This header is used to manage the sorting of an array.
 * @date July 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "precision.h"
#include "SWAP.h"

/**
 * @brief Sorts an array x1 in such a way that x1(i) <= x1(i+1).
 * @param [inout] x1 array to be sorted.
 * @param [in] n number of components of x1.
 */
template <typename Tx, typename Tn>
void heapsort(Tx* __restrict__ x1, const Tn n);
