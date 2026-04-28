/**
 * @file swap.h
 * @brief This template is used to swaps 2 TYPE variables.
 * @date July 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

/**
 * @brief Swaps 2 TYPE variables.
 */
template <typename TYPE>
inline void SWAP(TYPE &i1, TYPE &i2){TYPE tmp = i1; i1 = i2; i2 = tmp;}

// Instantiate template
template void SWAP<int>(int &, int &);
