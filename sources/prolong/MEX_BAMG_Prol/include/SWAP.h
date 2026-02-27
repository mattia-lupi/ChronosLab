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

#include "precision.h"

/**
 * @brief Swaps 2 TYPE variables.
 */
template <typename TYPE>
inline void SWAP(TYPE & __restrict__ i1, TYPE & __restrict__ i2){TYPE tmp = i1; i1 = i2; i2 = tmp;}

// Instantiate template
#if !IREG_LONG==IEXT_LONG
template void SWAP<iReg>(iReg & __restrict__ , iReg & __restrict__ );
#endif
template void SWAP<iExt>(iExt & __restrict__ , iExt & __restrict__ );
template void SWAP<rExt>(rExt & __restrict__ , rExt & __restrict__ );
