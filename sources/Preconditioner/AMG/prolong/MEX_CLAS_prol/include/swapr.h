/**
 * @file swapr.h
 * @brief This header is used to swaps 2 real variables.
 * @date August 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

typedef double rExt;

/**
 * @brief Swaps 2 real variables.
 */
inline void swapr(rExt & i1, rExt & i2){rExt tmp = i1; i1 = i2; i2 = tmp;}

