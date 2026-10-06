/**
 * @file swapi.h
 * @brief This header is used to swaps 2 integer variables.
 * @date July 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

typedef int iReg;

/**
 * @brief Swaps 2 integer variables.
 */
inline void swapi(iReg & i1, iReg & i2){iReg tmp = i1; i1 = i2; i2 = tmp;}

