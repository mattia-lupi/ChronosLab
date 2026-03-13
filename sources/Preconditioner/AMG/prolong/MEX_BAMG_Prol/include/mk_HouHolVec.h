/**
 * @file mk_HouHolVec.h
 * @brief This header is used to manage the Householder vector computation.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include <cmath> // to use: sqrt

#include "precision.h" // to use: iReg, rExt
#include "inl_blas1.h" // to use: inl_dnrm2, inl_dcopy

/**
 * @brief Creates from the input vector v, the Householder vector w such that (I - w*w^T) is a
 *        rotation matrix nullifying v from component 1 (included) on.
 * @param [in] n size of v
 * @param [in] v vector used as prototype for the rotation
 * @param [out] w vector for the rotation (same size as v)
 */
void mk_HouHolVec(const iReg n, const rExt *const v, rExt *w);
