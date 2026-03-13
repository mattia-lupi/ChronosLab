/**
 * @file Apply_HouHol_Rot.h
 * @brief This header is used to manage the appllication of the Householder rotation.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "precision.h" // to use: iReg, rExt
#include "inl_blas1.h" // to use: inl_dnrm2, inl_dcopy

/**
 * @brief Applies one Householder rotation to a vector.
 */
void Apply_HouHol_Rot(iReg n, const rExt *const w, rExt *v);
