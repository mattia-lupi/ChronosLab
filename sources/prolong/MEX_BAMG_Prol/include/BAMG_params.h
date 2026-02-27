/**
 * @file BAMG_params.h
 * @brief This header is used to manage the BAMG parameters.
 * @date March 2021
 * @version 1.1
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

// Levels of verbosity
#define VLEV_NONE 0
#define VLEV_LOW 1
#define VLEV_MEDIUM 2
#define VLEV_HIGH 3

// Factor used to relax the condition on the max allowed row norm
#define RELAX_FAC 1.5

#include "precision.h" // to use: iReg, rExt

/**
 * @brief Structure to store BAMG parameters.
 */
struct BAMG_params {int verbosity; rExt maxrownrm; rExt maxcond; iReg itmax_vol;
                    rExt tol_vol; rExt eps; iReg dist_min; iReg dist_max; iReg mmax; };
