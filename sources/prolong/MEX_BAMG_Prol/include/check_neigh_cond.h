/**
 * @file check_neigh_cond.h
 * @brief This header is used to manage the conditioning of the neighbouring nodes.
 * @date March 2021
 * @version 1.1
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

////////////////
#include <omp.h>
////////////////

#include "precision.h"    // to use: iReg, rExt
#include "linsol_error.h" // to throw errors

#include "add_new_neighs.h"
#include "get_cond.h"

/**
 * @brief Evaluates the conditioning of the test space vectors assiciated to the
 *        the neighbouring nodes of inod
 * @param [in] maxcond used as an indicator for numerical singularity
 * @param [in] inod central node considered to evaluate local conditioning
 * @param [in] ntv number of test vectors
 * @param [in] nn_S number of rows of the SoC matrix
 * @param [in] iat_S array of pointers of the SoC matrix
 * @param [in] ja_S array of column indices of the SoC matrix
 * @param [in] TV test vectors matrix of size [nn_S,ntv]
 * @param [out] TVcomp compressed matrix of test vectors
 * @param [out] neigh list of neighbouring nodes
 * @param [out] WI indicator for the presence of a neighbor
 * @param [out] local_cond local conditioning
 */
void check_neigh_cond(const rExt maxcond,const iReg inod, const iReg ntv, const iReg nn_S,
                      const iExt *const iat_S, const iReg *const ja_S,
                      const rExt *const *const TV,
                      rExt **TVcomp, iReg *neigh, iReg *WI, rExt &local_cond);

