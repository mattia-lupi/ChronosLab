/**
 * @file add_new_neighs.h
 * @brief This header is used to manage the addition of a new belt of neighbour.
 * @date March 2021
 * @version 1.1
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include "precision.h"    // to use: iReg
#include "linsol_error.h" // to throw errors

/**
 * @brief Adds a new belt of neighbour to a given set of nodes performing a level
 *        set traversal.
 * @param [in] iat_S array of pointers of the SoC matrix.
 * @param [in] ja_S array of column indices of the SoC matrix.
 * @param [in] istart_neigh starting point of the last level in neigh.
 * @param [in] iend_neigh ending point of the last level in neigh.
 * @param [in] nmax maximum storage for neigh.
 * @param [inout] n_neigh number of neighbors in the list.
 * @param [inout] neigh list of neighbors.
 * @param [inout] WI indicator for the presence of a neighbor.
 */
void add_new_neighs(const iExt *const iat_S, const iReg *const ja_S,
                    const iReg istart_neigh, const iReg iend_neigh,
                    const iReg nmax, iReg &n_neigh, iReg *neigh, iReg *WI);
