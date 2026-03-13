/**
 * @file ProlStripe_BAMG.h
 * @brief This header is used to manage BAMG computation among threads.
 * @date March 2021
 * @version 1.1
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

/**
 * @brief Inner part of BAMG prolongation that is used to distributed work among threads.
 */
void ProlStripe_BAMG(const BAMG_params& params, iReg firstrow_0, iReg firstrow, iReg lastrow,
                     iReg nn_S, iReg ntvecs, const iExt *const iat_S, const iReg *const ja_S,
                     const iReg *const fcnodes,
                     const rExt *const *const TV, iExt &nt_P, iExt *iat_P, iReg *ja_P,
                     rExt *coef_P, iReg *c_mark, iReg *dist_count);
