/**
 * @file BAMG.h
 * @brief This header is used to manage BAMG computation.
 * @date March 2021
 * @version 1.1
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include <vector>
using namespace std;

#include "precision.h"    // to use: iReg, iExt, rExt
#include "BAMG_params.h"  // to use: BAMG_params

/**
 * @brief Computes BAMG prolongation.
 * @param [in] params
 * @param [in] nthreads number of threads used in set-up
 * @param [in] nn_S number of rows of the SoC matrix
 * @param [in] iat_S array of pointers of the SoC matrix
 * @param [in] ja_S array of column indices of the SoC matrix
 * @param [in] ntv number of test vectors
 * @param [in] fcnodes F/C indicator (<  0 ==> FINE node, >= 0 ==> COARSE node with number fcnodes)
 * @param [in] TV test vectors matrix of size [nn_S,ntv]
 * @param [out] nt_I number of entries of the interpolation matrix
 * @param [out] vec_iat_I array of pointers of the interpolation matrix
 * @param [out] vec_ja_I array of column indices of the interpolation matrix
 * @param [out] vec_coef_I array of coefficients of the interpolation matrix
 */
int BAMG ( const BAMG_params &params, const iReg nthreads, iReg nn_L, iReg nn_C,
           iReg nn_S, const iExt *const iat_S, const iReg *const ja_S, const iReg ntv,
           const iReg *const fcnodes, const rExt *const *const TV,
           iExt &nt_I, vector<iExt>& vec_iat_I, vector<iReg>& vec_ja_I,
           vector<rExt>& vec_coef_I, vector<iReg> &vec_c_mark );
