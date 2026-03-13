/**
 * @file mkiat_Tloc.h
 * @brief This header is used to manage an utility to transpose a CSRMAT.
 *        Each thread computes the pointers to the beginning of each row of its part of the matrix
 * @date December 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

using namespace std;

/**
 * @brief Each thread computes the pointers to the beginning of each row of its part of the matrix
 * @param [in] nrows number of rows of the current thread.
 * @param [in] nequ number of equations.
 * @param [in] nthreads number of threads.
 * @param [in] firstrow first row of the current thread.
 * @param [in] WI work array WI[j+1,i] number of non-zeros going in T in i-th row and handled by j-th thread.
 * @param [out] iat_T iat array of the transposed matrix T.
 * @param [out] nnz number of non-zeros.
 */
void mkiat_Tloc (const int nrows, const int nequ, const int nthreads, const int firstrow,
                 int** __restrict__ WI, int* __restrict__ iat_T, int& nnz  );
