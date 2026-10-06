/**
 * @file mkiat_Tglo.h
 * @brief This header is used to manage an utility to transpose a CSRMAT.
 *        Computes the global iat_T and update WI pointers.
 * @date December 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

/**
 * @brief Computes the global iat_T and update WI pointers.
 * @param [in] current thread ID.
 * @param [in] nrows number of rows of the current thread.
 * @param [in] nthreads number of threads.
 * @param [in] firstrow first row of the current thread.
 * @param [inout] WI work array.
 * @param [in] iat_T iat array of the transposed matrix T.
 * @param [inout] nnz number of non-zeros for each thread.
 */
void mkiat_Tglo (const int myid, const int nrows, const int nthreads,
                 const int firstrow, int** RESTRICT WI,
                 const int* RESTRICT  nnz, int* RESTRICT iat_T);
