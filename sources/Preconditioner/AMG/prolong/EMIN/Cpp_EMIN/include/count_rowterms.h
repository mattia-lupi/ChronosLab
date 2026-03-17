/**
 * @file count_rowterms.h
 * @brief This header is used to manage an utility to transpose a CSRMAT.
 *        It counts how many terms assigned to a thread in each row.
 * @date December 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

/**
 * @brief Counts the number of non-zeroes assigned to a thread in each row.
 * @param [in] nequ number of equations.
 * @param [in] nterm number of non-zeros of the thread.
 * @param [in] ja column indexes.
 * @param [out] WI number of non-zeros for each row.
 */
void count_rowterms( const int nequ, const int nterm, const int* __restrict__ ja,
                     int* __restrict__ WI );
