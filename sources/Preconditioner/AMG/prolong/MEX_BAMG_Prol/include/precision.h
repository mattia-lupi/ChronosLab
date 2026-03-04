/**
 * @file precision.h
 * @brief This header is used to manage the computation precision.
 * @date November 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent
 */

#pragma once

#ifdef __cplusplus
#include <climits>
#else
#include <limits.h>
#endif

// Define types for iReg, iExt and iGlo
#define IREG_LONG 0
#define IEXT_LONG 0
#define IGLO_LONG 1

/**
 * @brief Regular integer.
 */
#if IREG_LONG
   typedef long int iReg;
#else
   typedef int iReg;
#endif

/**
 * @brief Extended integer.
 */
#if IEXT_LONG
   typedef long int iExt;
#else
   typedef int iExt;
#endif

/**
 * @brief Integer for global row/col indexing of a matrix
 */
#if IGLO_LONG
   typedef long int iGlo;
#else
   typedef int iGlo;
#endif

/**
 * @brief Regular real.
 */
typedef float rReg;

/**
 * @brief Extended real.
 */
typedef double rExt;

/**
 * @brief MPI regular integer.
 */
typedef int type_MPI_iReg;

/**
 * @brief MPI extended integer.
 */
typedef long int type_MPI_iExt;

/**
 * @brief OMP regular integer.
 */
typedef int type_OMP_iReg;

/**
 * @brief BLAS regular integer.
 */
typedef int type_BLAS_iReg;

/**
 * @brief BLAS extended integer.
 */
typedef long int type_BLAS_iExt;

/**
 * @brief LAPACK regular integer.
 */
typedef int type_LAPACK_iReg;

/**
 * @brief LAPACK extended integer.
 */
typedef long int type_LAPACK_iExt;

/**
 * @brief LAPACK extended real.
 */
typedef double type_LAPACK_rExt;

/**
 * @brief MKL regular integer.
 */
typedef int type_MKL_iReg;

/**
 * @brief Variable for the maximum MPI regular integer.
 */
#define MPI_MAX_iReg INT_MAX

/**
 * @brief Variable for MPI regular integer.
 */
#if IREG_LONG
   #define MPI_iReg MPI_LONG
#else
   #define MPI_iReg MPI_INT
#endif

/**
 * @brief Variable for MPI extended integer.
 */
#if IEXT_LONG
   #define MPI_iExt MPI_LONG
#else
   #define MPI_iExt MPI_INT
#endif

/**
 * @brief Variable for MPI global row/col indexing of a matrix.
 */
#if IGLO_LONG
   #define MPI_iGlo MPI_LONG
#else
   #define MPI_iGlo MPI_INT
#endif

/**
 * @brief Variable for MPI regular real.
 */
#define MPI_rReg MPI_FLOAT

/**
 * @brief Variable for MPI extended real.
 */
#define MPI_rExt MPI_DOUBLE
