/**
 * @file inl_blas1.h
 * @brief This header is used to manage the LEVEL 1 BLAS.
 * @date September 2019
 * @version 1.0
 * @author Carlo Janna
 * @author Giovanni Isotton
 * @par License
 *      This program is intended for private use only and can not be distributed
 *      elsewhere without authors' consent.
 */

#pragma once

#include <cmath> // to use: abs
//#include <math.h> // to use: sqrt

/**
 * @brief DCOPY
 */
inline void inl_dcopy(const int n, const double *const v1, const int k1,
                      double *const v2, const int k2){
   int j = 0;
   for (int i = 0; i < n; i += k1){
      v2[i] = v1[j];
      j += k2;
   }
   return;
}

/**
 * @brief DDOT
 */
inline double inl_ddot(const int n, const double *const v1, const int k1,
                     const double *const v2, const int k2){
   double ddot = 0.0;
   int j = 0;
   for (int i = 0; i < n; i += k1){
      ddot += v1[i]*v2[j];
      j += k2;
   }
   return ddot;
}

/**
 * @brief DNRM1
 */
inline double inl_dnrm1(const int n, const double *const v1, const int k1){
   double dnrm1 = 0.0;
   for (int i = 0; i < n; i += k1) dnrm1 += std::abs(v1[i]);
   return dnrm1;
}

/**
 * @brief DNRM2
 */
inline double inl_dnrm2(const int n, const double *const v1, const int k1){
   double dnrm2 = 0.0;
   for (int i = 0; i < n; i += k1) dnrm2 += (v1[i])*(v1[i]);
   return sqrt(dnrm2);
}

/**
 * @brief DAXPY
 */
inline void inl_daxpy(const int n, const double alpha, const double *const x, 
                      const int k1, double *const y, const int k2){
   int j = 0;
   for (int i = 0; i < n; i += k1){
      y[j] = y[j] + alpha*x[i];
      j += k2;
   }
   return;
}


