#pragma once

#include "precision.h"
#include <iostream>
#include <cstring>

void computeFirstQR(double *Ahat, ptrdiff_t sizeI, ptrdiff_t sizeJ, double *R, double *Rtriang, double *tau, double *work, ptrdiff_t lwork, ptrdiff_t &info);

void applyFirstQt(double *Ahat, ptrdiff_t sizeI, ptrdiff_t sizeJ, double *tau, double *a0k, double *work, ptrdiff_t lwork, ptrdiff_t &info);

void applyR(ptrdiff_t sizeJ, double *R, double *a0k, ptrdiff_t &info);

void applyQt(ptrdiff_t t, ptrdiff_t *sizeJ, ptrdiff_t *sizeI, ptrdiff_t *qStart, double *Ahat, 
             double *tau, double *a0k, ptrdiff_t nrowsA0k, ptrdiff_t ncolsA0k, double *work, ptrdiff_t lwork, ptrdiff_t &info);

void computeNewQR(ptrdiff_t t, ptrdiff_t *sizeI, ptrdiff_t *sizeJ, ptrdiff_t *qStart, double *Ahat, double *tau, double *R, 
                  double *Rtriang, double *work, ptrdiff_t lwork, ptrdiff_t &info);
