#pragma once

#include <iostream>
#include <cstring>
#include "precision.h"

void fullA0k(ptrdiff_t nn_A, ptrdiff_t *iat0, ptrdiff_t *ja0, double *coef0, ptrdiff_t k, double *A0k);

void findNonZeroInColJ(ptrdiff_t *J, ptrdiff_t *iatk, ptrdiff_t *jak, ptrdiff_t n2, ptrdiff_t *I, ptrdiff_t &sizeI);

void getA0k(double *a0k, ptrdiff_t *I, ptrdiff_t sizeI, ptrdiff_t oldSizeI, ptrdiff_t *iat0, ptrdiff_t *ja0, double *coef0, ptrdiff_t k);

void getAhat( ptrdiff_t *I, ptrdiff_t sizeI, ptrdiff_t *J, ptrdiff_t Jstart, ptrdiff_t Jend,
             ptrdiff_t *iatk, ptrdiff_t *jak, double *coefk, double *Ahat, ptrdiff_t &Astart);

void getAJ(ptrdiff_t *J, ptrdiff_t Jsize, ptrdiff_t nn_A, ptrdiff_t *iatk, ptrdiff_t *jak, double *coefk, double *AJ);

void fillL(ptrdiff_t *L, double *res, ptrdiff_t nn_A, ptrdiff_t &usedL);

void findJtilde(ptrdiff_t *Jtilde, ptrdiff_t &JtildeSize, ptrdiff_t *L, ptrdiff_t sizeL, ptrdiff_t *iatk, ptrdiff_t *jak, ptrdiff_t *J, ptrdiff_t sizeJ);

void fullAJtilde(ptrdiff_t nn_A, ptrdiff_t *iatk, ptrdiff_t *jak, double *coefk, ptrdiff_t *Jtilde, ptrdiff_t JtildeSize, double *AJtilde);
