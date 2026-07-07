#include "precision.h"
#include <cstring>
#include <iostream>

void cpt_sam_adaptive_left(ptrdiff_t *iatk, ptrdiff_t *jak,double *coefk,
                           ptrdiff_t *iat0, ptrdiff_t *ja0,double *coef0,
                           ptrdiff_t nthread, ptrdiff_t n_step, ptrdiff_t step_size, double eps, ptrdiff_t nn_A,
                           ptrdiff_t *&iatN, ptrdiff_t *&jaN, double *&coefN, double &avg_resRelNorm);