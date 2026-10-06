#pragma once
#include "precision.h"
#include <cstring>
#include <iostream>

void cpt_sam_adaptive_left(iExt *iatk, iReg *jak,double *coefk, double *coefkT,
                           iExt *iat0, iReg *ja0,double *coef0,
                           iReg nthread, iReg n_step, iReg step_size, double eps, iExt nn_A,
                           iExt *&iatN, iReg *&jaN, double *&coefN, double &avg_resRelNorm);
