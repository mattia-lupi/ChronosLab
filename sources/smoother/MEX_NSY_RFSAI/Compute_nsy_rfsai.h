#pragma once

#include <iostream>  // to use: cout,endl
#include <cstring>   // to use: memcpy
//#include <fstream>   // to use: ifstream,ofstream
//#include <math.h>
using namespace std;

#include "extract_diag.h"
#include "cpt_nsy_rfsai.h"
#include "compress_nsy_fsai.h"
#include "SymmetrizePattern.h"
#include "wrCSRmat.h"

int Compute_nsy_rfsai(const int nstep, const int step_size, const double eps,
                      const int nn_A, const int *iat_A, const int *ja_A,
                      const double *coef_A, int *&iat_FL, int *&ja_FL,
                      double *&coef_FL, int *&iat_FU, int *&ja_FU, double *&coef_FU);
