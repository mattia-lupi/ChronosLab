#pragma once
#include <iostream>
#include <algorithm>

#include "precision.h"

int transpose(int nthread, const int nrows, const int *const iat, const int *const ja,     
              const double *const coef, int *iat_T, int *ja_T, double *coef_T);
