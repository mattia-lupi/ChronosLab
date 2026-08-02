#pragma once
#include <iostream>
#include <algorithm>

#include "precision.h"

int transpose(const int nrows, const int *const iat, const int *const ja,     
              const double *const coef, double *coef_T);
