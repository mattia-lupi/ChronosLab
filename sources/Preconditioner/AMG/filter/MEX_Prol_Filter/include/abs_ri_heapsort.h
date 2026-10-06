#pragma once

#include <algorithm> // tu use: max

#include "swapi.h"
#include "swapr.h"

#if defined(_MSC_VER)
    #define RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#else
    #define RESTRICT
#endif

void abs_ri_heapsort(double* RESTRICT x1, int* RESTRICT x2, const int n);
