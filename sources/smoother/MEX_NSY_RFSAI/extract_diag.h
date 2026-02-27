#include <iostream>
#include <stdlib.h>
using namespace std;

void extract_diag(const int nrows, const int *const iat, const int *const ja,
                  const double *const coef, double *&diag);
