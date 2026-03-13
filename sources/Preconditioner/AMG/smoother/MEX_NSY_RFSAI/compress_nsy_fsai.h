#include <iostream>  // to use: cout,endl
#include <algorithm> // to use: fill_n
using namespace std;

int compress_nsy_fsai(const int nn, int &nt_FL, int &nt_FU,
                      int *&iat_FL, int *&ja_FL, int *&iat_FU, int *&ja_FU,
                      double *&coef_FL, double *&coef_FU);
