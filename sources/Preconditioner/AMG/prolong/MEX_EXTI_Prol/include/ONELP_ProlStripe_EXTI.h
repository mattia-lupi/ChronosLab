#include <vector>

int ONELP_ProlStripe_EXTI(const int firstrow, const int lastrow, const int nn_loc,
                          const int nn_A, const int nt_A,
                          const int *const iat_A, const int *const ja_A,
                          const double *const coef_A, const int *const coef_S,
                          const int *const fcnodes,
                          int &nn_P, int &nt_P, int *iat_P,
                          std::vector<int> &ja_P, std::vector<double> &coef_P);
