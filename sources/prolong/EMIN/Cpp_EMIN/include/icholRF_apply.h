void icholRF_apply(const int np, const int nblk, const int* __restrict__ pt_blk,
                   const int* __restrict__ iU, const int* __restrict__ jU,
                   const double* __restrict__ coef_U, const double* __restrict__ D_inv,
                   const double* __restrict__ vec, double* __restrict__ pvec);
