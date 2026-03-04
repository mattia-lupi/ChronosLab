void spmatv_blk(const int np, const int nblk, const int* __restrict__ pt_blk,
                const int* __restrict__ iat, const int* __restrict__ ja,
                const double* __restrict__ coef, const double* __restrict__ v_in,
                double* __restrict__ v_out);
