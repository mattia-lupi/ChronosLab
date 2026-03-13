int Transp_Patt(const int nthreads, const int nrows, const int ncols, const int nterm,
                const int* __restrict__ iat, const int* __restrict__ ja,
                int* __restrict__ iat_T, int* __restrict__ ja_T, int* __restrict__ perm,
                int* __restrict__ iperm);
