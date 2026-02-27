int gather_K(const int np, const int nn, const int nn_C,
             const int *iat_A, const int *ja_A, const double *coef_A,
             const int *iat_Pcol, const int *ja_Pcol, int &max_nrows_blk,
             int &max_nterm_blk, int &nn_K, int *&iat_K, int *&ja_K,
             double *&coef_K);
