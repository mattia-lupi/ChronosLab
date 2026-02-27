int gather_B_Z(const int np, const int nn, const int nn_C, const int ntv,
               const int *fcnode, const int *iat_patt, const int *ja_patt,
               const int *iat_Tpatt, const int *ja_Tpatt,
               const int *iat_A, const int *ja_A, const double *coef_A,
               const double *const *TV, int *&pt_Z, int *&pt_col_Z, double *&mat_Z,
               double *&vec_f, double *coef_P0);
