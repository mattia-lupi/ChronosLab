int NLSP_PCG(const int np, const int prec_type, const int nn, const int nn_C,
             const int nn_K, const int ntv, const int *perm, const int *iperm,
             const double Tr_A, const int *iat_K, const int *ja_K, const double *coef_K,
             const int *it_U, const int *jcol_U, const double *coef_U,
             const double *D_inv, const int *iat_patt, const int *ja_patt,
             const int *iat_Tpatt, const int *ja_Tpatt,
             const int *pt_Z, const int *pt_col_Z, const double *mat_Z,
             const double *vec_f, const int itmax, int &iter, double &relres,
             double *vec_P);
