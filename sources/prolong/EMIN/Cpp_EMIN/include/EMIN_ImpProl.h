#define EMIN_INFO_SZ 12

int EMIN_ImpProl(const int np, const int itmax, const double en_tol, const double condmax,
                 const double maxwgt, const int prec_type, const int sol_type,
                 const int min_lfil, const int max_lfil, const int D_lfil, const int nn,
                 const int nn_C, const int ntv, const int nt_A, const int nt_P, const int nt_patt,
                 const int *fcnode, const int *iat_A, const int *ja_A, const double *coef_A,
                 const int *iat_Pin, const int *ja_Pin, const double *coef_Pin,
                 const int *iat_patt, const int *ja_patt, const double *const *TV,
                 int *&iat_Pout, int *&ja_Pout, double *&coef_Pout, double *info);
