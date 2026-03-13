int cpt_ichol_K(const int np,  const int min_lfil, const int max_lfil, const int D_lfil,
                const int n_blk, const int *pt_blk, const int nnmax_blk,
                const int ntmax_blk, const int nn_K, const int *iat_K, const int *ja_K,
                const double *coef_K, double &avg_lfil, int *&it_U, int *&jcol_U,
                double *&coef_U, double *&D_inv);
