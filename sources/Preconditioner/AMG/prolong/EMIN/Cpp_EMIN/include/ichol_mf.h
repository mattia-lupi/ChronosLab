struct ichol_mf{
   int min_lfil;
   int max_lfil;
   int D_lfil;
   int max_nrows;
   int max_nnzr;
   int I_dim = 0;
   int R_dim = 0;
   int *I_scr = nullptr;
   double *R_scr = nullptr;
};
