void load_Jacobi(const int np, const int nn, const int nt, const int *iat_Tpatt,
                 const int *ja_Tpatt, const int *iat, const int *ja, const double *coef,
                 double *scr, double *D_inv){

   // Load diagonal entries in scr
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn; i++){
      int ind = iat[i];
      while (ja[ind] < i) ind++;
      scr[i] = 1.0 / coef[ind];
   }

   // Load diagonal entries in scr
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nt; i++){
      D_inv[i] = scr[ja_Tpatt[i]];
   }

}
