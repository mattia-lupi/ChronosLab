void spmatv(const int np, const int nn, const int* __restrict__ iat,
            const int* __restrict__ ja, const double* __restrict__ coef,
            const double* __restrict__ v_in, double* __restrict__ v_out){

   #pragma omp parallel for num_threads(np)
   for( int i = 0; i < nn; i++){
      int istart = iat[i];
      int iend = iat[i+1];
      v_out[i] = 0.;
      for( int j = istart; j < iend; j++) v_out[i] += coef[j]*v_in[ja[j]];
   }

}
