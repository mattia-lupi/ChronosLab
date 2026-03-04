inline double abs_norm(const int nn, const double *v){

   double nrm = 0.0;
   for (int i = 0; i < nn; i++) nrm += abs(v[i]);
   return nrm;

}
