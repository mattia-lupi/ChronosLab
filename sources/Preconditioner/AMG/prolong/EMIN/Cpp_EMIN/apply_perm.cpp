// Applies permutation perm to the input vector vec_in producing vec_out

void apply_perm(const int np, const int nn, const int *perm, const double *vec_in,
                double *vec_out){

   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn; i++) vec_out[i] = vec_in[perm[i]];

}
