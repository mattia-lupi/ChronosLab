#include <stdlib.h>
#include <omp.h>

int Prol_add_Cnodes(const int np, const int nn, const int nn_C,
                    const int* __restrict__ fcnode, const int* __restrict__ iat_in,
                    const int* __restrict__ ja_in, const double* __restrict__ coef_in,
                    int *&iat_out, int *&ja_out, double *&coef_out){

   // Allocate new prolongation
   int nt_out = iat_in[nn] + nn_C;
   iat_out   = (int*) malloc( (nn+1)*sizeof(int) );
   ja_out    = (int*) malloc( nt_out*sizeof(int) );
   coef_out  = (double*) malloc( nt_out*sizeof(double) );
   if (iat_out == nullptr || ja_out == nullptr || coef_out == nullptr ) return 1;

   // Allocate scratch
   int *pt_loc = (int*) malloc( (np+1)*sizeof(int) );
   if (pt_loc == nullptr) return 2;

   #pragma omp parallel num_threads(np)
   {
       // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstrow, nrowth, lastrow;
      if (mythid <= resto) {
         nrowth = bsize+1;
         firstrow = mythid*nrowth;
         if (mythid == resto) nrowth--;
      } else {
         nrowth = bsize;
         firstrow = mythid*bsize + resto;
      }
      lastrow = firstrow + nrowth;

      // Count the number of coarse nodes in each partition
      int count = 0;
      for (int i = firstrow; i < lastrow; i++) if (fcnode[i] >= 0) count++;
      pt_loc[mythid+1] = count;
      #pragma omp barrier

      // Reduce pt_loc
      #pragma omp single
      {
         pt_loc[0] = 0;
         for (int i = 0; i < np; i++ ){
            pt_loc[i+1] += pt_loc[i];
         }
      }

      // Complete each chuck of prolongation
      int iend = iat_in[firstrow];
      int ind_out = iend + pt_loc[mythid];
      for (int irow = firstrow; irow < lastrow; irow++){
         iat_out[irow] = ind_out;
         int istart = iend;
         iend = iat_in[irow+1];
         int node = fcnode[irow];
         if (node >= 0){
            // This is a COARSE node
            ja_out[ind_out] = node;
            coef_out[ind_out] = 1.0;
            ind_out++;
         } else {
            // This is a FINE node
            for (int j = istart; j < iend; j++){
               ja_out[ind_out] = ja_in[j];
               coef_out[ind_out] = coef_in[j];
               ind_out++;
            }
         }
      }
      if (mythid == np-1) iat_out[nn] = nt_out;

   }

   // Free scratch
   free(pt_loc);

   return 0;

}
