#include <stdlib.h>
#include <omp.h>
#include "FilterComp_stripe.h"

int FilterComp(const int nthreads, const double tau, const int nn_A,
               int *iat_A, int *ja_A, double *coef_A, const int nt_patt,
               const int *iat_patt, const int *ja_patt, const int ntv,
               const double *const *TV, int &nt_AC, int *&iat_AC,
               int *&ja_AC, double *&coef_AC){

   // Init error code
   int ierr = 0;

   // Allocate scratch
   int *ridv_i = (int*) malloc((nthreads+1)*sizeof(int));
   if (ridv_i == nullptr) return 1;

   // Distribute rows among processes
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   # pragma omp parallel num_threads(nthreads)
   {
      int istart_AC, shift, istart_scr, iend_scr;

      // Init local error code
      int ierr_L = 0;

      // Retrieve processor ID
      int myid = omp_get_thread_num();

      // Partition rows between threads
      int resto = nn_A%nthreads;
      int nrows = nn_A/nthreads;
      int firstrow;
      if (myid < resto){
         nrows++;
         firstrow = myid*nrows;
      } else {
         firstrow = myid*nrows + resto;
      }

      // Allocate scratches for the filtered matrix
      int ntmax = iat_A[firstrow+nrows] - iat_A[firstrow] + nn_A;
      int     *iat_scr = (int*) malloc((nrows+1)*sizeof(int));
      int      *ja_scr = (int*) malloc((ntmax)*sizeof(int));
      double *coef_scr = (double*) malloc((ntmax)*sizeof(double));
      if (iat_scr == nullptr || ja_scr == nullptr || coef_scr == nullptr){
         ierr_L = 1; goto exit_pragma;
      }

      int nt_AC_loc;
      ierr_L = FilterComp_stripe(tau,firstrow,nrows,nn_A,&(iat_A[firstrow]),ja_A,coef_A,
                                 nt_patt,&(iat_patt[firstrow]),ja_patt,ntv,TV,
                                 nt_AC_loc,iat_scr,ja_scr,coef_scr);
      if (ierr_L != 0) {ierr_L = 2; goto exit_pragma;}

      // Store local number of entries
      ridv_i[myid+1] = nt_AC_loc;
      #pragma omp barrier

      // Allocate the prolongation matrix
      #pragma omp single
      {
         ridv_i[0] = 0;
         for (int ip = 0; ip < nthreads; ip++){
             ridv_i[ip+1] = ridv_i[ip] + ridv_i[ip+1];
         }
         nt_AC = ridv_i[nthreads];
         iat_AC = (int*) malloc((nn_A+1)*sizeof(int));
         ja_AC = (int*) malloc((nt_AC)*sizeof(int));
         coef_AC = (double*) malloc((nt_AC)*sizeof(double));
         if (iat_AC == nullptr || ja_AC == nullptr || coef_AC == nullptr) ierr = 1;
      }
      if (ierr != 0) goto exit_pragma;

      // Store iat_scr, ja_scr and coef_scr in the filtered matrix
      istart_AC = ridv_i[myid];
      shift = firstrow;
      istart_scr;
      iend_scr = iat_scr[0];
      for (int irow = 0; irow < nrows; irow++){
         iat_AC[shift+irow] = istart_AC;
         istart_scr = iend_scr;
         iend_scr = iat_scr[irow+1];
         int len_AC = iend_scr - istart_scr;
         for (int k = 0; k < len_AC; k++){
            ja_AC[istart_AC+k] = ja_scr[istart_scr+k];
            coef_AC[istart_AC+k] = coef_scr[istart_scr+k];
         }
         istart_AC += len_AC;
      }
      if (myid == nthreads-1) iat_AC[nn_A] = ridv_i[nthreads];

      // Free scratches
      free(iat_scr);
      free(ja_scr);
      free(coef_scr);

      // Exit point
      exit_pragma: ;

      // Check error codes
      #pragma omp atomic
      ierr += ierr_L;

   }
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Free scratches
   free(ridv_i);

   return ierr;

}
