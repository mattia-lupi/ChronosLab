///////////////////////////////
#include <iostream>
using namespace std;
/////////////////////////////////

#include <stdlib.h>
#include <omp.h>
#include "FilterProl.h"
#include "FilterProl_stripe.h"

int FilterProl(const int nthreads, const double perc, const double tol,
               const int nn_P, int *iat_P, int *ja_P, double *coef_P, 
               const int ntv, const double *const *TV, int &nt_PF, int *&iat_PF,
               int *&ja_PF, double *&coef_PF){

   // Init error code
   int ierr = 0;

   // Allocate scratch
   int *ridv_i = (int*) malloc((nthreads+1)*sizeof(int));
   if (ridv_i == nullptr) return 1;

   // Distribute rows among processes
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   # pragma omp parallel num_threads(nthreads)
   {
      int istart_PF, shift, istart_scr, iend_scr;

      // Init local error code
      int ierr_L = 0;

      // Retrieve processor ID
      int myid = omp_get_thread_num();

     // Partition rows between threads
      int resto = nn_P%nthreads;
      int nrows = nn_P/nthreads;
      int firstrow;
      if (myid < resto){
         nrows++;
         firstrow = myid*nrows;
      } else {
         firstrow = myid*nrows + resto;
      }

      // Allocate scratches for the filtered matrix
      int ntmax = iat_P[firstrow+nrows] - iat_P[firstrow] + nn_P;
      int     *iat_scr = (int*) malloc((nrows+1)*sizeof(int));
      int      *ja_scr = (int*) malloc((ntmax)*sizeof(int));
      double *coef_scr = (double*) malloc((ntmax)*sizeof(double));
      if (iat_scr == nullptr || ja_scr == nullptr || coef_scr == nullptr){
         ierr_L = 1; goto exit_pragma;
      }

      int nt_PF_loc;
      ierr_L = FilterProl_stripe(perc,tol,firstrow,nrows,nn_P,&(iat_P[firstrow]),
                                 ja_P,coef_P,ntv,TV,nt_PF_loc,iat_scr,ja_scr,coef_scr);
      if (ierr_L != 0) {ierr_L = 2; goto exit_pragma;}

      // Store local number of entries
      ridv_i[myid+1] = nt_PF_loc;
      #pragma omp barrier

      // Allocate the prolongation matrix
      #pragma omp single
      {
         ridv_i[0] = 0;
         for (int ip = 0; ip < nthreads; ip++){
             ridv_i[ip+1] = ridv_i[ip] + ridv_i[ip+1];
         }
         nt_PF = ridv_i[nthreads];
         iat_PF = (int*) malloc((nn_P+1)*sizeof(int));
         ja_PF = (int*) malloc((nt_PF)*sizeof(int));
         coef_PF = (double*) malloc((nt_PF)*sizeof(double));
         if (iat_PF == nullptr || ja_PF == nullptr || coef_PF == nullptr) ierr = 1;
      }
      if (ierr != 0) goto exit_pragma;

      // Store iat_scr, ja_scr and coef_scr in the filtered matrix
      istart_PF = ridv_i[myid];
      shift = firstrow;
      istart_scr;
      iend_scr = iat_scr[0];
      for (int irow = 0; irow < nrows; irow++){
         iat_PF[shift+irow] = istart_PF;
         istart_scr = iend_scr;
         iend_scr = iat_scr[irow+1];
         int len_PF = iend_scr - istart_scr;
         for (int k = 0; k < len_PF; k++){
            ja_PF[istart_PF+k] = ja_scr[istart_scr+k];
            coef_PF[istart_PF+k] = coef_scr[istart_scr+k];
         }
         istart_PF += len_PF;
      }
      if (myid == nthreads-1) iat_PF[nn_P] = ridv_i[nthreads];

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
