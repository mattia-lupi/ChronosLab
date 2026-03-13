#include "wrCSRmat.h"

void wrCSRmat(FILE *ofile, const bool patt, const int nn_A, const int *const iat_A,
              const int *const ja_A, const double *const coef_A){

   if (patt){
      for (int i = 0; i < nn_A; i++){
         for (int j = iat_A[i]; j < iat_A[i+1]; j++){
            fprintf(ofile,"%10d %10d %1d\n",i+1,ja_A[j]+1,1);
         }
      }
   } else {
      for (int i = 0; i < nn_A; i++){
         for (int j = iat_A[i]; j < iat_A[i+1]; j++){
            fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_A[j]+1,coef_A[j]);
         }
      }
   }

}
