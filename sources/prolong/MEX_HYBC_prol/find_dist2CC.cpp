#include "find_dist2CC.h"


void find_dist2CC(const int n_bad_FS, const int* const ind_bad_FS, const int* const ja_FS, const int* const iat_A,
                  const int* const ja_A, const double* const coef_A, const int* const coef_S, const int* const fcnodes, 
                  int &n_CC_2, int* const ja_CC_2, double* const coef_CC_2, int* const deg_CC_2){

   // Find distance two coarse node
   n_CC_2=0;
   for(int k=0; k<n_bad_FS; k++){
      int knod=ja_FS[ind_bad_FS[k]];
      for (int l=iat_A[knod]; l<iat_A[knod+1]; l++){
         if (coef_S[l] > 0){
            int lcol = ja_A[l];
            if (fcnodes[lcol] >= 0){
               if (deg_CC_2[lcol]==0){
                  ja_CC_2[n_CC_2] = lcol;
                  coef_CC_2[n_CC_2] = coef_A[l];
                  n_CC_2++;
               } 
               deg_CC_2[lcol]++;
            }
         }
      }
   }
}
