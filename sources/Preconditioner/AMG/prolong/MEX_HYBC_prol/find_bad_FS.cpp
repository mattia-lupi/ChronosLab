#include "find_bad_FS.h"

void find_bad_FS( const int n_FS, const int* const ja_FS, const int* const WI, 
                  const int * const iat_A, const int * const ja_A, int &n_bad_FS, int * const ind_bad_FS){
   n_bad_FS=0;
   for (int k=0; k<n_FS; k++){
      int knod=ja_FS[k];
      bool flag_bad_FS=true;
      for (int l=iat_A[knod]; l<iat_A[knod+1]; l++){
         int lcol=ja_A[l];
         if ( WI[lcol] > 0){
            flag_bad_FS=false;
            break;
         } 
      }
      if (flag_bad_FS){
         ind_bad_FS[n_bad_FS]=k;
         n_bad_FS++;
      }
   }
}
