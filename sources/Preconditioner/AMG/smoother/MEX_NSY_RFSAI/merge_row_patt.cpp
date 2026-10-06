#include "merge_row_patt.h"

int count_merged_row_patt(const int len_1, const int *const ja_1,
                          const int len_2, const int *const ja_2){
   int ind_1 = 0;
   int ind_2 = 0;
   int len = 0;

   while (ind_1 < len_1 && ind_2 < len_2){
      if (ja_1[ind_1] < ja_2[ind_2]){
         len++;
         ind_1++;
      } else if (ja_2[ind_2] < ja_1[ind_1]){
         len++;
         ind_2++;
      } else {
         len++;
         ind_1++;
         ind_2++;
      }
   }
   while (ind_1 < len_1){
      len++;
      ind_1++;
   }
   while (ind_2 < len_2){
      len++;
      ind_2++;
   }
   return len;
}

void merge_row_patt(const int len_1, const int *const ja_1, const double *const coef_1,
                    const int len_2, const int *const ja_2,  const double *const coef_2,
                    int &len_3, int *const ja_3, double *const coef_3){
   (void)coef_2;
   int ind_1 = 0;
   int ind_2 = 0;
   int ind_3 = 0;

   while (ind_1 < len_1 && ind_2 < len_2){
      if (ja_1[ind_1] < ja_2[ind_2]){
         ja_3[ind_3] = ja_1[ind_1];
         coef_3[ind_3] = coef_1[ind_1];
         ind_1++;
         ind_3++;
      } else if (ja_2[ind_2] < ja_1[ind_1]){
         ja_3[ind_3] = ja_2[ind_2];
         coef_3[ind_3] = 0.0;
         ind_2++;
         ind_3++;
      } else {
         ja_3[ind_3] = ja_1[ind_1];
         coef_3[ind_3] = coef_1[ind_1];
         ind_1++;
         ind_2++;
         ind_3++;
      }
   }

   while (ind_1 < len_1){
      ja_3[ind_3] = ja_1[ind_1];
      coef_3[ind_3] = coef_1[ind_1];
      ind_1++;
      ind_3++;
   }

   while (ind_2 < len_2){
      ja_3[ind_3] = ja_2[ind_2];
      coef_3[ind_3] = 0.0;
      ind_2++;
      ind_3++;
   }

   len_3 = ind_3;
}
