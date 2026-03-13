void merge_row_patt(const int len_1, const int *const ja_1, const double *const coef_1,
                    const int len_2, const int *const ja_2,  const double *const coef_2,
                    int &len_out, int *ja_out, double *coef_out){

   int i = 0;
   int j = 0;
   int k = 0;
   while ( (i < len_1) && (j < len_2)){
      if (ja_1[i] < ja_2[j]){
         ja_out[k]   = ja_1[i];
         coef_out[k] = coef_1[i];
         k++;
         i++;
      } else if (ja_1[i] == ja_2[j]){
         ja_out[k]   = ja_1[i];
         coef_out[k] = coef_1[i];
         k++;
         i++;
         j++;
      } else {
         ja_out[k]   = ja_2[j];
         coef_out[k] = 0.0;
         k++;
         j++;
      }
   }
   while (i < len_1){
      ja_out[k] = ja_1[i];
      coef_out[k] = coef_1[i];
      k++;
      i++;
   }
   while(j < len_2){
      ja_out[k] = ja_2[j];
      coef_out[k] = 0.0;
      k++;
      j++;
   }
   len_out = k;

}
