void gather_f(const int nn, const int *indices, const int len, const int *jcol,
              const double *coef, double *vec_f){

   int ii = 0;
   int jj = 0;
   if (len > 0){
      while (ii < nn){
         // Advance jj until jcol[jj] < indices[ii]
         while (jcol[jj] < indices[ii]){
            jj++;
            if (jj == len) goto exit_loop;
         }
         if (jcol[jj] == indices[ii]){
            vec_f[ii] = -coef[jj];
         } else {
            vec_f[ii] = 0.0;
         }
         ii++;
      }
   }
   exit_loop: ;
   // Nullify remaining entries
   for (int k = ii; k < nn; k++) vec_f[k] = 0.0;

/******************************************************************************
 * CICLO INVERTITO
   while (jj < len){

      // Make sure that indices[ii] >= jcol[jj]
      while (indices[ii] < jcol[jj]){
         vec_f[n_added] = 0.0;
         n_added++;
         ii++;
         // Exit if there are no more indices
         if (ii == nn) return;
      }

      if (indices[ii] == jcol[jj]){
         // If jcol(jj) == indices(ii), load the term
         vec_f[n_added] = -coef[jj];
         n_added++;
      }
      jj++;

   }
   // Complete vec_f with zeros
   for (int j = n_added; j < nn; j++) vec_f[j] = 0.0;
***********************************************************************************/

}
