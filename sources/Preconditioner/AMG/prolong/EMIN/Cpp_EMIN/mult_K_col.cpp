double mult_K_col(const int nn, const int *indices, const int len, const int *jcol,
                  const double *coef, const double *vec_x){
   ///////////////////////////////////////////////////////
   /*
   cout << "# indices " << nn << endl;
   cout << "# jcol " << len << endl;
   cout << " INDICES ";
   for (int i = 0; i < nn; i++) cout << " " << indices[i];
   cout << endl << endl;
   cout << " JCOL ";
   for (int i = 0; i < len; i++) cout << " " << jcol[i];
   cout << endl << endl;
   cout << " COEF ";
   for (int i = 0; i < len; i++) cout << " " << coef[i];
   cout << endl << endl;
   */
   ///////////////////////////////////////////////////////

   double product = 0.0;
   int ii = 0;
   int jj = 0;
   while (jj < len){

      // Make sure that indices[ii] >= jcol[jj]
      while (indices[ii] < jcol[jj]){
         ii++;
         // Exit if there are no more indices
         if (ii == nn) return product;
      }

      if (indices[ii] == jcol[jj]){
         // If jcol(jj) == indices(ii), load the term
         product += coef[jj]*vec_x[ii];
      }
      jj++;

   }

   return product;

}
