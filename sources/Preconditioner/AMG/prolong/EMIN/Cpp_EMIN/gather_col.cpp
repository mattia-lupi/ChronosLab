void gather_col(const int nn, const int offset, const int *indices, const int len,
                const int *jcol, const double *coef, int &n_added, int *ja_out,
                double *coef_out){
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

   n_added = 0;
   int ii = 0;
   int jj = 0;
   while (jj < len){

      // Make sure that indices[ii] >= jcol[jj] (we are loading only upper part)
      while (indices[ii] < jcol[jj]){
         ii++;
         // Exit if there are no more indices
         if (ii == nn) return;
      }

      if (indices[ii] == jcol[jj]){
         // If jcol(jj) == indices(ii), load the term
         ja_out[n_added] = ii + offset;
         coef_out[n_added] = coef[jj];
         n_added++;
      }
      jj++;

   }

}
