// Transpose pattern indices.
void mvjcol( const int firstrow, const int nrows, const int nequ, const int nterm,
             const int nterm_T, const int* __restrict__ iat, const int* __restrict__ ja,
             int* __restrict__ ja_T, int* __restrict__ punt, int* __restrict perm){

   // Transpose local stripe
   int shift = firstrow;
   for ( int i = 0; i < nrows; i++ ) {
      for ( int j = iat[i]; j < iat[i+1]; j++ ) {
         int irow = ja[j];
         perm[j] = punt[irow];
         ja_T[punt[irow]] = i + shift;
         punt[irow]++;
      }
   }

}
