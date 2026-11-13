inline int get_diagpos(const int irow, const int *iat, const int *ja){
   int const istart = iat[irow];
   int const iend   = iat[irow + 1];
   for (int ind = istart; ind < iend; ind++) {
      if (ja[ind] == irow) return ind;
      if (ja[ind] > irow) break;
   }
   return -1; // diagonal not found
}
