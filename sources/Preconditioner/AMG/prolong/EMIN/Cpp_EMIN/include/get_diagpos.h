inline int get_diagpos(const int irow, const int *iat, const int *ja){
   int ind = iat[irow];
   while (ja[ind] < irow) ind++;
   return ind;
}
