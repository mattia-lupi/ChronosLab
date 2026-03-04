#include "ri_sortsplit_nsy.h"

//----------------------------------------------------------------------------------------

// Select the ncut integers of I_vec corresponding to the absolute sum of R_vec and T_vec
// R_vec(i) >= R_vec(ncut) for i < ncut
// R_vec(i) <= R_vec(ncut) for i > ncut
// And perfroms the same permutations on integer vector I_vec
void ri_sortsplit_nsy(const int n, const int ncut, double* R_vec, int* I_vec){

   int first = 1;
   int last = n;

   // Check ncut consistency
   if (ncut < first || ncut > last) return;

   while(true){

      int mid = first;
      double val = R_vec[mid-1];

      for (int j = first+1; j < last+1; j ++){
         if (R_vec[j-1] > val){
            mid += 1;
            // Exchange
            double tmp  = R_vec[mid-1];
            int itmp = I_vec[mid-1];
            R_vec[mid-1] = R_vec[j-1];
            I_vec[mid-1] = I_vec[j-1];
            R_vec[j-1]   = tmp;
            I_vec[j-1]   = itmp;
         }
      }

      // Exchange
      double tmp      = R_vec[mid-1];
      R_vec[mid-1]    = R_vec[first-1];
      R_vec[first-1]  = tmp;

      int itmp      = I_vec[mid-1];
      I_vec[mid-1]   = I_vec[first-1];
      I_vec[first-1] = itmp;

      // Exit test
      if (mid == ncut) return;
      if (mid >  ncut){
         last = mid-1;
      }else{
         first = mid+1;
      }

   }

}
