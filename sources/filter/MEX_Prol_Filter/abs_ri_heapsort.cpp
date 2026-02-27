
#include "abs_ri_heapsort.h"

//----------------------------------------------------------------------------------------

// Sorts an integer array x1 in such a way that x1(i) <= x1(i+1)
// and a real array in the same way
void abs_ri_heapsort(double* __restrict__ x1, int* __restrict__ x2, const int n){

   for (int node = 2; node < n+1; node ++){
      int i = node;
      int j = i/2;
      while( abs(x1[j-1]) < abs(x1[i-1]) ){
         swapr(x1[j-1],x1[i-1]);
         swapi(x2[j-1],x2[i-1]);
         i = j;
         j = i/2;
         if (i == 1) break;
      }
   }

   for (int i = n; i > 1; i --){
      swapr(x1[i-1],x1[0]);
      swapi(x2[i-1],x2[0]);
      int k = i - 1;
      int ik = 1;
      int jk = 2;
      if (k >= 3){
         if (abs(x1[2]) > abs(x1[1])) jk = 3;
      }
      bool cont_cycle = false;
      if (jk <= k){
         if (abs(x1[jk-1]) > abs(x1[ik-1])) cont_cycle = true;
      }
      while (cont_cycle){
         swapr(x1[jk-1],x1[ik-1]);
         swapi(x2[jk-1],x2[ik-1]);
         ik = jk;
         jk = ik*2;
         if (jk+1 <= k){
            if (abs(x1[jk]) > abs(x1[jk-1])) jk = jk+1;
         }
         cont_cycle = false;
         if (jk <= k){
            if (abs(x1[jk-1]) > abs(x1[ik-1])) cont_cycle = true;
         }
      }
   }

}

//----------------------------------------------------------------------------------------
