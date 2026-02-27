
#include "heapsort.h"

//----------------------------------------------------------------------------------------

// Sorts an array x1 in such a way that x1(i) <= x1(i+1) 
template <typename Tx, typename Tn>
void heapsort(Tx* __restrict__ x1, const Tn n){

   for (Tn node = 2; node < n+1; node ++){
      Tn i = node;
      Tn j = i/2;
      while( x1[j-1] < x1[i-1] ){
         SWAP(x1[j-1],x1[i-1]);
         i = j;
         j = i/2;
         if (i == 1) break;
      }
   }

   for (Tn i = n; i > 1; i --){
      SWAP(x1[i-1],x1[0]);
      Tn k = i - 1;
      Tn ik = 1;
      Tn jk = 2;
      if (k >= 3){
         if (x1[2] > x1[1]) jk = 3;
      }
      bool cont_cycle = false;
      if (jk <= k){
         if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
      }
      while (cont_cycle){
         SWAP(x1[jk-1],x1[ik-1]);
         ik = jk;
         jk = ik*2;
         if (jk+1 <= k){
            if (x1[jk] > x1[jk-1]) jk = jk+1;
         }
         cont_cycle = false;
         if (jk <= k){
            if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
         }
      }
   }

}

//----------------------------------------------------------------------------------------

// Instantiate template
template void heapsort<int,int>(int* __restrict__ x1, const int n);
