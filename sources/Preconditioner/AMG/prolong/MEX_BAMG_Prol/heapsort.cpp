
#include "heapsort.h"

//----------------------------------------------------------------------------------------

// Sorts an array x1 in such a way that x1(i) <= x1(i+1)
template <typename Tx, typename Tn>
void heapsort(Tx* __restrict__ x1, const Tn n){

   for (Tn node = 1; node < n; node ++){
      Tn i = node;
      Tn j = i/2;
      while( x1[j] < x1[i] ){
         SWAP(x1[j],x1[i]);
         i = j;
         j = i/2;
         if (i == 0) break;
      }
   }

   for (Tn i = n-1; i > 0; i --){
      SWAP(x1[i],x1[0]);
      Tn k = i-1;
      Tn ik = 0;
      Tn jk = 1;
      if (k >= 2){
         if (x1[2] > x1[1]) jk = 2;
      }
      bool cont_cycle = false;
      if (jk <= k){
         if (x1[jk] > x1[ik]) cont_cycle = true;
      }
      while (cont_cycle){
         SWAP(x1[jk],x1[ik]);
         ik = jk;
         jk = ik*2;
         if (jk+1 <= k){
            if (x1[jk+1] > x1[jk]) jk = jk+1;
         }
         cont_cycle = false;
         if (jk <= k){
            if (x1[jk] > x1[ik]) cont_cycle = true;
         }
      }
   }

}

//----------------------------------------------------------------------------------------

// Instantiate template
#if !IREG_LONG==IEXT_LONG
   template void heapsort<iReg,iReg>(iReg* __restrict__ x1, const iReg n);
   template void heapsort<iExt,iReg>(iExt* __restrict__ x1, const iReg n);
   template void heapsort<rExt,iReg>(rExt* __restrict__ x1, const iReg n);
   template void heapsort<iReg,iExt>(iReg* __restrict__ x1, const iExt n);
#endif
#if IREG_LONG & IEXT_LONG & IGLO_LONG
   template void heapsort<type_MPI_iReg,type_MPI_iReg>(type_MPI_iReg* __restrict__ x1,
                                                       const type_MPI_iReg n);
#endif
template void heapsort<iExt,iExt>(iExt* __restrict__ x1, const iExt n);
template void heapsort<rExt,iExt>(rExt* __restrict__ x1, const iExt n);
