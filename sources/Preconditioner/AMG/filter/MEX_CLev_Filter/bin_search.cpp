#include "bin_search.h"

template<typename INT, typename TYPE>
INT bin_search(const TYPE ii, INT iend, const TYPE *const v){
   INT istart = 0;
   INT pos = iend/2;
   while (iend-istart > 1){
      if (v[pos] > ii){
         // Get the left interval
         iend = pos;
      } else {
         // Get the right interval
         istart = pos;
      }
      pos = (istart+iend)/2;
   }
   return pos;
}

template int bin_search<int,int>(const int, int, const int *const);
