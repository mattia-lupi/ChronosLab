
#include "add_new_neighs.h"

//----------------------------------------------------------------------------------------

// Adds a new belt of neighbour to a given set of nodes performing a level set traversal.
void add_new_neighs(const iExt *const iat_S, const iReg *const ja_S,
                    const iReg istart_neigh, const iReg iend_neigh,
                    const iReg nmax, iReg &n_neigh, iReg *neigh, iReg *WI){

   // Loop over the last added level
   for (iReg i_neigh = istart_neigh; i_neigh < iend_neigh; i_neigh++){ // Last_level_loop

      iReg jnod   = neigh[i_neigh];
      iExt istart = iat_S[jnod];
      iExt iend   = iat_S[jnod+1];
      // Check storage
      if ( n_neigh + istart-iend > nmax) {
         throw linsol_error ("add_new_neighs","storage problem");
      }
      for (iExt i = istart; i < iend; i++){
         iReg jcol = ja_S[i];
         if (WI[jcol] == 0){
            neigh[n_neigh] = jcol;
            n_neigh++;
            WI[jcol] = 1;
         }
      }

   } // Last_level_loop end

}
