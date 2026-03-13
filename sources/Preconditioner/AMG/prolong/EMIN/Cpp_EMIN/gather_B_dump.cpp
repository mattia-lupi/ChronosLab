#include <stdlib.h>
#include <omp.h>
#include <limits>
using namespace std;

#include "DebEnv.h"
#include "parm_EMIN.h"

int gather_B_dump(const int np, const int nn, const int nn_C, const int ntv,
                  const int *fcnode, const int *iat_patt, const int *ja_patt,
                  const double *const *TV, double *&mat_B){

   // Init error code
   int ierr = 0;

   // Allocate B
   int nrows_B = iat_patt[nn];
   int nterm_B = ntv*nrows_B;
   mat_B = (double*) malloc( nterm_B*sizeof(double) );
   if (mat_B == nullptr) return ierr = 3;

   // Allocate shared scratches
   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   if (c2glo == nullptr)  return ierr = 1;

   #pragma omp parallel num_threads(np)
   {
      // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstcol, ncolth, lastcol;
      if (mythid <= resto) {
         ncolth = bsize+1;
         firstcol = mythid*ncolth;
         if (mythid == resto) ncolth--;
      } else {
         ncolth = bsize;
         firstcol = mythid*bsize + resto;
      }
      lastcol = firstcol + ncolth;

      // Estimate number of rows and first position for this chunk of columns
      int pos_g = iat_patt[firstcol];
      int pos_B = pos_g*ntv;

      // Set proper position in B
      double *BB_scr = &(mat_B[pos_B]);

      // Create mapping from coarse node numbering to global (original) numbering
      #pragma omp for
      for (int i = 0; i < nn; i++){
         int k = fcnode[i];
         if (k >= 0) c2glo[k] = i;
      }

      // Loop over the current chunk of columns
      int ind_BB, nnz_BB;
      ind_BB = 0;
      int istart_patt, iend_patt;
      iend_patt = iat_patt[firstcol];
      for (int icol = firstcol; icol < lastcol; icol++){
         // Check that this is a FINE node
         if (fcnode[icol] < 0){

            istart_patt = iend_patt;
            iend_patt = iat_patt[icol+1];
            int nr_BB_loc = iend_patt-istart_patt;

            // Check that the row is not empty
            if (nr_BB_loc > 0){
               // Copy TV into BB_scr
               int k = ind_BB;
               for (int i = istart_patt; i < iend_patt; i++){
                  int i_F = c2glo[ja_patt[i]];
                  int kk = k;
                  for (int j = 0; j < ntv; j++){
                     BB_scr[kk] = TV[i_F][j];
                     kk += nr_BB_loc;
                  }
                  k++;
               }

               // Compute the number of entries of Q
               nnz_BB = nr_BB_loc*ntv;

               // Update pointers in BB, g and tau
               ind_BB += nnz_BB;

            }

         }

      } // End loop over columns

   } // End of parallel region

   // Free shared scratches
   free(c2glo);

   return ierr;

}
