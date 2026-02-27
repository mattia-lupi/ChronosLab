#include "gather_fullsys.h"

void gather_fullsys(const int irow, const int mrow, const int *const vecinc,
                    const int nn, const int *const iat, const int *const ja,
                    const double *const coef, double *full_A, double *rhs_L,
                    double *rhs_U, bool &null_L, bool &null_U){

   // Loop over the rows(+1) of full_A
   null_L = true;
   null_U = true;
   int ind_row = 0;
   for (int i = 0; i < mrow; i++){
      // Load i-th row of full_A exploring the irow-th row of matrix A
      int ii = 0;
      int krow = vecinc[i];
      int jj = iat[krow];
      int endrow = iat[krow+1];
      while (ii < mrow){
         // Check ja[jj] >= vecinc[ii]
         while (ja[jj] < vecinc[ii]){
            jj++;
            if (jj == endrow){
               // End of row reached, set to zero the remaining entries and go to next row
               for (int k = ii; k < mrow; k++) full_A[ind_row+k] = 0.0;
               rhs_U[i] = 0.0;
               goto next_row;
            }
         }
         if (vecinc[ii] == ja[jj]){
            // Add this entry to full_A
            full_A[ind_row+ii] = coef[jj];
            ii++;
         } else {
            // This entry is null
            full_A[ind_row+ii] = 0.0;
            ii++;
         }
      } // End external while loop
      // Set the i-th entry of rhs_U
      ////////////////////////////////////////////
      //if (DEBUG){
      //   fprintf(dbfile,"jj %d ii %d ja[jj] %d vecinc[ii] %d\n",jj,ii,ja[jj],vecinc[ii]);
      //   fflush(dbfile);
      //}
      ////////////////////////////////////////////
      while(ja[jj] < irow){
         jj++;
         if (jj == endrow){
            rhs_U[i] = 0.0;
            goto next_row;
         }
      }
      if (irow == ja[jj]){
         // Add this entry to rhs_U
         rhs_U[i] = -coef[jj];
         null_U = false;
      } else {
         // This entry is null
         rhs_U[i] = 0.0;
      }
      next_row:;
      ind_row += mrow;
   } // End Row loop

   // Load rhs_L only exploring the vencinc[mrow]-th row of matrix A
   int ii = 0;
   int diag_row = irow;
   int jj = iat[diag_row];
   int endrow = iat[diag_row+1];
   while (ii < mrow){
      // Check ja[jj] >= vecinc[ii]
      while (ja[jj] < vecinc[ii]){
         jj++;
         if (jj == endrow){
            // End of row reached, set to zero the remaining entries and go to next row
            for (int k = ii; k < mrow; k++) rhs_L[k] = 0.0;
            goto end_rhs_L;
         }
      }
      if (vecinc[ii] == ja[jj]){
         // Add this entry to rhs_L
         rhs_L[ii] = -coef[jj];
         null_L = false;
         ii++;
      } else {
         // This entry is null
         rhs_L[ii] = 0.0;
         ii++;
      }
   } // End external while loop
   end_rhs_L:;

}
