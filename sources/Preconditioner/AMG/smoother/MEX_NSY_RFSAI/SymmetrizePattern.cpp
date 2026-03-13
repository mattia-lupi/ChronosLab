#include "SymmetrizePattern.h"

// Symmetrizes the pattern of a non-symmetric matrix by padding with zeroes missing enries
int SymmetrizePattern(const int nrows, int *& iat, int *& ja, double *& coef,
                      double *&coef_T){

   // Transpose input pattern
   int *iat_T = nullptr;
   int *ja_T = nullptr;
   int ierr = transpose(nrows,nrows,iat,ja,coef,iat_T,ja_T,coef_T);
   if (ierr != 0) return 1;

   // Allocate room for the merged matrix
   int nterm = iat[nrows];
   int nterm_tmp = 2*nterm - nrows;
   int *iat_tmp = (int*) malloc((nrows+1) * sizeof(int));
   int *ja_tmp = (int*) malloc((nterm_tmp) * sizeof(int));
   double *coef_tmp = (double*) malloc((nterm_tmp) * sizeof(double));
   if (iat_tmp == nullptr || ja_tmp == nullptr || coef_tmp == nullptr) return 1;
   
   // Merge original and transposed patterns
   iat_tmp[0] = 0;
   int ind = 0;
   int iend   = iat[0];
   int iend_T = iat_T[0];
   for (int i = 0; i < nrows; i++){

      int len_out;

      int istart   = iend;
      iend     = iat[i+1];
      int len      = iend-istart;
      int istart_T = iend_T;
      iend_T   = iat_T[i+1];
      int len_T    = iend_T-istart_T;
      merge_row_patt(len,&(ja[istart]),&(coef[istart]),
                     len_T,&(ja_T[istart_T]),&(coef_T[istart_T]),
                     len_out,&(ja_tmp[ind]),&(coef_tmp[ind]));

      // Update pointer
      ind += len_out;
      iat_tmp[i+1] = ind;

   }

   // Free transposed matrix
   free(iat_T);
   free(ja_T);
   free(coef_T);

   // Backup input matrix
   int *iat_sav = iat;
   int *ja_sav = ja;
   double *coef_sav = coef;

   // Resize temporary storage
   nterm_tmp = ind;
   iat = iat_tmp;
   ja    = (int*) realloc( ja_tmp, (nterm_tmp) * sizeof(int));
   coef  = (double*) realloc( coef_tmp , (nterm_tmp) * sizeof(double));
   if ( ja == nullptr || coef == nullptr ) return 1;

   // Free initial storage of the input matrix
   free(iat_sav);
   free(ja_sav);
   free(coef_sav);

   // Transpose again to get coef_T
   ierr = transpose(nrows,nrows,iat,ja,coef,iat_T,ja_T,coef_T);
   if (ierr != 0) return 1;

   // Free unused data
   free(iat_T);
   free(ja_T);

   return 0;

}
