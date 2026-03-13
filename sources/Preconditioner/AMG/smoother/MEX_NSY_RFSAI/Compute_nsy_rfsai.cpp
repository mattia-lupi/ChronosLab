#include "Compute_nsy_rfsai.h"

int Compute_nsy_rfsai(const int nstep, const int step_size, const double eps,
                      const int nn_A, const int *iat_A, const int *ja_A,
                      const double *coef_A, int *&iat_FL, int *&ja_FL,
                      double *&coef_FL, int *&iat_FU, int *&ja_FU, double *&coef_FU){

   // Init error code
   int ierr = 0;

   // Print the input matrix
   if (false){
      cout << "Printing input Matrix" << endl;
      FILE *ofile = fopen("mat_input.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_A[i]; j < iat_A[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_A[j]+1,coef_A[j]);
          }
      }
      fclose(ofile);
   }

   // --- Extract diagonal entries -------------------------------------------------------

   double *diag_A = nullptr;
   extract_diag(nn_A,iat_A,ja_A,coef_A,diag_A);

   // --- Symmetrize input matrix --------------------------------------------------------

   int nt_A = iat_A[nn_A];
   int *iat_Asym     = (int*) malloc((nn_A+1) * sizeof(int));
   int *ja_Asym      = (int*) malloc(nt_A * sizeof(int));
   double *coef_Asym = (double*) malloc(nt_A * sizeof(double));
   if (iat_Asym == nullptr || ja_Asym == nullptr || coef_Asym == nullptr){
      cout << "Error in allocating symmetrize pattern A" << endl;
      return 1;
   }
   memcpy(iat_Asym,iat_A,(nn_A+1)*sizeof(int)); 
   memcpy(ja_Asym,ja_A,nt_A*sizeof(int)); 
   memcpy(coef_Asym,coef_A,nt_A*sizeof(double)); 

   double *coef_AT  = nullptr;
   ierr = SymmetrizePattern(nn_A,iat_Asym,ja_Asym,coef_Asym,coef_AT);
   if (ierr != 0){
      cout << "Error in SymmetrizePattern" << endl;
      return 2;
   }

   // Print the input matrix
   if (false){
      cout << "Printing transposed input Matrix" << endl;
      FILE *ofile = fopen("mat_transposed.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_Asym[i]; j < iat_Asym[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_Asym[j]+1,coef_AT[j]);
          }
      }
      fclose(ofile);
      cout << "Printing symmetrized input Matrix" << endl;
      ofile = fopen("mat_symmetrized.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_Asym[i]; j < iat_Asym[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_Asym[j]+1,coef_Asym[j]);
          }
      }
      fclose(ofile);
   }

   // --- Compute the non-symmetric FSAI -------------------------------------------------

   ierr = cpt_nsy_rfsai(nstep,step_size,eps,nn_A,nt_A,diag_A,iat_Asym,ja_Asym,coef_Asym,
                        coef_AT,iat_FL,ja_FL,coef_FL,coef_FU);
   if (ierr != 0){
      cout << "Error in cpt_nsy_rfsai" << endl;
      return 3;
   }

   if (false){
      FILE *pfile = fopen("mat_FL.csr","w"); if (!pfile) exit(1);
      wrCSRmat(pfile,false,nn_A,iat_FL,ja_FL,coef_FL);
      fclose(pfile);
      pfile = fopen("mat_FU.csr","w"); if (!pfile) exit(1);
      wrCSRmat(pfile,false,nn_A,iat_FL,ja_FL,coef_FU);
      fclose(pfile);
   }

   // --- Compress the preconditioner (transposing FU) -----------------------------------

   int nt_FL;
   int nt_FU;
   ierr = compress_nsy_fsai(nn_A,nt_FL,nt_FU,iat_FL,ja_FL,iat_FU,ja_FU,coef_FL,coef_FU);
   if (ierr != 0){
      cout << "Error in compress_nsy_fsai" << endl;
      return 4;
   }

   // Delete scratch vectors
   free(diag_A);
   free(iat_Asym);
   free(ja_Asym);
   free(coef_Asym);
   free(coef_AT);

   return 0;

}
