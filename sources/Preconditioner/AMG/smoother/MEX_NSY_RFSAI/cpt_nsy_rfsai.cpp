// NOTES: there is the need to access the rows of the lower part of A and the columns of
// the upper part of A. The ideal storage is with low(A) in CSR and upp(A) in CSC

#include "cpt_nsy_rfsai.h"

int cpt_nsy_rfsai(const int nstep, const int step_size, const double eps, const int nn_A,
                  const int nt_A, const double *diag_A, const int *iat_A, const int *ja_A,
                  const double *coef_A, const double *coef_AT, int *&iat_FL, int *&ja_FL,
                  double *&coef_FL, double *&coef_FUT){

   // Open DEBUG log
   Open_DebugLog();

   // Allocate room for the preconditioner
   int mmax = nstep*step_size;
   int ntmax_F = nn_A*(mmax+1) + nn_A;
   iat_FL   = (int*) malloc((nn_A+1) * sizeof(int));
   ja_FL    = (int*) malloc((ntmax_F) * sizeof(int));
   coef_FL  = (double*) malloc((ntmax_F) * sizeof(double));
   coef_FUT = (double*) malloc((ntmax_F) * sizeof(double));
   if ( iat_FL == nullptr ||  ja_FL == nullptr ||
        coef_FL == nullptr || coef_FUT == nullptr ) return 1;

   // Allocate scratches
   // JW: Non-zero indicator for retained entries
   int *JWN = (int*) malloc(nn_A * sizeof(int));
   double *WR_L = (double*) malloc(nn_A * sizeof(double));
   double *WR_U = (double*) malloc(nn_A * sizeof(double));
   if ( JWN == nullptr || WR_L == nullptr || WR_U == nullptr ) return 2;
   // Scratches for local dense systems
   double *full_A = (double*) malloc( (mmax*mmax) * sizeof(double));
   double *rhs_L = (double*) malloc( (mmax+1) * sizeof(double));
   double *rhs_U = (double*) malloc( (mmax+1) * sizeof(double));
   lapack_int *ipvt = (lapack_int*) malloc( mmax * sizeof(lapack_int));
   if (full_A == nullptr || rhs_L == nullptr || rhs_U == nullptr || ipvt == nullptr)
      return 2;
   double *rhs_L_sav = (double*) malloc( (mmax+1) * sizeof(double));
   double *rhs_U_sav = (double*) malloc( (mmax+1) * sizeof(double));
   if ( rhs_L_sav == nullptr || rhs_U_sav == nullptr ) return 2;

   // Init JWN
   std::fill_n(JWN,nn_A,0);

   // Initialize pointer to the beginning of the row
   int ind_FL = 0;
   iat_FL[0] = ind_FL;

   // Loop over the rows of the current processor
   for( int irow = 0; irow < nn_A; irow++){

      //////////////////////////////////////////////////////////
      if (DEBUG){
         fprintf(dbfile,"-------------------------------\n");
         fprintf(dbfile,"IROW: %d\n",irow);
      }
      //////////////////////////////////////////////////////////

      int mrow = 0;

      // Loop for the refinement of the row pattern
      int istep = 0;
      double DKap_old = 0.0;
      bool Refine = (nstep >= 1);
      while (Refine){

         istep++;
         //////////////////////////////////////////////////////////
         if (DEBUG) fprintf(dbfile,"istep %6d mroww %6d\n",istep,mrow);
         //////////////////////////////////////////////////////////

         // Compute the Kaporin gradient
         int mrow_old = mrow;
         KapGrad_NSY(istep,irow,mrow,irow,step_size,iat_A,ja_A,coef_A,coef_AT,rhs_L,rhs_U,
                     &(ja_FL[ind_FL]),JWN,WR_L,WR_U);

         // Compute the F_L and F_U rows if the pattern is not null
         if (mrow > mrow_old){

            // Gather the coefficients of the full local systems
            bool null_L;
            bool null_U;
            gather_fullsys(irow,mrow,&(ja_FL[ind_FL]),nn_A,iat_A,ja_A,coef_A,full_A,
                           rhs_L,rhs_U,null_L,null_U);
            //////////////////////////////////////////////////////////
            if (DEBUG){
               fprintf(dbfile,"full_A:\n");
               for (int i = 0; i < mrow; i++){
                  for (int j = 0; j < mrow; j++) fprintf(dbfile," %15.6e",full_A[j*mrow+i]);
                  fprintf(dbfile,"\n");
               }
               fprintf(dbfile,"JCOLS: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15d",ja_FL[ind_FL+i]);
               fprintf(dbfile,"\n");
               fprintf(dbfile,"RHS_L: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L[i]);
               fprintf(dbfile,"\n");
               fprintf(dbfile,"RHS_U: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U[i]);
               fprintf(dbfile,"\n");
            }
            //////////////////////////////////////////////////////////

            // Factorize the dense matrix
            if (!null_L || !null_U){
               lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,mrow,mrow,full_A,mrow,
                                                ipvt);
               if (info != 0) return 3; //@@@@@@@@@@@@@@@@@
               //@@@@@@@@@@@@@@@@@@@@@@@@@@
               // GESTIONE ERRORE LAPACK
               /*
               if (info < 0) {
                   printf("cpt_aFSAIcoef: LAPACK ERROR %d FOR ROW %d\n",info,irow);
                  throw linsol_error ("cpt_aFSAIcoef","error in LAPACKE_dpotrf");
               } else if(info > 0) {
                  if (DEBUG){
                     type_OMP_int myid = omp_get_thread_num();
                     fprintf(DebEnv.t_logfile[myid],"LAPACK ERROR %d FOR ROW %d\n",info,irow);
                  }
                  // Solve with a smaller number of non-zeroes
                  int k = 0;
                  for (int i = 0; i < mrow; i++){
                     // Remove entries added in the last step
                     IWN[i-k] = IWN[i]; // IWN ora e SOVRAPPOSTO &(ja_FL[ind_FL])
                     if (JWN[IWN[i]-1] == -istep){
                        k++;
                        JWN[IWN[i]-1] = 0;
                     }
                  }
                  mrow -= k;
                  // Gather the system again
                  GatherFullSys(nulrhs,irow,mrow,irow,nequ,nterm,mmax,iat,ja,IWN,
                                coef_A,full_A,rhs);
                  // Factorize
                  info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',mrow,full_A,mmax);
                  // Save rhs
                  vec_rhs_sav = vec_rhs;
                  // Backward and forward substitution
                  info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',mrow,1,full_A,mmax,rhs,
                                        max(int(1),mrow));
                  // Exit the refinement loop
                  goto exit_Refinement;
      
               }
               */
               //@@@@@@@@@@@@@@@@@@
      
            }

            // Backup system and rhs
            //for (int k = 0; k < mrow*mrow; k++) full_A_sav[k] = full_A[k];
            for (int k = 0; k < mrow; k++) rhs_L_sav[k] = rhs_L[k];
            for (int k = 0; k < mrow; k++) rhs_U_sav[k] = rhs_U[k];

            // Compute coefficients of the irow-th row of FL / column of FU
            if (!null_L){
               lapack_int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'T',mrow,1,full_A,mrow,
                                                ipvt,rhs_L,1);
               if (info != 0) return 3; //@@@@@@@@@@@@@@@@@
            }
            if (!null_U){
               lapack_int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',mrow,1,full_A,mrow,
                                                ipvt,rhs_U,1);
               if (info != 0) return 3; //@@@@@@@@@@@@@@@@@
            }
            //////////////////////////////////////////////////////////
            if (DEBUG){
               fprintf(dbfile,"SOL_L: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L[i]);
               fprintf(dbfile,"\n");
               fprintf(dbfile,"SOL_U: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U[i]);
               fprintf(dbfile,"\n");
            }
            //////////////////////////////////////////////////////////

            // Compute the Kaporin number decrease
            double DKap_L_new = inl_ddot(mrow,rhs_L,1,rhs_U_sav,1);
            double DKap_U_new = inl_ddot(mrow,rhs_U,1,rhs_L_sav,1);
            double DKap_new = fabs(DKap_L_new + DKap_U_new);

            // Exit check
            if (istep == nstep){
               Refine = false;
            } else {
               Refine = (fabs(DKap_new-DKap_old) >= eps*DKap_old) && (DKap_new != 0.);
               DKap_old = fabs(DKap_new);
            } 

         } else {

            // If the pattern is empty the row is uncoupled
            Refine = false;

         }

      } // end refinement loop

      // Exit point
      exit_Refinement: ;

      // Compute the scaling factor for this row
      double diag_entry = diag_A[irow];
      double scal_fac = diag_entry - 0.5*( inl_ddot(mrow,rhs_U,1,rhs_L_sav,1) +
                                           inl_ddot(mrow,rhs_L,1,rhs_U_sav,1) );
      //////////////////////////////////////////////////////////
      if (DEBUG){
         fprintf(dbfile,"mrow: %d\n",mrow);
         fprintf(dbfile,"RHS_L: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"RHS_U: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"RHS_L_SAV: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_L_sav[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"RHS_U_SAV: ");
         for (int i = 0; i < mrow; i++) fprintf(dbfile," %15.6e",rhs_U_sav[i]);
         fprintf(dbfile,"\n");
         fprintf(dbfile,"diag_entry: %15.6e\n",diag_entry);
         fprintf(dbfile,"f*b: %15.6e\n",inl_ddot(mrow,rhs_U,1,rhs_L_sav,1));
         fprintf(dbfile,"c*g: %15.6e\n",inl_ddot(mrow,rhs_L,1,rhs_U_sav,1));
         fprintf(dbfile,"SCAL FACTOR: %15.6e\n",scal_fac);
      }
      //////////////////////////////////////////////////////////

      // Check zero diagonal
      double check_val = fabs(scal_fac / diag_entry);
      if ( check_val < 1.0e-10 ){
         std::cout << "SMALL DIAGONAL = " << check_val << " IN ROW: " << irow << std::endl;
      }

      // Scale lower part
      double fac = 1.0 / sqrt(fabs(scal_fac));
      for (int k = 0; k < mrow; k++) coef_FL[ind_FL+k] = fac*rhs_L[k];
      // Store diagonal entry
      coef_FL[ind_FL+mrow] = fac;
      ja_FL[ind_FL+mrow] = irow;

      // Scale upper part
      if (scal_fac < 0.0) fac = -fac;
      for (int k = 0; k < mrow; k++) coef_FUT[ind_FL+k] = fac*rhs_U[k];
      // Store diagonal entry
      coef_FUT[ind_FL+mrow] = fac;

      // Reset the non-zero indicator
      for (int k = 0; k < mrow; k++) JWN[ja_FL[ind_FL+k]] = 0;

      // Set the pointer to the beginning of next row
      ind_FL += mrow + 1;
      iat_FL[irow+1] = ind_FL;

   } // end row loop

   // Free scratches
   free(JWN);
   free(WR_L);
   free(WR_U);
   free(full_A);
   free(rhs_L);
   free(rhs_U);
   free(ipvt);
   free(rhs_L_sav);
   free(rhs_U_sav);

   // Reallocate ja_FL, coef_FL and coef_FUT with their true length
   int nt_F = ind_FL;
   ja_FL    = (int*) realloc( ja_FL, (nt_F) * sizeof(int));
   coef_FL  = (double*) realloc( coef_FL , (nt_F) * sizeof(double));
   coef_FUT = (double*) realloc( coef_FUT , (nt_F) * sizeof(double));
   if ( iat_FL == nullptr ||  ja_FL == nullptr || 
        coef_FL == nullptr || coef_FUT == nullptr ) return 1;

   // Close DEBUG log
   Close_DebugLog();

   return 0;
}
