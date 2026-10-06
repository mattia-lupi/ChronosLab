#include "cpt_nsy_rfsai.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <omp.h>
#include "KapGrad_NSY.h"
#include "gather_fullsys.h"
#include "inl_blas1.h"
#include "DEBUG.h"

#if defined(USE_MKL)
   #include "mkl_lapacke.h"
#elif defined(USE_OPENBLAS) || defined(__APPLE__)
   #include "lapacke.h"
#else
   #include "lapacke.h"
#endif

namespace {
struct FSAIThreadWorkspace {
    std::vector<int> JWN;
    std::vector<double> WR_L;
    std::vector<double> WR_U;
    std::vector<int> IWN_local;
    std::vector<double> full_A;
    std::vector<double> rhs_L;
    std::vector<double> rhs_U;
    std::vector<double> rhs_L_sav;
    std::vector<double> rhs_U_sav;
    std::vector<lapack_int> ipvt;

    void init(int nn, size_t max_m) {
        JWN.assign(nn, 0);
        WR_L.resize(nn);
        WR_U.resize(nn);
        IWN_local.resize(nn);
        full_A.resize(max_m * max_m);
        rhs_L.resize(max_m);
        rhs_U.resize(max_m);
        rhs_L_sav.resize(max_m);
        rhs_U_sav.resize(max_m);
        ipvt.resize(max_m);
    }
};
}

int cpt_nsy_rfsai(const int nstep, const int step_size, const double eps,
                  const int nn_A, const int nt_A, const double *diag_A,
                  const int *iat_A, const int *ja_A, const double *coef_A,
                  const double *coef_AT, int *&iat_FL, int *&ja_FL,
                  double *&coef_FL, double *&coef_FUT,
                  const int num_threads){
   (void)nt_A;
   Open_DebugLog();

   int mmax = nstep * step_size;

   // Compute exact maximum capacity per row: deg(row_i) + mmax + 1
   size_t max_mrow = (size_t)mmax + 1;
   std::vector<size_t> row_start_offset(nn_A + 1, 0);
   for (int i = 0; i < nn_A; i++) {
      size_t max_row_cap = (size_t)(iat_A[i+1] - iat_A[i]) + (size_t)mmax + 1;
      row_start_offset[i+1] = row_start_offset[i] + max_row_cap;
      if (max_row_cap > max_mrow) max_mrow = max_row_cap;
   }
   size_t total_temp_cap = row_start_offset[nn_A];

   // Pre-allocate temporary row buffers for lock-free parallel computation
   int *row_nnz = (int*) calloc(nn_A, sizeof(int));
   int *ja_FL_temp = (int*) malloc(total_temp_cap * sizeof(int));
   double *coef_FL_temp = (double*) malloc(total_temp_cap * sizeof(double));
   double *coef_FUT_temp = (double*) malloc(total_temp_cap * sizeof(double));

   if (row_nnz == nullptr || ja_FL_temp == nullptr ||
       coef_FL_temp == nullptr || coef_FUT_temp == nullptr) {
      if (row_nnz) free(row_nnz);
      if (ja_FL_temp) free(ja_FL_temp);
      if (coef_FL_temp) free(coef_FL_temp);
      if (coef_FUT_temp) free(coef_FUT_temp);
      return 1;
   }

   std::vector<FSAIThreadWorkspace> workspaces(num_threads);
   for (int t = 0; t < num_threads; ++t) {
      workspaces[t].init(nn_A, max_mrow);
   }

   int global_ierr = 0;

   // Loop over the rows in parallel
   #pragma omp parallel for schedule(guided) num_threads(num_threads)
   for (int irow = 0; irow < nn_A; irow++){
      if (global_ierr != 0) continue;

      int tid = omp_get_thread_num();
      FSAIThreadWorkspace &ws = workspaces[tid];

      int *JWN = ws.JWN.data();
      double *WR_L = ws.WR_L.data();
      double *WR_U = ws.WR_U.data();
      int *IWN = ws.IWN_local.data();
      double *full_A = ws.full_A.data();
      double *rhs_L = ws.rhs_L.data();
      double *rhs_U = ws.rhs_U.data();
      double *rhs_L_sav = ws.rhs_L_sav.data();
      double *rhs_U_sav = ws.rhs_U_sav.data();
      lapack_int *ipvt = ws.ipvt.data();

      // Loop for the refinement of the row pattern
      int mrow = 0;
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
         KapGrad_NSY(istep, irow, mrow, irow, step_size, iat_A, ja_A, coef_A, coef_AT,
                     rhs_L, rhs_U, IWN, JWN, WR_L, WR_U);

         // Compute the F_L and F_U rows if the pattern is not null
         if (mrow > mrow_old){

            // Gather the coefficients of the full local systems
            bool null_L = true;
            bool null_U = true;
            gather_fullsys(irow,mrow,IWN,nn_A,iat_A,ja_A,coef_A,full_A,
                           rhs_L,rhs_U,null_L,null_U);
            //////////////////////////////////////////////////////////
            if (DEBUG){
               fprintf(dbfile,"full_A:\n");
               for (int i = 0; i < mrow; i++){
                  for (int j = 0; j < mrow; j++) fprintf(dbfile," %15.6e",full_A[j*mrow+i]);
                  fprintf(dbfile,"\n");
               }
               fprintf(dbfile,"JCOLS: ");
               for (int i = 0; i < mrow; i++) fprintf(dbfile," %15d",IWN[i]);
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
               lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, mrow, mrow, full_A, mrow, ipvt);
               if (info != 0) {
                  #pragma omp atomic write
                  global_ierr = 3;
                  break;
               }
            }

            // Backup rhs
            for (int k = 0; k < mrow; k++) rhs_L_sav[k] = rhs_L[k];
            for (int k = 0; k < mrow; k++) rhs_U_sav[k] = rhs_U[k];

            // Solve linear systems for F_L row and F_U column
            if (!null_L){
               lapack_int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'T', mrow, 1, full_A, mrow, ipvt, rhs_L, 1);
               if (info != 0) {
                  #pragma omp atomic write
                  global_ierr = 3;
                  break;
               }
            }
            if (!null_U){
               lapack_int info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', mrow, 1, full_A, mrow, ipvt, rhs_U, 1);
               if (info != 0) {
                  #pragma omp atomic write
                  global_ierr = 3;
                  break;
               }
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
            double DKap_L_new = inl_ddot(mrow, rhs_L, 1, rhs_U_sav, 1);
            double DKap_U_new = inl_ddot(mrow, rhs_U, 1, rhs_L_sav, 1);
            double DKap_new = fabs(DKap_L_new + DKap_U_new);

            // Exit check
            if (istep == nstep){
               Refine = false;
            } else {
               Refine = (fabs(DKap_new - DKap_old) >= eps * DKap_old) && (DKap_new != 0.0);
               DKap_old = fabs(DKap_new);
            }
         } else {
            // Pattern is empty: row is uncoupled
            Refine = false;
         }
      } // end refinement loop

      if (global_ierr != 0) continue;

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

      size_t row_dest = row_start_offset[irow];
      double fac = 1.0 / sqrt(fabs(scal_fac));

      // Scale lower part
      for (int k = 0; k < mrow; k++) {
         ja_FL_temp[row_dest + k] = IWN[k];
         coef_FL_temp[row_dest + k] = fac * rhs_L[k];
      }
      // Store diagonal entry
      ja_FL_temp[row_dest + mrow] = irow;
      coef_FL_temp[row_dest + mrow] = fac;

      // Scale upper part
      double fac_U = (scal_fac < 0.0) ? -fac : fac;
      for (int k = 0; k < mrow; k++) {
         coef_FUT_temp[row_dest + k] = fac_U * rhs_U[k];
      }
      // Store diagonal entry
      coef_FUT_temp[row_dest + mrow] = fac_U;

      // Reset the non-zero indicator
      for (int k = 0; k < mrow; k++) JWN[IWN[k]] = 0;

      row_nnz[irow] = mrow + 1;
   } // end parallel row loop

   if (global_ierr != 0) {
      free(row_nnz);
      free(ja_FL_temp);
      free(coef_FL_temp);
      free(coef_FUT_temp);
      return global_ierr;
   }

   // Build final contiguous CSR structure
   iat_FL = (int*) malloc((nn_A + 1) * sizeof(int));
   if (iat_FL == nullptr) {
      free(row_nnz); free(ja_FL_temp); free(coef_FL_temp); free(coef_FUT_temp);
      return 1;
   }

   iat_FL[0] = 0;
   for (int i = 0; i < nn_A; i++) {
      iat_FL[i + 1] = iat_FL[i] + row_nnz[i];
   }

   int nt_F = iat_FL[nn_A];
   ja_FL    = (int*) malloc(nt_F * sizeof(int));
   coef_FL  = (double*) malloc(nt_F * sizeof(double));
   coef_FUT = (double*) malloc(nt_F * sizeof(double));

   if (ja_FL == nullptr || coef_FL == nullptr || coef_FUT == nullptr) {
      free(row_nnz); free(ja_FL_temp); free(coef_FL_temp); free(coef_FUT_temp);
      return 1;
   }

   // Parallel copy from row buffers to packed output arrays
   #pragma omp parallel for schedule(guided) num_threads(num_threads)
   for (int i = 0; i < nn_A; i++) {
      size_t src_offset = row_start_offset[i];
      int dst_offset = iat_FL[i];
      int cnt = row_nnz[i];

      for (int k = 0; k < cnt; k++) {
         ja_FL[dst_offset + k] = ja_FL_temp[src_offset + k];
         coef_FL[dst_offset + k] = coef_FL_temp[src_offset + k];
         coef_FUT[dst_offset + k] = coef_FUT_temp[src_offset + k];
      }
   }

   // Free temporary buffers
   free(row_nnz);
   free(ja_FL_temp);
   free(coef_FL_temp);
   free(coef_FUT_temp);

   // Close DEBUG log
   Close_DebugLog();

   return 0;
}
