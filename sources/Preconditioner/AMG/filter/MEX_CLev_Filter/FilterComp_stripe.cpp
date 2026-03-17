////////////////////////////////////////////////
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>  // std::stringstream
#include <string>   // std::string
#define DEBUG false
//#define DEBUG true
////////////////////////////////////////////////

#include "FilterComp_stripe.h"

#define RCOND 1.e-12

int FilterComp_stripe(const double tau, const int shift, const int nrows,
                      const int nrows_A, int* iat_A, int* ja_A, double* coef_A,
                      const int nt_patt, const int* iat_patt, const int* ja_patt,
                      const int ntv, const double *const *TV, int &nt_AC_loc,
                      int *&iat_out, int *&ja_out, double *&coef_out){
   ////////////////////////////
   std::ofstream ofile;
   if (DEBUG){
   int myid = omp_get_thread_num();
   std::stringstream ss;
   ss << myid;
   std::string myid_label = ss.str();
   std::string logfile_name = "LOGFILE_" + myid_label + ".txt";
   ofile.open (logfile_name);
   }
   ////////////////////////////

   // Count the maximal number of entries per row
   int mmax = 0;
   int iend_A = iat_A[0];
   for (int irow = 0; irow < nrows; irow++){
      int istart_A = iend_A;
      iend_A = iat_A[irow+1];
      mmax = std::max(mmax,iend_A-istart_A);
   }

   // Set-up the length of the pattern
   bool patt_min = (nt_patt > 0);
   ////////////////////////////
   if (DEBUG){
      if (patt_min){
         std::cout << "PATT_MIN ERA TRUE" << std::endl;
      }
   }
   ////////////////////////////

   // Query workspace for dgelsy
   lapack_int *dummy_int;
   double *dummy_double;
   int optimal_lwork;
   lapack_int rank_out;
   double db_lwork;
   double rcond = RCOND;
   lapack_int ierr_lapack;
   ierr_lapack = LAPACKE_dgelsy_work(LAPACK_COL_MAJOR,static_cast<lapack_int>(ntv),
                 static_cast<lapack_int>(mmax),1,dummy_double,
                 static_cast<lapack_int>(ntv),dummy_double,
                 static_cast<lapack_int>(std::max(ntv,mmax)),dummy_int,rcond,&rank_out,
                 &db_lwork,-1);
   if (ierr_lapack != 0) return 2;
   optimal_lwork = static_cast<int>(db_lwork);
   ///////////////////////////////////////
   if (DEBUG){
   ofile << "mmax " << mmax << std::endl;
   ofile << "optimal_lwork " << optimal_lwork << std::endl;
   ofile.flush();
   }
   ///////////////////////////////////////

   // Allocate scratches
   lapack_int* JPVT = (lapack_int*) malloc(mmax*sizeof(lapack_int));
   double* full_SYS = (double*) malloc((mmax*ntv)*sizeof(double));
   double* lapack_WORK = (double*) malloc(optimal_lwork*sizeof(double));
   if (JPVT == nullptr || full_SYS == nullptr || lapack_WORK == nullptr) return 1;

   // Init pointer of the output matrix
   int ind_out = 0;
   iat_out[0] = 0;

   // Loop over rows
   iend_A = iat_A[0];
   for (int irow = 0; irow < nrows; irow++){

      int istart_A = iend_A;
      iend_A = iat_A[irow+1];
      int nt = iend_A - istart_A;
      // Compute the threshold for row filtering
      double threshold = 0.0;
      for (int i = istart_A; i < iend_A; i++) threshold += abs(coef_A[i]);
      threshold *= tau / 2.0;
      ///////////////////////////////////////
      if (DEBUG){
      ofile << "IROW " << irow << std::endl;
      ofile << "nt " << nt << std::endl;
      ofile << "threshold " << threshold << std::endl;
      ofile.flush();
      }
      ///////////////////////////////////////

      // Sort the row entries in absolute value
      abs_ri_heapsort(&(coef_A[istart_A]),&(ja_A[istart_A]),nt);
      ///////////////////////////////////////
      if (DEBUG){
      ofile << "ja_A  ";
      for (int i = 0; i < nt; i++) ofile << ja_A[istart_A+i] << " ";
      ofile << std::endl;
      for (int i = 0; i < nt; i++) ofile << coef_A[istart_A+i] << " ";
      ofile << std::endl;
      ofile.flush();
      }
      ///////////////////////////////////////
      
      // Count how many entries to remove
      int k_remove = -1;
      int ind_diag = -1;
      double trash = 0.0;
      while (trash < threshold){
         k_remove++;
         if (ja_A[istart_A+k_remove] == irow+shift){
            // Save the position of the diagonal entry
            ind_diag = istart_A+k_remove;
         } else {
            trash += abs(coef_A[istart_A+k_remove]);
         }
         //////////////////////////////////
         //ofile << "TRASH" << k_remove << " " << trash << std::endl;
         //////////////////////////////////
      }
      // Move the diagonal in the part to retain if it is not there yet
      if (ind_diag >= 0){
         k_remove--;
         swapi(ja_A[ind_diag],ja_A[k_remove]);
         swapr(coef_A[ind_diag],coef_A[k_remove]);
      }
      ///////////////////////////////////////
      if (DEBUG){
      if (ind_diag >= 0) ofile << "ind_diag POSITIVE" << ind_diag << std::endl;
      ofile << "k_remove " << k_remove << std::endl;
      ofile.flush();
      }
      ///////////////////////////////////////

      // Check the minimal pattern
      if (patt_min){
         // Sort by increasing column index
         ir_heapsort(&(ja_A[istart_A]),&(coef_A[istart_A]),k_remove);
         for (int i = iat_patt[irow]; i < iat_patt[irow+1]; i++){
            int jcol = ja_patt[i];
            int pos = bin_search<int,int>(jcol,k_remove,&(ja_A[istart_A]));
            /////////////////////////////////////////////
            if (DEBUG){
            ofile << "CERCO "<< jcol << "IN:" << std::endl;
            for (int j = istart_A; j < istart_A + k_remove; j++)
               ofile << " " << ja_A[j];
            ofile << std::endl;
            ofile << "TROVATO "<< pos << " " << ja_A[istart_A+pos] << std::endl;
            }
            /////////////////////////////////////////////
            if (ja_A[istart_A+pos] == jcol){
               // This entry belongs to the minimal pattern move it to the part to retain
               k_remove--;
               swapi(ja_A[istart_A+pos],ja_A[istart_A+k_remove]);
               swapr(coef_A[istart_A+pos],coef_A[istart_A+k_remove]);
            }
         }
      }

      // Sort entries in increasing column index order
      int k_retain = nt - k_remove;
      ir_heapsort(&(ja_A[istart_A+k_remove]),&(coef_A[istart_A+k_remove]),k_retain);
      ///////////////////////////////////////
      if (DEBUG){
      ofile << "k_remove / k_retain" << k_remove << " " << k_retain << std::endl;
      for (int i = 0; i < k_retain; i++) ofile << ja_A[istart_A+k_remove+i] << " ";
      ofile << std::endl;
      for (int i = 0; i < k_retain; i++) ofile << coef_A[istart_A+k_remove+i] << " ";
      ofile << std::endl;
      ofile.flush();
      ofile << "RHS CREATE: " << std::endl;
      }
      ///////////////////////////////////////

      // Create the rhs
      for (int i = 0; i < ntv; i++) coef_out[ind_out+i] = 0.0;
      for (int i = 0; i < k_remove; i++){
         int jcol = ja_A[istart_A+i];
         double fac = coef_A[istart_A+i];
         ///////////////////////////////////
         //ofile << "FAC: " << fac << " JCOL " << jcol << std::endl;
         //for (int i = 0; i < ntv; i++) ofile << TV[jcol][i] << " ";
         //ofile << std::endl;
         ///////////////////////////////////
         for (int i = 0; i < ntv; i++) coef_out[ind_out+i] += fac*TV[jcol][i];
      }

      // Create the dense overdetermined system
      double *pt_SYS = full_SYS;
      for (int i = 0; i < k_retain; i++){
         int jcol = ja_A[istart_A+k_remove+i];
         if (jcol == irow+shift){
            for (int i = 0; i < ntv; i++) pt_SYS[i] = 0.0;
         } else {
            for (int i = 0; i < ntv; i++) pt_SYS[i] = TV[jcol][i];
         }
         pt_SYS += ntv;
      }
      ///////////////////////////////////////
      if (DEBUG){
      ofile << "PRIMA DI LAPACK" << std::endl;
      ofile << "SYS:" << std::endl;
      for (int i = 0; i < ntv; i++){
         for (int j = 0; j < k_retain; j++){
            ofile << full_SYS[j*ntv+i] << " ";
         }
         ofile << std::endl;
      }
      ofile << "RHS:" << std::endl;
      for (int i = 0; i < ntv; i++) ofile << coef_out[ind_out+i] << std::endl;
      ofile << "optimal_lwork " << optimal_lwork << std::endl;
      ofile.flush();
      }
      ///////////////////////////////////////

      // Solve the overdetermined system
      for (int i = 0; i < k_retain; i++) JPVT[i] = 0;
      ierr_lapack = LAPACKE_dgelsy_work(LAPACK_COL_MAJOR,static_cast<lapack_int>(ntv),
                    static_cast<lapack_int>(k_retain),1,full_SYS,
                    static_cast<lapack_int>(ntv),&(coef_out[ind_out]),
                    static_cast<lapack_int>(std::max(k_retain,ntv)),JPVT,rcond,
                    &rank_out,lapack_WORK,optimal_lwork);
      ofile << "ierr_lapack " << ierr_lapack << "rank_out " << rank_out<< std::endl;
      if (ierr_lapack != 0) return 2;
      ////////////////////////////////////////////
      if (DEBUG){
      ofile << "SOL: " << k_retain << std::endl;
      for (int i = 0; i < k_retain; i++) ofile << coef_out[ind_out+i] << std::endl;
      ofile << "TRATTENUTI" << std::endl;
      }
      ////////////////////////////////////////////

      // Store this row
      for (int i = 0; i < k_retain; i++){
         ja_out[ind_out+i] = ja_A[istart_A+k_remove+i];
         coef_out[ind_out+i] += coef_A[istart_A+k_remove+i];
         ///////////////////////////////////////////////
         if (DEBUG){
         ofile << ja_out[ind_out+i] << " |  " << coef_out[ind_out+i] << std::endl;
         }
         ///////////////////////////////////////////////
      }
      ind_out += k_retain;
      iat_out[irow+1] = ind_out;

   }
   nt_AC_loc = ind_out;

   // Deallocate scratches
   free(JPVT);
   free(full_SYS);
   free(lapack_WORK);

   /////////////////////////////
   if (DEBUG){
   ofile << "----------------------------------------------" << std::endl;
   //for (int i = 0; i < nrows_A; i++){
      //for (int j = 0; j < ntv; j++){
         //ofile << TV[i][j] << " ";
      //}
      //ofile << std::endl;
   //}
   ofile.close();
   }
   /////////////////////////////

   return 0;
}
