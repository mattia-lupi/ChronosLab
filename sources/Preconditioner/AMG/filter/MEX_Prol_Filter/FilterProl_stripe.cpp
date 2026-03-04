////////////////////////////////////////////////
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>  // std::stringstream
#include <string>   // std::string
using namespace std;
#define DEBUG false
//#define DEBUG true
//////////////////////////////////////////////////

#include "FilterProl_stripe.h"

#define RCOND 1.e-12

int FilterProl_stripe(const double perc, const double tol, const int firstrow,
                      const int nrows, const int nn_P, int *iat_P, int *ja_P,
                      double *coef_P, const int ntv, const double *const *TV,
                      int &nt_PF_loc, int *&iat_PF, int *&ja_PF, double *&coef_PF){

   ////////////////////////////
   ofstream ofile;
   if (DEBUG){
   int myid = omp_get_thread_num();
   stringstream ss;
   ss << myid;
   string myid_label = ss.str();
   string logfile_name = "PLOGFILE_" + myid_label + ".txt";
   ofile.open (logfile_name);
   }
   ////////////////////////////

   // Count the maximal number of entries per row
   int mmax = 0;
   int iend_P = iat_P[0];
   for (int irow = 0; irow < nrows; irow++){
      int istart_P = iend_P;
      iend_P = iat_P[irow+1];
      mmax = max(mmax,iend_P-istart_P);
   }

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
                 static_cast<lapack_int>(max(ntv,mmax)),dummy_int,rcond,&rank_out,
                 &db_lwork,-1);
   if (ierr_lapack != 0) return 2;
   optimal_lwork = static_cast<int>(db_lwork);
   ///////////////////////////////////////
   if (DEBUG){
   ofile << "mmax " << mmax << endl;
   ofile << "optimal_lwork " << optimal_lwork << endl;
   ofile.flush();
   }
   ///////////////////////////////////////

   // Allocate scratches
   lapack_int* JPVT = (lapack_int*) malloc(mmax*sizeof(lapack_int));
   double* full_SYS = (double*) malloc((mmax*ntv)*sizeof(double));
   double* lapack_WORK = (double*) malloc(optimal_lwork*sizeof(double));
   if (JPVT == nullptr || full_SYS == nullptr || lapack_WORK == nullptr) return 1;

   // Init pointer of the output matrix
   int ind_PF = 0;
   iat_PF[0] = 0;

   // Loop over rows
   double perc_100 = perc / 100.0;
   iend_P = iat_P[0];
   for (int irow = 0; irow < nrows; irow++){

      int istart_P = iend_P;
      iend_P = iat_P[irow+1];
      int nt = iend_P - istart_P;

      // Sort the row entries in absolute value
      abs_ri_heapsort(&(coef_P[istart_P]),&(ja_P[istart_P]),nt);

      // Compute the absolute norm of the row
      double row_nrm = abs_norm(nt,&(coef_P[istart_P]));

      // Compute threshold (sum of the entries to retain)
      double threshold = perc_100*row_nrm;

      double cum_sum = 0.0;
      int k_remove = nt;
      while (cum_sum < threshold){
         k_remove--;
         cum_sum += abs(coef_P[istart_P+k_remove]);
      }
      int k_retain = nt - k_remove;
      if (k_retain < ntv){
         k_retain = min(nt,ntv);
         k_remove = nt - k_retain;
      }
      ///////////////////////////////////////
      if (DEBUG){
      ofile << "IROW " << irow << endl;
      ofile << "row_nrm " << row_nrm << endl;
      ofile << "threshold " << threshold << endl;
      ofile << "k_remove " << k_remove << endl;
      ofile << "k_retain " << k_retain << endl;
      ofile << "nt " << nt << endl;
      }
      ///////////////////////////////////////

      // Compute a correction whose norm is sufficiently small
      double abs_tol = tol*row_nrm;
      while (true){

         // Exit if there are no entries to remove
         if (k_remove == 0) {
            fill_n(&(coef_PF[ind_PF]),nt,0.0);
            break;
         }

         // Create the rhs
         for (int i = 0; i < ntv; i++) coef_PF[ind_PF+i] = 0.0;
         for (int i = 0; i < k_remove; i++){
            int jcol = ja_P[istart_P+i];
            double fac = coef_P[istart_P+i];
            ///////////////////////////////////
            //ofile << "FAC: " << fac << " JCOL " << jcol << endl;
            //for (int i = 0; i < ntv; i++) ofile << TV[jcol][i] << " ";
            //ofile << endl;
            ///////////////////////////////////
            for (int i = 0; i < ntv; i++) coef_PF[ind_PF+i] += fac*TV[jcol][i];
         }

         // Create the dense overdetermined system
         double *pt_SYS = full_SYS;
         for (int i = k_remove; i < nt; i++){
            int jcol = ja_P[istart_P+i];
            for (int i = 0; i < ntv; i++) pt_SYS[i] = TV[jcol][i];
            pt_SYS += ntv;
         }

         ///////////////////////////////////////
         if (DEBUG){
         ofile << "PRIMA DI LAPACK" << endl;
         ofile << "SYS:" << endl;
         for (int i = 0; i < ntv; i++){
            for (int j = 0; j < (k_retain); j++){
               ofile << full_SYS[j*ntv+i] << " ";
            }
            ofile << endl;
         }
         ofile << "RHS:" << endl;
         for (int i = 0; i < ntv; i++) ofile << coef_PF[ind_PF+i] << endl;
         ofile << "optimal_lwork " << optimal_lwork << endl;
         ofile.flush();
         }
         ///////////////////////////////////////

         // Solve the overdetermined system
         for (int i = 0; i < k_retain; i++) JPVT[i] = 0;
         ierr_lapack = LAPACKE_dgelsy_work(LAPACK_COL_MAJOR,static_cast<lapack_int>(ntv),
                       static_cast<lapack_int>(k_retain),1,full_SYS,
                       static_cast<lapack_int>(ntv),&(coef_PF[ind_PF]),
                       static_cast<lapack_int>(max(k_retain,ntv)),JPVT,rcond,
                       &rank_out,lapack_WORK,optimal_lwork);
         ofile << "ierr_lapack " << ierr_lapack << "rank_out " << rank_out<< endl;
         if (ierr_lapack != 0) return 2;
         ////////////////////////////////////////////
         if (DEBUG){
         ofile << "SOL: " << k_retain << endl;
         for (int i = 0; i < k_retain; i++) ofile << coef_PF[ind_PF+i] << endl;
         }
         ////////////////////////////////////////////

         // Compute the norm of the correction
         double nrm_h = abs_norm(k_retain,&(coef_PF[ind_PF]));

         ////////////////////////////////////////////
         if (DEBUG){
            ofile << "nrm_h: " << nrm_h << endl;
            ofile << "nrm_h " << nrm_h << " " << abs_tol << endl;
         }
         ////////////////////////////////////////////

         // Check if the solution is correct and do not alter previous row_norm
         if ( nrm_h > abs_tol || rank_out < ntv){
            int k_add = 1 + static_cast<int>(0.3*static_cast<double>(k_retain));
            k_retain = min(k_add+k_retain,nt);
            k_remove = nt - k_retain;
            ////////////////////////////////////////////
            if (DEBUG){
               ofile << "AUMENTO" << endl;
               ofile << "k_retain new: " << k_retain << endl;
            }
            ////////////////////////////////////////////
         } else {
            break;
         }

      }

      // Store this row
      for (int i = 0; i < k_retain; i++){
         ja_PF[ind_PF+i] = ja_P[istart_P+k_remove+i];
         coef_PF[ind_PF+i] += coef_P[istart_P+k_remove+i];
         ///////////////////////////////////////////////
         if (DEBUG){
         ofile << ja_PF[ind_PF+i] << " |  " << coef_PF[ind_PF+i] << endl;
         }
         ///////////////////////////////////////////////
      }
      // Sort the row
      ir_heapsort(&(ja_PF[ind_PF]),&(coef_PF[ind_PF]),k_retain);
      ind_PF += k_retain;
      iat_PF[irow+1] = ind_PF;

   }
   nt_PF_loc = ind_PF; 

   // Free scratch vectors
   free(JPVT);
   free(full_SYS);
   free(lapack_WORK);

   return 0;

}
