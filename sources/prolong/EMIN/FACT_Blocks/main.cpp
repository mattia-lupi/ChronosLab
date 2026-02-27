#include <iostream>  // to use: cout,endl
#include <iomanip>
#include <stdlib.h>  // to use: exit
#include <fstream>   // to use: ifstream,ofstream
#include <sstream>   // to use: stringstream
#include <math.h>
#include <libgen.h>  // to use: basename
#include <chrono>    // to use: system_clock,duration
#include "lapacke.h"  // to use: dpotrf,dpotrs
using namespace std;

#include "readMat.h"
#include "wrCSRmat.h"
#include "ICHOL_wrapper.h"
#include "icholRF_apply.h"
#define charStrLen 2000

#define LAPACK_SOLVE false
#define ICHOL_SOLVE true

// MAIN PROGRAM
int main(int argc, const char* argv[]){

FILE *catFILE;

// Other variables
bool BINREAD;
int ncat, cat_size;
int ierr, kk;
int nn_A, nt_A;
int *iat_A, *ja_A;
double *coef_A;

char line[1024];

   // Check arguments
   if (argc < 3) {
      printf("Too few arguments.\n Usage: [cat file] [matrix file]\n");
      exit(1);
   }

   // Read categories
   int nblk;
   int* pt_cat;
   int* pt_blk;
   if ( (catFILE = fopen(argv[1], "r")) != NULL) {
      bool read_err = false;
      char *fget_ret;
      fget_ret = fgets(line, sizeof line, catFILE); sscanf(line, "%d",  &cat_size);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, catFILE); sscanf(line, "%d",  &ncat);
      read_err = read_err || (fget_ret == nullptr);
      // Allocate pt_cat
      pt_cat = (int*) malloc((ncat+1)*sizeof(int));
      fget_ret = fgets(line, sizeof line, catFILE);
      int kk = 0;
      int nlines = (ncat+1)/10; if ((ncat+1)%10 > 0) nlines++;
      for (int i = 0; i < nlines; i ++){
         fget_ret = fgets(line, sizeof line, catFILE);
         string s(line);
         stringstream ss(s);
         int next_number;
         while (ss >> next_number) pt_cat[kk++] = next_number;
      }
      nblk =  pt_cat[ncat] - 1;
      // Allocate pt_blk
      pt_blk = (int*) malloc((nblk+1)*sizeof(int));
      fget_ret = fgets(line, sizeof line, catFILE);
      kk = 0;
      nlines = (nblk+1)/10; if ((nblk+1)%10 > 0) nlines++;
      for (int i = 0; i < nlines; i ++){
         fget_ret = fgets(line, sizeof line, catFILE);
         string s(line);
         stringstream ss(s);
         int next_number;
         while (ss >> next_number) pt_blk[kk++] = next_number;
      }
      if (read_err){
         printf(" Error while opening file %s\n", argv[1]);
         exit(1);
      }
   } else {
      printf(" Error while opening file %s\n", argv[1]);
      exit(0);
   }

   // Dump input
   cout << "cat_size:  " << cat_size << endl;
   cout << "ncat:      " << ncat << endl;
   cout << "pt_cat:" << endl;
   for (int i = 0; i < ncat+1; i++){
      cout << " " << pt_cat[i];
      if ( (i+1)%10 == 0 ) cout << endl;
   }
   cout << endl;
   cout << "pt_blk:" << endl;
   for (int i = 0; i < nblk+1; i++){
      cout << " " << pt_blk[i];
      if ( (i+1)%10 == 0 ) cout << endl;
   }
   cout << endl;

   // Adapt pointers to C
   for (int i = 0; i < ncat+1; i++) pt_cat[i]--;
   for (int i = 0; i < nblk+1; i++) pt_blk[i]--;

   // Read the matrix
   BINREAD = false;
   ierr = readCSRmat(&nn_A, &kk, &nt_A, &iat_A, &ja_A, &coef_A, argv[2], BINREAD);
   if (ierr != 0) exit(ierr);

   cout << "Matrix read" << endl;

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

   // Extract max size and max nnz
   int nnmax = 0;
   int ntmax = 0;
   int istart, iend, kstart, kend;
   iend = pt_blk[0];
   kend = iat_A[iend];
   for (int iblk = 0; iblk < nblk; iblk++){
      istart = iend;
      iend = pt_blk[iblk+1];
      kstart = kend;
      kend = iat_A[iend];
      nnmax = max(nnmax,iend-istart);
      ntmax = max(ntmax,kend-kstart);
   }
   cout << "nnmax: " << nnmax << endl;
   cout << "ntmax: " << ntmax << endl;

   // Allocate rhs
   double *rhs    = (double*) malloc(nnmax*sizeof(double));
   double *sol    = (double*) malloc(nnmax*sizeof(double));
   for (int i = 0; i < nnmax; i++) rhs[i] = 1.0;

   // --- Local variables for time printing ----------------------------------------------
   chrono::duration<double> elapsed_seconds;

   // *** Solve all the blocks with LAPACK ***********************************************
   double *DPOTRF_secs;
   double *DPOTRS_secs;
   if (LAPACK_SOLVE){

      cout << endl;
      cout << "SOLVING WITH LAPACK" << endl;

      // Allocate dense matrix and rhs
      double *full_A = (double*) malloc((nnmax*nnmax)*sizeof(double));

      DPOTRF_secs = (double*) malloc(ncat*sizeof(double));
      DPOTRS_secs = (double*) malloc(ncat*sizeof(double));
      for (int i = 0; i < ncat; i++){
         DPOTRF_secs[i] = 0.0;
         DPOTRS_secs[i] = 0.0;
      }
      iend = pt_cat[0];
      for (int icat = 0; icat < ncat; icat++){
         istart = iend;
         iend = pt_cat[icat+1];
         kend = pt_blk[istart];
         for (int iblk = istart; iblk < iend; iblk++){

            if (iblk%100 == 0) cout << "Block " << setw(10) << iblk << " out of " << nblk << endl;

            kstart = kend;
            kend = pt_blk[iblk+1];

            //---START-------------------------------
            auto start = chrono::system_clock::now();
            //---------------------------------------

            // 1) Load the sparse block into the dense block
            int nn_loc = kend - kstart;
            fill_n(full_A,nn_loc*nn_loc,0.0);
            double *pt_full_A = full_A;
            for (int irow = kstart; irow < kend; irow++){
               for (int j = iat_A[irow]; j < iat_A[irow+1]; j++){
                  int jcol = ja_A[j] - kstart;
                  pt_full_A[jcol] = coef_A[j];
               }
               pt_full_A += nn_loc;
            }

            // 2) Factorize the dense system
            lapack_int info;
            info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',nn_loc,full_A,nn_loc);
            //if(info != 0){throw linsol_error ("cpt_aFSAIcoef","error in LAPACKE_dpotrf");}

            //---STOP--------------------------------
            auto end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            DPOTRF_secs[icat] += elapsed_seconds.count();
            //---------------------------------------

            //---START-------------------------------
            start = chrono::system_clock::now();
            //---------------------------------------

            // 3) Solve the dense system
            info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',nn_loc,1,full_A,nn_loc,rhs,nn_loc);
            //if(info != 0){throw linsol_error ("cpt_aFSAIcoef","error in LAPACKE_dpotrs");}

            //---STOP--------------------------------
            end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            DPOTRS_secs[icat] += elapsed_seconds.count();
            //---------------------------------------

         }
      }

      // Free dense A
      free(full_A);
   }

   // *** Solve all the blocks with ICHOL ************************************************
   double *ICHOLF_secs;
   double *ICHOLS_secs;
   if (ICHOL_SOLVE){

      cout << endl;
      cout << "SOLVING WITH ICHOL" << endl;

      // Allocate local matrix and incomplete factorization
      int *iat_loc     = (int*) malloc((nnmax+1)*sizeof(int));
      int *ja_loc      = (int*) malloc(ntmax*sizeof(int));
      double *coef_loc = (double*) malloc(ntmax*sizeof(double));
      int lfil = 0;
      int iwk_U = ntmax + nnmax * (lfil+1);
      int *it_U     = (int*) malloc((nnmax+1)*sizeof(int));
      int *jcol_U      = (int*) malloc(iwk_U*sizeof(int));
      double *coef_U = (double*) malloc(iwk_U*sizeof(double));
      double *D_inv = (double*) malloc(nnmax*sizeof(double));
      int ireg_scr_size = nnmax + iwk_U;
      int *ireg_scr = (int*) malloc(ireg_scr_size*sizeof(int));
      int iext_scr_size = 4*nnmax + iwk_U;
      int *iext_scr = (int*) malloc(iext_scr_size*sizeof(int));

      ICHOLF_secs = (double*) malloc(ncat*sizeof(double));
      ICHOLS_secs = (double*) malloc(ncat*sizeof(double));
      for (int i = 0; i < ncat; i++){
         ICHOLF_secs[i] = 0.0;
         ICHOLS_secs[i] = 0.0;
      }
      iend = pt_cat[0];
      for (int icat = 0; icat < ncat; icat++){
         istart = iend;
         iend = pt_cat[icat+1];
         kend = pt_blk[istart];
         for (int iblk = istart; iblk < iend; iblk++){

            if (iblk%100 == 0) cout << "Block " << setw(10) << iblk << " out of " << nblk << endl;

            kstart = kend;
            kend = pt_blk[iblk+1];

            //---START-------------------------------
            auto start = chrono::system_clock::now();
            //---------------------------------------

            // 1) Load the sparse block into the local sparse block
            //    (moving for C to Fortran style)
            int nn_loc = kend - kstart;
            int ind = 0;
            iat_loc[0] = 1;
            for (int irow = kstart; irow < kend; irow++){
               for (int j = iat_A[irow]; j < iat_A[irow+1]; j++){
                  int jcol = ja_A[j];
                  if (jcol >= irow){
                     // Load only if it belongs to the upper part
                     ja_loc[ind] = jcol - kstart + 1;
                     coef_loc[ind] = coef_A[j];
                     ind++;
                  }
               }
               iat_loc[irow-kstart+1] = ind + 1;
            }
            int nt_loc = ind;
            ////////////////////////////////////////
            //if (iblk == 0){
            //        cout << nt_loc << endl;
            //   FILE *otmp = fopen("XXXX","w");
            //   for (int irow = 0; irow < nn_loc; irow++){
            //      for (int j = iat_loc[irow]; j < iat_loc[irow+1]; j++){
            //         fprintf(otmp,"%6d %6d %15.6e\n",irow+1,ja_loc[j],coef_loc[j]);
            //      }
            //   }
            //   fflush(otmp);
            //   fclose(otmp);
            //}
            ////////////////////////////////////////

            // 2) Factorize the local sparse system
            int jcol_offset = 0;
            /////////////////////////////////////////
            //cout << "PRIMA DI CHIAMARE:" << endl;
            //cout << nn_loc << " " << nt_loc << " max iwk_U " << iwk_U << " curr iwk_U "
            //     << nt_loc + (lfil+1)*nn_loc << endl;
            /////////////////////////////////////////
            ierr = ICHOL_wrapper(lfil,jcol_offset,nn_loc,nt_loc,ireg_scr_size,iext_scr_size,
                                 iat_loc,ja_loc,coef_loc,it_U,jcol_U,coef_U,D_inv,
                                 ireg_scr,iext_scr);

            //---STOP--------------------------------
            auto end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            ICHOLF_secs[icat] += elapsed_seconds.count();
            //---------------------------------------

            //---START-------------------------------
            start = chrono::system_clock::now();
            //---------------------------------------

            // 3) Solve the sparse system
            icholRF_apply(nn_loc,it_U,jcol_U,coef_U,D_inv,rhs,sol);

            //---STOP--------------------------------
            end = chrono::system_clock::now();
            elapsed_seconds = end - start;
            ICHOLS_secs[icat] += elapsed_seconds.count();
            //---------------------------------------

         }
      }

      // Free scratches
      free(iat_loc);
      free(ja_loc);
      free(coef_loc);
      free(it_U);
      free(jcol_U);
      free(coef_U);
      free(D_inv);
      free(ireg_scr);
      free(iext_scr);

   }

   // Compute total times
   double DPOTRF_tot = 0.0;
   double DPOTRS_tot = 0.0;
   if (LAPACK_SOLVE){
      for (int icat = 0; icat < ncat; icat++) DPOTRF_tot += DPOTRF_secs[icat];
      for (int icat = 0; icat < ncat; icat++) DPOTRS_tot += DPOTRS_secs[icat];
   }
   double ICHOLF_tot = 0.0;
   double ICHOLS_tot = 0.0;
   if (ICHOL_SOLVE){
      for (int icat = 0; icat < ncat; icat++) ICHOLF_tot += ICHOLF_secs[icat];
      for (int icat = 0; icat < ncat; icat++) ICHOLS_tot += ICHOLS_secs[icat];
   }

   // Dump results on file
   FILE *times = fopen("time.log","w");
   fprintf(times,"  %s  %s","Size up to","# of Blocks");
   if (LAPACK_SOLVE) fprintf(times,"  %s  %s","LAPACK fact [s]","LAPACK solve [s]");
   if (ICHOL_SOLVE) fprintf(times,"  %s  %s","ICHOL fact [s]","ICHOL solve [s]");
   fprintf(times,"\n");
   for (int icat = 0; icat < ncat; icat++){
      fprintf(times," %11d %12d",(icat+1)*cat_size,pt_cat[icat+1]-pt_cat[icat]);
      if (LAPACK_SOLVE) fprintf(times,"  %15.6f   %15.6f",DPOTRF_secs[icat],DPOTRS_secs[icat]);
      if (ICHOL_SOLVE) fprintf(times," %15.6f  %15.6f",ICHOLF_secs[icat],ICHOLS_secs[icat]);
      fprintf(times,"\n");
   }
   fprintf(times,"%25s","********* TOTAL *********");
   if (LAPACK_SOLVE) fprintf(times,"  %15.6f   %15.6f",DPOTRF_tot,DPOTRS_tot);
   if (ICHOL_SOLVE) fprintf(times," %15.6f  %15.6f",ICHOLF_tot,ICHOLS_tot);
   fprintf(times,"\n");
   fclose(times);


   if (LAPACK_SOLVE){
      cout << endl << "TIMES WITH LAPACK" << endl;
      for (int icat = 0; icat < ncat; icat++)
         cout << DPOTRF_secs[icat] << " " << DPOTRS_secs[icat] << endl;
   }


   if (ICHOL_SOLVE){
      cout << endl << "TIME WITH ICHOL" << endl;
      for (int icat = 0; icat < ncat; icat++)
         cout << ICHOLF_secs[icat] << " " << ICHOLS_secs[icat] << endl;
   }

   exit(0);

}
