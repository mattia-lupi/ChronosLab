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

#include "DebEnv.h"
#include "readMat.h"
#include "wrCSRmat.h"
#include "EMIN_ImpProl.h"
#define charStrLen 2000
#define DUMP_MATIN false
/*
 * This driver takes the following input:
 *
 * np: number of openMP threads
 * BINREAD:                 read matrix and pattern as binary or ascii file
 * itmax:                   max number of CG iterations for energy minimization
 * System matrix:           iat_A, ja_A, coef_A
 * Initial prolongation:    iat_P, ja_P, coef_P
 * Prolongation pattern:    iat_patt, ja_patt
 * Test space:              TV(nn,ntv)
 * fcnode:                  fine/coarse partition
 */

// MAIN PROGRAM
int main(int argc, const char* argv[]){

FILE *parmFILE;

// Input variables
bool BINREAD;
int np;
int itmax;
double condmax;
double en_tol;
double maxwgt;
int prec_type, sol_type;
int min_lfil, max_lfil, D_lfil;
int ierr, k, kk;
char MAT_file[charStrLen]="";
char PROL_file[charStrLen]="";
char PATT_file[charStrLen]="";
char TV_file[charStrLen]="";
char FC_file[charStrLen]="";

char line[1024];

   // Check arguments
   if (argc < 2) {
      printf("Too few arguments.\n Usage: [parm file]\n");
      exit(1);
   }

   // Read parameters
   if ( (parmFILE = fopen(argv[1], "r")) != NULL) {
      bool read_err = false;
      char *fget_ret;
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &kk);
      read_err = read_err || (fget_ret == nullptr);
      BINREAD = kk != 0;
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &np);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &itmax);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%lf",  &en_tol);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%lf",  &condmax);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%lf",  &maxwgt);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &prec_type);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &sol_type);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &min_lfil);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &max_lfil);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%d",  &D_lfil);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%s",  MAT_file);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%s",  PROL_file);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%s",  PATT_file);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%s",  TV_file);
      read_err = read_err || (fget_ret == nullptr);
      fget_ret = fgets(line, sizeof line, parmFILE); sscanf(line, "%s",  FC_file);
      read_err = read_err || (fget_ret == nullptr);
      if (read_err){
         printf(" Error while opening file %s\n", argv[1]);
         exit(1);
      }
   } else {
      printf(" Error while opening file %s\n", argv[1]);
      exit(0);
   }

   // Dump input
   cout << "BINREAD:   " << BINREAD << endl;
   cout << "NP:        " << np << endl;
   cout << "ITMAX:     " << itmax << endl;
   cout << "EN. Toll.: " << en_tol << endl;
   cout << "CONDMAX:   " << condmax << endl;
   cout << "MAXWGT:    " << maxwgt << endl;
   cout << "PREC_TYPE: " << prec_type << endl;
   cout << "SOL_TYPE:  " << sol_type << endl;
   cout << "MIN_LFIL:  " << min_lfil << endl;
   cout << "MAX_LFIL:  " << max_lfil << endl;
   cout << "D_LFIL:    " << D_lfil << endl;
   cout << "MAT_file:  " << MAT_file << endl;
   cout << "PROL_file: " << PROL_file << endl;
   cout << "PATT_file: " << PATT_file << endl;
   cout << "TV_file:   " << TV_file << endl;
   cout << "FC_file:   " << FC_file << endl;
   cout << endl;

   // Set variables for DEBUG
   DebEnv.SetDebEnv(np,"w");

   // Read the matrix
   int nn_A, nt_A;
   int *iat_A, *ja_A;
   double *coef_A;
   ierr = readCSRmat(&nn_A, &kk, &nt_A, &iat_A, &ja_A, &coef_A, MAT_file, BINREAD);
   if (ierr != 0) exit(ierr);

   cout << "Matrix read" << endl;

   // Print the input matrix
   if (DUMP_MATIN){
      cout << "Printing input Matrix" << endl;
      FILE *ofile = fopen("mat_input.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_A; i++){
          for (int j = iat_A[i]; j < iat_A[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_A[j]+1,coef_A[j]);
          }
      }
      fclose(ofile);
   }

   // Read initial prolongation
   int nn_P, nt_P;
   int *iat_P, *ja_P;
   double *coef_P;
   ierr = readCSRmat(&nn_P, &kk, &nt_P, &iat_P, &ja_P, &coef_P, PROL_file, BINREAD);
   if (ierr != 0) exit(ierr);

   cout << "Prolongation read" << endl;

   // Print the input prolongation
   if (DUMP_MATIN){
      cout << "Printing input Prolongation" << endl;
      FILE *ofile = fopen("prol_input.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_P; i++){
          for (int j = iat_P[i]; j < iat_P[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_P[j]+1,coef_P[j]);
          }
      }
      fclose(ofile);
   }

   // Read patten
   int nn_patt, nt_patt;
   int *iat_patt, *ja_patt;
   double *coef_tmp;
   ierr = readCSRmat(&nn_patt, &kk, &nt_patt, &iat_patt, &ja_patt, &coef_tmp,
                     PATT_file, BINREAD);
   if (ierr != 0) exit(ierr);

   cout << "Pattern read" << endl;

   // Print the input pattern
   if (DUMP_MATIN){
      cout << "Printing input pattern" << endl;
      FILE *ofile = fopen("patt_input.csr","w"); if (!ofile) exit(1);
      for (int i = 0; i < nn_patt; i++){
          for (int j = iat_patt[i]; j < iat_patt[i+1]; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_patt[j]+1,coef_tmp[j]);
          }
      }
      fclose(ofile);
   }
   free(coef_tmp);

   // Read the test space
   int ntv;
   // Open the input file
   ifstream fileTV(TV_file);
   // Read header
   k = 0;
   fileTV >> kk >> ntv;
   if (kk != nn_A) exit(1);
   // Allocate and read TV matrix
   double **TV   = (double**) malloc( nn_A*sizeof(double*) );
   double *TVbuf = (double*) malloc( (ntv*nn_A)*sizeof(double) );
   for (int i = 0; i < nn_A; i++){
      TV[i] = &(TVbuf[k]);
      k += ntv;
      for (int j = 0; j < ntv; j++) fileTV >> TV[i][j];
   }
   // Close the input fileTV
   fileTV.close();
   cout << "TSPACE read " << ntv << endl;

   // Read fcnode
   // Open the input file
   ifstream fileFC(FC_file);
   // Allocate fcnode
   int *fcnode   = (int*) malloc( nn_A*sizeof(int) );
   for (int i = 0; i < nn_A; i++) fileFC >> fcnode[i];
   // Close the input file
   fileFC.close();
   cout << "FCNODE read " << endl;

   // Count coarse nodes and adjust indices of coarse nodes
   int nn_C = 0;
   for (int i = 0; i < nn_A; i++){
      if (fcnode[i] > 0) {
         fcnode[i]--;
         nn_C++;
      }
   }
   cout << "nn_C: " << nn_C << endl;

   // --- Local variables for time printing ----------------------------------------------
   chrono::duration<double> elapsed_seconds;

   // *** Improve prolongation through Energy Minimization *******************************

   //---START-------------------------------
   auto start = chrono::system_clock::now();
   //---------------------------------------

   int *iat_Pnew;
   int *ja_Pnew;
   double *coef_Pnew;
   double emin_info[EMIN_INFO_SZ];
   ierr = EMIN_ImpProl(np,itmax,en_tol,condmax,maxwgt,prec_type,sol_type,min_lfil,max_lfil,
                       D_lfil,nn_A,nn_C,ntv,nt_A,nt_P,nt_patt,fcnode,iat_A,ja_A,coef_A,
                       iat_P,ja_P,coef_P,iat_patt,ja_patt,TV,iat_Pnew,ja_Pnew,coef_Pnew,
                       emin_info,false);

   //---STOP--------------------------------
   auto end = chrono::system_clock::now();
   elapsed_seconds = end - start;
   //---------------------------------------

   cout << endl << endl;
   cout << "Time for Enengy Minimization [sec]: " << elapsed_seconds.count() << endl;
   cout << endl;

   //cout << "Printing Final Prolongation" << endl;
   //FILE *pfile = fopen("Prol_EMIN.csr","w"); if (!pfile) exit(1);
   //wrCSRmat(pfile,false,nn_P,iat_Pnew,ja_Pnew,coef_Pnew);
   //fclose(pfile);

   exit(0);

}
