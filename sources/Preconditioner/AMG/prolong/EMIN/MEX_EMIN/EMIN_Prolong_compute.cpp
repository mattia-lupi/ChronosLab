//----------------------------------------------------------------------------------------

#if defined PRINT
#define dump true
#else
#define dump false
#endif

#include <iostream>  // to use: cout,endl
using namespace std;

#include "mex.h"
#include "EMIN_ImpProl.h"
#include "DebEnv.h"

//----------------------------------------------------------------------------------------
//
// MATLAB function
//
// function [iat_Pout,ja_Pout,coef_Pout,info] =
//          EMIN_Prolong_compute(np,itmax,ener_tol,condmax,prec,sol_type,min_lfil,
//                               max_lfil,D_lfil,nn,nn_C,ntv,nt_A,nt_P,nt_patt,fcnode,
//                               iat_A,ja_A,coef_A,iat_Pin,ja_Pin,coef_Pin,iat_patt,
//                               ja_patt,TV);
//
//----------------------------------------------------------------------------------------

// Gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

   // Get the value of the input scalars
   if (dump) mexPrintf("*** EMIN_Prolong_compute ***\n");
   if (dump) mexPrintf("- Get input scalars\n");
   int level      = mxGetScalar(prhs[ 0]);
   int np         = mxGetScalar(prhs[ 1]);
   int itmax      = mxGetScalar(prhs[ 2]);
   double en_tol  = mxGetScalar(prhs[ 3]);
   double condmax = mxGetScalar(prhs[ 4]);
   double maxwgt  = mxGetScalar(prhs[ 5]);
   int prec       = mxGetScalar(prhs[ 6]);
   int sol_type   = mxGetScalar(prhs[ 7]);
   int min_lfil   = mxGetScalar(prhs[ 8]);
   int max_lfil   = mxGetScalar(prhs[ 9]);
   int D_lfil     = mxGetScalar(prhs[10]);
   int nn         = mxGetScalar(prhs[11]);
   int nn_C       = mxGetScalar(prhs[12]);
   int ntv        = mxGetScalar(prhs[13]);
   int nt_A       = mxGetScalar(prhs[14]);
   int nt_P       = mxGetScalar(prhs[15]);
   int nt_patt    = mxGetScalar(prhs[16]);

   // Get the value of the input arrays
   if (dump) mexPrintf("- Get input arrays\n");
   int *fcnode       = (int*)    mxGetData(prhs[17]);
   int *iat_A        = (int*)    mxGetData(prhs[18]);
   int *ja_A         = (int*)    mxGetData(prhs[19]);
   double *coef_A    = (double*) mxGetData(prhs[20]);
   int *iat_Pin      = (int*)    mxGetData(prhs[21]);
   int *ja_Pin       = (int*)    mxGetData(prhs[22]);
   double *coef_Pin  = (double*) mxGetData(prhs[23]);
   int *iat_patt     = (int*)    mxGetData(prhs[24]);
   int *ja_patt      = (int*)    mxGetData(prhs[25]);
   double *TVbuf     = (double*) mxGetData(prhs[26]);
   
   // Allocate TV
   double **TV = (double**) malloc( nn*sizeof(double*) );
   if (TV == nullptr){
      cout << "ERROR in allocating TV" << endl;
      return;
   }
   int k = 0;
   for (int i = 0; i < nn; i++){
      TV[i] = &(TVbuf[k]);
      k += ntv;
   }

   // Init debug
   if (level == 1){
      DebEnv.SetDebEnv(np,"w");
   } else {
      DebEnv.OpenDebugLog("a");
   }
   if (DEBUG){
      fprintf(DebEnv.r_logfile,"\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n",level);
      fflush(DebEnv.r_logfile);
      for (int i = 0; i < np; i++){
         fprintf(DebEnv.t_logfile[i],"\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n",level);
         fflush(DebEnv.t_logfile[i]);
      }
   }

   // Compute Pout
   if (dump) mexPrintf("- Compute Pout entries\n");
   double info[EMIN_INFO_SZ];
   int *iat_Pout, *ja_Pout;
   double *coef_Pout;
   int ierr = EMIN_ImpProl(np,itmax,en_tol,condmax,maxwgt,prec,sol_type,min_lfil,max_lfil,
                           D_lfil,nn,nn_C,ntv,nt_A,nt_P,nt_patt,fcnode,iat_A,ja_A,
                           coef_A,iat_Pin,ja_Pin,coef_Pin,iat_patt,ja_patt,TV,
                           iat_Pout,ja_Pout,coef_Pout,info);
   if (ierr != 0){
      cout << "ERROR in EMIN_ImpProl with code: " << ierr << endl;
      return;
   }
   free(TV);

   // Close debug log
   DebEnv.CloseDebugLog();

   // Resize output arrays
   if (dump) mexPrintf("- Store Pout into the output arrays\n");

   int *pt_i;
   double *pt_r;

   int nt_Pout = iat_Pout[nn];
   // iat_Pout
   mxArray *iat_final = mxCreateNumericMatrix((mwSize) 1,(mwSize) (nn+1),
                                                mxINT32_CLASS,mxREAL);
   pt_i = (int*) mxGetData(iat_final);
   for (int k = 0; k < nn+1; k++) pt_i[k] = iat_Pout[k]+1;
   free(iat_Pout);
   // ja_Pout
   mxArray *ja_final   = mxCreateNumericMatrix((mwSize) 1,(mwSize) nt_Pout,
                                                mxINT32_CLASS,mxREAL);
   pt_i = (int*) mxGetData(ja_final);
   for (int k = 0; k < nt_Pout; k++) pt_i[k] = ja_Pout[k]+1;
   free(ja_Pout);
   // coef_Pout
   mxArray *coef_final = mxCreateNumericMatrix((mwSize) 1,(mwSize) nt_Pout,
                                                mxDOUBLE_CLASS,mxREAL);
   pt_r = (double*) mxGetData(coef_final);
   for (int k = 0; k < nt_Pout; k++) pt_r[k] = coef_Pout[k];
   free(coef_Pout);

   // Store info
   mxArray *info_mex = mxCreateNumericMatrix((mwSize) 1,(mwSize) EMIN_INFO_SZ,
                                              mxDOUBLE_CLASS,mxREAL);
   pt_r = (double*) mxGetData(info_mex);
   pt_r[0]  = info[0];
   pt_r[1]  = info[1];
   pt_r[2]  = info[2];
   pt_r[3]  = info[3];
   pt_r[4]  = info[4];
   pt_r[5]  = info[5];
   pt_r[6]  = info[6];
   pt_r[7]  = info[7];
   pt_r[8]  = info[8];
   pt_r[9]  = info[9];
   pt_r[10] = info[10];
   pt_r[11] = info[11];

   // Return pointers
   plhs[0] = iat_final;
   plhs[1] = ja_final;
   plhs[2] = coef_final;
   plhs[3] = info_mex;

   if (dump) mexPrintf("\n");

}

//----------------------------------------------------------------------------------------
