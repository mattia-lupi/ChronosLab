/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 * outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include <stdlib.h>
////////////////////////////////
#include <iostream>
////////////////////////////////
using namespace std;

#include "mex.h"
#include "BAMG.h"
#include "BAMG_params.h"
#include "DebEnv.h"

/* The computational routine wrapper*/
int cpt_Prolongation_BAMG(const int level, const BAMG_params params, const int nthreads,
                          const int nn_S, const int nt_S, const int *const iat_S,
                          const int *const ja_S, const int *const coef_S, const int ntv,
                          const int *const fcnodes, const double *const *const TV,
                          const int nc_I, int &nt_I, vector<int> &vec_iat_I,
                          vector<int> &vec_ja_I, vector<double> &vec_coef_I,
                          vector<int> &vec_c_mark){

     // Init debug
     if (level == 1){
        DebEnv.SetDebEnv(nthreads,"w");
     } else {
        DebEnv.OpenDebugLog("a");
     }
     if (DEBUG && BAMG_DEBUG){
        fprintf(DebEnv.r_logfile,"\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n",level);
        fflush(DebEnv.r_logfile);
        for (int i = 0; i < nthreads; i++){
           fprintf(DebEnv.t_logfile[i],"\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n",level);
           fflush(DebEnv.t_logfile[i]);
        }
     }

     // Set nn_L to 0 as it is useless without MPI
     iReg nn_L = 0;
     iReg nn_C = nn_S;
     int ierr = BAMG(params,nthreads,nn_L,nn_C,nn_S,iat_S,ja_S,ntv,fcnodes,TV,
                     nt_I,vec_iat_I,vec_ja_I,vec_coef_I,vec_c_mark);

     // Close debug log
     DebEnv.CloseDebugLog();

   return ierr;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

// Prototype of MATLAB calling
// [nt_I,iat_I,ja_I,coef_I,c_mark] = cpt_Prolongation_BAMG(level,np,stronly,itmax_vol,...
//                                   dist_min,dist_max,mmax,maxcond,maxrownrmd,tol_vol,...
//                                   eps_prol,nn_S,nt_S,iat_S,ja_S,coef_S,ntv,fcnodes,TV,
//                                   nn_I,nc_I);

   // Check for proper number of arguments (Inputs and Outputs)
   if (nrhs != 20) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","22 inputs required.");
   }
   if(nlhs!=5) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongations","5 outputs required.");
   }

   // Translate inputs
   int    level      =           mxGetScalar(prhs[ 0]);
   int    np         =           mxGetScalar(prhs[ 1]);
   int    itmax_vol  =           mxGetScalar(prhs[ 2]);
   int    dist_min   =           mxGetScalar(prhs[ 3]);
   int    dist_max   =           mxGetScalar(prhs[ 4]);
   int    mmax       =           mxGetScalar(prhs[ 5]);
   double maxcond    =           mxGetScalar(prhs[ 6]);
   double maxrownrm  =           mxGetScalar(prhs[ 7]);
   double tol_vol    =           mxGetScalar(prhs[ 8]);
   double eps_prol   =           mxGetScalar(prhs[ 9]);
   int    nn_S       =           mxGetScalar(prhs[10]);
   int    nt_S       =           mxGetScalar(prhs[11]);
   int    *iat_S     = (int*)    mxGetPr(    prhs[12]);
   int    *ja_S      = (int*)    mxGetPr(    prhs[13]);
   int    *coef_S    = (int*)    mxGetPr(    prhs[14]);
   int    ntv        =           mxGetScalar(prhs[15]);
   int    *fcnodes   = (int*)    mxGetPr(    prhs[16]);
   double *TV        = (double*) mxGetPr(    prhs[17]);
   int    nn_I       =           mxGetScalar(prhs[18]);
   int    nc_I       =           mxGetScalar(prhs[19]);

   // Store input in the structure
   BAMG_params params;
   params.verbosity = VERB_LEV;
   params.itmax_vol = itmax_vol;
   params.dist_min = dist_min;
   params.dist_max = dist_max;
   params.mmax = mmax;
   params.maxcond = maxcond;
   params.maxrownrm = maxrownrm;
   params.tol_vol = tol_vol;
   params.eps = eps_prol;

   // Load TV in a 2D buffer
   double **TV_2D = (double**) malloc(nn_S*sizeof(double*));
   if (TV_2D == NULL) mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","TV_2D allocation");
   int kk = 0;
   for (int i = 0; i < nn_S; i++){
      TV_2D[i] = (double*) malloc(ntv*sizeof(double));
      if (TV_2D[i] == NULL) mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","TV_2D allocation");
      for (int j = 0; j < ntv; j++) {
         TV_2D[i][j] = TV[kk];
         kk++;
      }
   }
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
   //cout << "DENTRO LEVEL   " << level << endl;
   //cout << "PARAMETRI: " << 
   //" VERB: " << params.verbosity <<
   //" ITMAX_VOL: " << params.itmax_vol <<
   //" DIST_MIN: " << params.dist_min <<
   //" DIST_MAX: " << params.dist_max <<
   //" MMAX: " << params.mmax <<
   //" MAXCOND: " << params.maxcond <<
   //" MAXROWNRM: " << params.maxrownrm <<
   //" TOL_VOL: " << params.tol_vol <<
   //" EPS: " << params.eps << endl;
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
   
   // Call the computational routine
   int nt_I;
   vector<int> vec_iat_I,vec_ja_I;
   vector<double> vec_coef_I;
   vector<int> vec_c_mark;
   int ierr = cpt_Prolongation_BAMG(level,params,np,nn_S,nt_S,iat_S,ja_S,coef_S,
                                    ntv,fcnodes,TV_2D,nc_I,nt_I,vec_iat_I,vec_ja_I,
                                    vec_coef_I,vec_c_mark);
 
   // Copy results in MATLAB (@@@ Typed Data Access NOT working on RUSSEL)

   // Set some pointers
   int    *iat_I  = vec_iat_I.data();
   int    *ja_I   = vec_ja_I.data();
   double *coef_I = vec_coef_I.data();
   int    *c_mark = vec_c_mark.data();

   // nt_I
   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
   double *out_nt_I = mxGetPr(plhs[0]);
   *out_nt_I = (double) nt_I;

   // iat_I
   plhs[1] = mxCreateDoubleMatrix((mwSize) nn_I+1,1,mxREAL);
   double *out_iat_I = mxGetPr(plhs[1]);
   for (int i = 0; i < nn_I+1; i++) out_iat_I[i] = (double) (iat_I[i]+1);

   // ja_I
   plhs[2] = mxCreateDoubleMatrix((mwSize) nt_I,1,mxREAL);
   double *out_ja_I = mxGetPr(plhs[2]);
   for (int i = 0; i < nt_I; i++) out_ja_I[i] = (double) (ja_I[i]+1);

   // coef_I
   plhs[3] = mxCreateDoubleMatrix((mwSize) nt_I,1,mxREAL);
   double *out_coef_I = mxGetPr(plhs[3]);
   for (int i = 0; i < nt_I; i++) out_coef_I[i] = coef_I[i];

   // c_mark
   plhs[4] = mxCreateDoubleMatrix((mwSize) nn_I,1,mxREAL);
   double *out_c_mark = mxGetPr(plhs[4]);
   for (int i = 0; i < nn_I; i++) out_c_mark[i] = (double) (c_mark[i]);

   // Free temporarily allocated arrays
   for (int i = 0; i < nn_S; i++) free(TV_2D[i]); free(TV_2D);

}
