/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include <stdlib.h>

#include "mex.h"
#include "EXTI_prolongation.h"

/* The computational routine wrapper*/
int cpt_Prolongation_EXTI(const int level, const int nthreads,
                          const int *const vecstart, const int nn_A, const int nt_A,
                          const int *const iat_A, const int *const ja_A,
                          const double *const coef_A,const int *const coef_S,
                          const int *const iat_C, const int *const ja_C,
                          const double *const coef_C, const int *const fcnodes,
                          const int nr_I, const int nc_I, int &nt_I, int *&iat_I,
                          int *&ja_I, double *&coef_I){

   int ierr = 0;
   ierr = EXTI_prolongation(level,nthreads,vecstart,nn_A,nt_A,iat_A,ja_A,
                            coef_A,coef_S,iat_C,ja_C,coef_C,fcnodes,nr_I,nc_I,
                            nt_I,iat_I,ja_I,coef_I);

   return ierr;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

// Prototype of MATLAB calling
//
// [nt_I,iat_I,ja_I,coef_I] = cpt_Prolongation_EXTI(level,np,exti_type,distance,...
//                            vecstart,nn_A,nt_A,iat_A,ja_A,coef_A,coef_S,...
//                            iat_C,ja_C,coef_C,fcnodes,nn_I,nc_I);

   // Check for proper number of arguments (Inputs and Outputs)
   if (nrhs != 15) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongation","14 inputs required.");
   }
   if(nlhs!=4) {
       mexErrMsgIdAndTxt("ERROR in cpt_Prolongations","4 outputs required.");
   }

   // Translate inputs
   int    level      =           mxGetScalar(prhs[ 0]);
   int    np         =           mxGetScalar(prhs[ 1]);
   int    *vecstart  = (int*)    mxGetPr(    prhs[ 2]);
   int    nn_A       =           mxGetScalar(prhs[ 3]);
   int    nt_A       =           mxGetScalar(prhs[ 4]);
   int    *iat_A     = (int*)    mxGetPr(    prhs[ 5]);
   int    *ja_A      = (int*)    mxGetPr(    prhs[ 6]);
   double *coef_A    = (double*) mxGetPr(    prhs[ 7]);
   int    *coef_S    = (int*)    mxGetPr(    prhs[ 8]);
   int    *iat_C     = (int*)    mxGetPr(    prhs[ 9]);
   int    *ja_C      = (int*)    mxGetPr(    prhs[10]);
   double *coef_C    = (double*) mxGetPr(    prhs[11]);
   int    *fcnodes   = (int*)    mxGetPr(    prhs[12]);
   int    nn_I       =           mxGetScalar(prhs[13]);
   int    nc_I       =           mxGetScalar(prhs[14]);

   // Call the computational routine
   int nt_I;
   int *iat_I,*ja_I;
   double *coef_I;
   int ierr = cpt_Prolongation_EXTI(level,np,vecstart,nn_A,nt_A,iat_A,ja_A,coef_A,
		                            coef_S,iat_C,ja_C,coef_C,fcnodes,nn_I,nc_I,nt_I,
                                    iat_I,ja_I,coef_I);
 
   // Copy results in MATLAB (@@@ Typed Data Access NOT working on RUSSEL)

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

   // Free temporarily allocated arrays
   free(iat_I);
   free(ja_I);
   free(coef_I);

}
