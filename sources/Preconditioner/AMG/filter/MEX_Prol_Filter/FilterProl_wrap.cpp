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
#include <string.h>
////////////////////////////////
//#include <iostream>
//#include <fstream>
//using namespace std;
////////////////////////////////

#include "mex.h"
#include "FilterProl.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

// Prototype of MATLAB calling
//[nt_PF,iat_PF,ja_PF,coef_PF] = FilterProl_wrap(np,perc,tol,...
//                                               nn_P,iat_P,ja_P,coef_P,ntv,TV);

   // Check for proper number of arguments (Inputs and Outputs)
   if (nrhs != 10) {
       mexErrMsgIdAndTxt("ERROR in FilterComp_wrap","(10 inputs required.");
   }
   if(nlhs!=4) {
       mexErrMsgIdAndTxt("ERROR in FilterComp_wrap","4 outputs required.");
   }

   // Translate inputs
   int    np         =           mxGetScalar(prhs[ 0]);
   double perc       =           mxGetScalar(prhs[ 1]);
   double tol        =           mxGetScalar(prhs[ 2]);
   int    nn_P       =           mxGetScalar(prhs[ 3]);
   int    *iat_P     = (int*)    mxGetPr(    prhs[ 4]);
   int    *ja_P      = (int*)    mxGetPr(    prhs[ 5]);
   double *coef_P    = (double*) mxGetPr(    prhs[ 6]);
   int    nr_TV      =           mxGetScalar(prhs[ 7]);
   int    ntv        =           mxGetScalar(prhs[ 8]);
   double *TV        = (double*) mxGetPr(    prhs[ 9]);

   // @@@@@FORSE SI PUO METTERE TUTTO DIRETTAMENTE IN TV SENZA SPRECARE MEMORIA
   // Load TV in a 2D buffer
   double **TV_2D = (double**) malloc(nr_TV*sizeof(double*));
   if (TV_2D == NULL) mexErrMsgIdAndTxt("ERROR in FilterProl_wrap","TV_2D allocation");
   double *buffer = (double*) malloc(nr_TV*ntv*sizeof(double));
   if (buffer == NULL) mexErrMsgIdAndTxt("ERROR in FilterProl_wrap","buffer allocation");
   double *ptr = buffer;
   for (int i = 0; i < nr_TV; i++){
      TV_2D[i] = ptr;
      ptr += ntv;
   }
   memcpy(buffer,TV,nr_TV*ntv*sizeof(double));
   
   // Call the computational routine
   int nt_PF;
   int *iat_PF,*ja_PF;
   double *coef_PF;
   int ierr = FilterProl(np,perc,tol,nn_P,iat_P,ja_P,coef_P,ntv,TV_2D,
                         nt_PF,iat_PF,ja_PF,coef_PF);
 
   // Copy results in MATLAB (@@@ Typed Data Access NOT working on RUSSEL)
   //////////////////////////////////////////
   //cout << "COPIO RISULTATI FINALI"<< endl;
   //cout << "nn_PF:   " << nn_P << endl;
   //cout << "nt_PF:   " << nt_PF << endl;
   ////////////////////////////////////////

   // nt_PF
   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
   double *out_nt_PF = mxGetPr(plhs[0]);
   *out_nt_PF = (double) nt_PF;
   ////////////////////////////////////////
   //cout << "COPIATO nt_I" << endl;
   ////////////////////////////////////////

   // iat_PF
   plhs[1] = mxCreateDoubleMatrix((mwSize) nn_P+1,1,mxREAL);
   double *out_iat_PF = mxGetPr(plhs[1]);
   for (int i = 0; i < nn_P+1; i++) out_iat_PF[i] = (double) (iat_PF[i]+1);
   ////////////////////////////////////////
   //cout << "COPIATO iat_I" << endl;
   ////////////////////////////////////////

   // ja_PF
   plhs[2] = mxCreateDoubleMatrix((mwSize) nt_PF,1,mxREAL);
   double *out_ja_PF = mxGetPr(plhs[2]);
   for (int i = 0; i < nt_PF; i++) out_ja_PF[i] = (double) (ja_PF[i]+1);
   ////////////////////////////////////////
   //cout << "COPIATO ja_I" << endl;
   ////////////////////////////////////////

   // coef_PF
   plhs[3] = mxCreateDoubleMatrix((mwSize) nt_PF,1,mxREAL);
   double *out_coef_PF = mxGetPr(plhs[3]);
   for (int i = 0; i < nt_PF; i++) out_coef_PF[i] = coef_PF[i];

   // Free temporarily allocated arrays
   free(TV_2D);
   free(buffer);
   free(iat_PF);
   free(ja_PF);
   free(coef_PF);

}
