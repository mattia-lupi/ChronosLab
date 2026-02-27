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
#include "FilterComp.h"

/* The computational routine wrapper*/
int FilterComp_wrap(const int nthreads, const double tau, const int nn_A,
                    int *iat_A, int *ja_A, double *coef_A, const int nt_patt,
                    const int *iat_patt, const int *ja_patt, const int ntv,
                    const double *const *TV, int &nt_AC, int *&iat_AC,
                    int *&ja_AC, double *&coef_AC){

   int ierr = FilterComp(nthreads,tau,nn_A,iat_A,ja_A,coef_A,nt_patt,iat_patt,ja_patt,
                         ntv,TV,nt_AC,iat_AC,ja_AC,coef_AC);

   return ierr;
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

// Prototype of MATLAB calling
// [nt_AC,iat_AC,ja_AC,coef_AC] = FilterComp_wrap(np,tau,nn_A,iat_A,ja_A,coef_A,...
//                                                nt_patt,iat_patt,ja_patt,ntv,TV);

   // Check for proper number of arguments (Inputs and Outputs)
   if (nrhs != 11) {
       mexErrMsgIdAndTxt("ERROR in FilterComp_wrap","11 inputs required.");
   }
   if(nlhs!=4) {
       mexErrMsgIdAndTxt("ERROR in FilterComp_wrap","4 outputs required.");
   }

   // Translate inputs
   int    np         =           mxGetScalar(prhs[ 0]);
   double tau        =           mxGetScalar(prhs[ 1]);
   int    nn_A       =           mxGetScalar(prhs[ 2]);
   int    *iat_A     = (int*)    mxGetPr(    prhs[ 3]);
   int    *ja_A      = (int*)    mxGetPr(    prhs[ 4]);
   double *coef_A    = (double*) mxGetPr(    prhs[ 5]);
   int    nt_patt    =           mxGetScalar(prhs[ 6]);
   int    *iat_patt  = (int*)    mxGetPr(    prhs[ 7]);
   int    *ja_patt   = (int*)    mxGetPr(    prhs[ 8]);
   int    ntv        =           mxGetScalar(prhs[ 9]);
   double *TV        = (double*) mxGetPr(    prhs[10]);

   // @@@@@FORSE SI PUO METTERE TUTTO DIRETTAMENTE IN TV SENZA SPRECARE MEMORIA
   // Load TV in a 2D buffer
   double **TV_2D = (double**) malloc(nn_A*sizeof(double*));
   if (TV_2D == NULL) mexErrMsgIdAndTxt("ERROR in FilterComp_wrap","TV_2D allocation");
   double *buffer = (double*) malloc(nn_A*ntv*sizeof(double));
   if (buffer == NULL) mexErrMsgIdAndTxt("ERROR in FilterComp_wrap","buffer allocation");
   double *ptr = buffer;
   for (int i = 0; i < nn_A; i++){
      TV_2D[i] = ptr;
      ptr += ntv;
   }
   memcpy(buffer,TV,nn_A*ntv*sizeof(double));

   // Call the computational routine
   int nt_AC;
   int *iat_AC,*ja_AC;
   double *coef_AC;
   ////////////////////////////////////////
   //cout << "ENTRO"<< endl;
   ////////////////////////////////////////
   int ierr = FilterComp_wrap(np,tau,nn_A,iat_A,ja_A,coef_A,nt_patt,iat_patt,ja_patt,
                              ntv,TV_2D,nt_AC,iat_AC,ja_AC,coef_AC);
 
   // Copy results in MATLAB (@@@ Typed Data Access NOT working on RUSSEL)
   ////////////////////////////////////////
   //cout << "COPIO RISULTATI FINALI"<< endl;
   //cout << "nn_I:   " << nn_I << endl;
   //cout << "nt_I:   " << nt_I << endl;
   ////////////////////////////////////////

   // nt_AC
   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
   double *out_nt_AC = mxGetPr(plhs[0]);
   *out_nt_AC = (double) nt_AC;
   ////////////////////////////////////////
   //cout << "COPIATO nt_I" << endl;
   ////////////////////////////////////////

   // iat_AC
   plhs[1] = mxCreateDoubleMatrix((mwSize) nn_A+1,1,mxREAL);
   double *out_iat_AC = mxGetPr(plhs[1]);
   for (int i = 0; i < nn_A+1; i++) out_iat_AC[i] = (double) (iat_AC[i]+1);
   ////////////////////////////////////////
   //cout << "COPIATO iat_I" << endl;
   ////////////////////////////////////////

   // ja_AC
   plhs[2] = mxCreateDoubleMatrix((mwSize) nt_AC,1,mxREAL);
   double *out_ja_AC = mxGetPr(plhs[2]);
   for (int i = 0; i < nt_AC; i++) out_ja_AC[i] = (double) (ja_AC[i]+1);
   ////////////////////////////////////////
   //cout << "COPIATO ja_I" << endl;
   ////////////////////////////////////////

   // coef_AC
   plhs[3] = mxCreateDoubleMatrix((mwSize) nt_AC,1,mxREAL);
   double *out_coef_AC = mxGetPr(plhs[3]);
   for (int i = 0; i < nt_AC; i++) out_coef_AC[i] = coef_AC[i];

   // Free temporarily allocated arrays
   free(TV_2D);
   free(buffer);
   free(iat_AC);
   free(ja_AC);
   free(coef_AC);

}
