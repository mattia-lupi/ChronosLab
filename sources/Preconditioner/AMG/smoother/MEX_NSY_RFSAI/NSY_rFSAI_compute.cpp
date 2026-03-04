//----------------------------------------------------------------------------------------

#if defined PRINT
   bool dump = true;
#else
   bool dump = false;
#endif

#include "mex.h"
#include "Compute_nsy_rfsai.h"

//----------------------------------------------------------------------------------------
//
// MATLAB function
//
// function [iat_FL,ja_FL,coef_FL,iat_FU,ja_FU,coef_FU] = ...
//          NSY_rFSAI_compute(nstep,step_size,epsilon,nn_A,iat_A,ja_A,coef_A)
//
//----------------------------------------------------------------------------------------

// Gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

   // Get the value of the input scalars
   if (dump) mexPrintf("*** NSY_rFSAI_compute ***\n");
   if (dump) mexPrintf("- Get input scalars\n");
   int nstep      = mxGetScalar(prhs[0]);
   int step_size  = mxGetScalar(prhs[1]);
   double epsilon = mxGetScalar(prhs[2]);
   int nn_A       = mxGetScalar(prhs[3]);

   // Get the value of the input arrays
   if (dump) mexPrintf("- Get input arrays\n");
   int *iat_A      = (int*) mxGetData(prhs[4]);
   int *ja_A       = (int*) mxGetData(prhs[5]);
   double *coef_A  = (double*) mxGetData(prhs[6]);

   // Compute FL and FU entries
   if (dump) mexPrintf("- Compute FL and FU entries\n");
   int *iat_FL, *ja_FL;
   int *iat_FU, *ja_FU;
   double *coef_FL, *coef_FU;
   int ierr = Compute_nsy_rfsai(nstep,step_size,epsilon,nn_A,iat_A,ja_A,coef_A,
                                iat_FL,ja_FL,coef_FL,iat_FU,ja_FU,coef_FU);
   if (ierr != 0){
      cout << "ERROR in cpt_nsy_sfsai_coef with code: " << ierr << endl;
      return;
   }

   // Resize output arrays
   if (dump) mexPrintf("- Store FL and FU into the output arrays\n");

   int *pt_i;
   double *pt_r;

   int nt_FL = iat_FL[nn_A];
   // iat_FL
   mxArray *iat_FL_out  = mxCreateNumericMatrix((mwSize) 1,(mwSize) (nn_A+1),
                                                mxINT32_CLASS,mxREAL);
   pt_i = (int*) mxGetData(iat_FL_out);
   for (int k = 0; k < nn_A+1; k++) pt_i[k] = iat_FL[k]+1;
   free(iat_FL);
   // ja_FL
   mxArray *ja_FL_out   = mxCreateNumericMatrix((mwSize) 1,(mwSize) nt_FL,
                                                mxINT32_CLASS,mxREAL);
   pt_i = (int*) mxGetData(ja_FL_out);
   for (int k = 0; k < nt_FL; k++) pt_i[k] = ja_FL[k]+1;
   free(ja_FL);
   // coef_FL
   mxArray *coef_FL_out = mxCreateNumericMatrix((mwSize) 1,(mwSize) nt_FL,
                                                mxDOUBLE_CLASS,mxREAL);
   pt_r = (double*) mxGetData(coef_FL_out);
   for (int k = 0; k < nt_FL; k++) pt_r[k] = coef_FL[k];
   free(coef_FL);

   int nt_FU = iat_FU[nn_A];
   // iat_FU
   mxArray *iat_FU_out  = mxCreateNumericMatrix((mwSize) 1,(mwSize) (nn_A+1),
                                                mxINT32_CLASS,mxREAL);
   pt_i = (int*) mxGetData(iat_FU_out);
   for (int k = 0; k < nn_A+1; k++) pt_i[k] = iat_FU[k]+1;
   free(iat_FU);
   // ja_FU
   mxArray *ja_FU_out   = mxCreateNumericMatrix((mwSize) 1,(mwSize) nt_FU,
                                                mxINT32_CLASS,mxREAL);
   pt_i = (int*) mxGetData(ja_FU_out);
   for (int k = 0; k < nt_FU; k++) pt_i[k] = ja_FU[k]+1;
   free(ja_FU);
   // coef_FU
   mxArray *coef_FU_out = mxCreateNumericMatrix((mwSize) 1,(mwSize) nt_FU,
                                                mxDOUBLE_CLASS,mxREAL);
   pt_r = (double*) mxGetData(coef_FU_out);
   for (int k = 0; k < nt_FU; k++) pt_r[k] = coef_FU[k];
   free(coef_FU);

   // Return pointers
   plhs[0] = iat_FL_out;
   plhs[1] = ja_FL_out;
   plhs[2] = coef_FL_out;
   plhs[3] = iat_FU_out;
   plhs[4] = ja_FU_out;
   plhs[5] = coef_FU_out;

   if (dump) mexPrintf("\n");

}

//----------------------------------------------------------------------------------------
