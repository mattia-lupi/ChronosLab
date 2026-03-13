#include "mex.h"

//----------------------------------------------------------------------------------------
//
// MATLAB function
//
// function DumpChronos(nrows,nterm,ntv,irow,jcol,coef,TV0);
//
//----------------------------------------------------------------------------------------

// Gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

   // Get the value of the input scalars
   mexPrintf("*** DUMPING CHRONOS BINARY DATA ***\n");
   long int nrows   = mxGetScalar(prhs[0]);
   long int nterm   = mxGetScalar(prhs[1]);
   long int ntv     = mxGetScalar(prhs[2]);
   
   // Get the value of the input arrays
   long int *irow    = (long int*)    mxGetData(prhs[3]);
   long int *jcol    = (long int*)    mxGetData(prhs[4]);
   double *coef      = (double*)      mxGetData(prhs[5]);
   double *TV0       = (double*)      mxGetData(prhs[6]);

   FILE *of_matrix = fopen("ChronosMatrix.Ext_bin","wb");
   fwrite(&nrows,sizeof(long int),1,of_matrix);
   fwrite(&nterm,sizeof(long int),1,of_matrix);
   for (long int i = 0; i < nterm; i++){
      fwrite(irow+i,sizeof(long int),1,of_matrix);
      fwrite(jcol+i,sizeof(long int),1,of_matrix);
      fwrite(coef+i,sizeof(double),1,of_matrix);
   }
   fclose(of_matrix);

   FILE *of_tv0 = fopen("ChronosTV0.Ext_bin","wb");
   fwrite(&nrows,sizeof(long int),1,of_tv0);
   fwrite(&ntv,sizeof(long int),1,of_tv0);
   for (long int i = 0; i < nrows; i++){
      for (long int j = 0; j < ntv; j++){
         fwrite(TV0+j*nrows+i,sizeof(double),1,of_matrix);
      }
   }
   fclose(of_matrix);

   mexPrintf("*** DUMPING CHRONOS BINARY DATA COMPLETED ***\n");

}
