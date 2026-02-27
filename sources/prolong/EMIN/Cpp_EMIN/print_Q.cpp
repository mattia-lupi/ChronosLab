#include <stdio.h>
////////////////////////////////////////////////////
#include <iostream>
using namespace std;
////////////////////////////////////////////////////

void print_Q(FILE *ofile, const int nblk, const int ntv, const int *pt_blk,
             const double *mat_Q){

    int ind_Q = 0;
    int irow = 0;
    int icol = 0;
    int iend = pt_blk[0];
    for (int iblk = 0; iblk < nblk; iblk++){
       int istart = iend;
       iend = pt_blk[iblk+1];
       int nr = iend - istart;
       ///////////////////////////////////////////
       //cout << "IBLK " << iblk << " ISTART " << istart << " IEND " << iend << endl;
       //cout << "IROW " << irow << " ICOL " << icol << endl;
       //cout.flush();
       ///////////////////////////////////////////
       int k = ind_Q;
       for (int i = istart; i < iend; i++){
          for (int j = 0; j < ntv; j++){
             fprintf(ofile,"%10d %10d %25.15e\n",irow+1,icol+j+1,mat_Q[k+j*nr]);
          }
          k++;
          irow++;
       }
       ind_Q += nr*ntv;
       if (nr > 0) icol += ntv;
    }

}
