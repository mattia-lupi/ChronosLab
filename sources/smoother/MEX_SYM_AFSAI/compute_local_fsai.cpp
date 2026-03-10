
#include "compute_local_fsai.h" 

//----------------------------------------------------------------------------------------

void compute_local_fsai(iReg nthread, iReg n_step, iReg step_size, rExt tau, rExt eps,
                        iReg nrows, iReg nrows_M, iExt nterm_M, const iExt * const iat_M,
                        const iReg * const ja_M, const rExt * const coef_M, iExt *nterm_G,
                        iExt *iat_G, iReg *ja_G, rExt *coef_G){

// Set the chunk size for the parallel do (this way each thread should execute 20 loops)
iReg chunk_size = nrows / (20*nthread); chunk_size = max(1,chunk_size);

// Allocate scratches
iExt kmax    = 1 + (iExt) n_step * (iExt) step_size;
iExt nzmax_G = (iExt) nrows * kmax;

vector<iExt> nt_slice; nt_slice.resize(nthread);
vector<iExt> istart_scr; istart_scr.resize(nrows);
vector<iExt> istop_scr; istop_scr.resize(nrows);
vector<iReg> ja_scr; ja_scr.resize(nzmax_G);
vector<rExt> coef_scr; coef_scr.resize(nzmax_G);

// Set pointers to the beginning of each row
iExt ind = 1;
for ( iReg irow = 0; irow < nrows; irow++ ){
   istart_scr[irow] = ind;
   ind += kmax;
}

// Parallel computation of G
nterm_G[0] = 0;
iReg shift = nrows_M - nrows;

#pragma omp parallel num_threads(nthread)
{

   iReg myid = omp_get_thread_num();

   iExt loc_nt_G = 0;
   double time_1 = omp_get_wtime();
   cpt_afsai_coef(chunk_size,n_step,step_size,tau,eps,shift,nrows,nrows_M,
                  nterm_M,loc_nt_G,iat_M,ja_M,istart_scr.data(),istop_scr.data(),
                  ja_scr.data(),coef_M,coef_scr.data());
   double time_2 = omp_get_wtime();

   // Set the local number of non-zeroes
   #pragma omp atomic
   nterm_G[0] += loc_nt_G;

// if (myid == 0) {cout << time_2-time_1 << endl;}

} // end parallel region

//////////////////////////////////////////////////////////////////////////////////////////
//cout << "istart - istop" << endl;
//for (iExt i = 0; i < nrows; i++){
//   cout << istart_scr[i] << "    " << istop_scr[i] << endl;
//}
//////////////////////////////////////////////////////////////////////////////////////////

#pragma omp parallel num_threads(nthread)
{

   // Retrieve processor ID
   iReg myid = omp_get_thread_num();

   // Get row indices
   iReg my_nrows    = nrows / nthread;
   iReg my_firstrow = myid * my_nrows + 1;
   iReg my_lastrow  = my_firstrow + my_nrows - 1;
   if ( myid + 1 == nthread ){
      iReg resto = nrows%nthread;
      my_lastrow += resto;
      my_nrows   += resto;
   }
   iExt kcount = my_nrows;
   for ( iReg irow = my_firstrow-1; irow < my_lastrow; irow++ ) {
      kcount += istop_scr[irow] - istart_scr[irow];
   }
   nt_slice[myid] = kcount;
   // Sync threads
   #pragma omp barrier

   // Set pointer 
   iExt ind_G  = 1;
   for ( iReg id = 0; id < myid; id++ ) {
      ind_G += nt_slice[id];
   }
   for ( iReg irow = my_firstrow-1; irow < my_lastrow; irow++ ) {

      iat_G[irow] = ind_G;
      iExt istart = istart_scr[irow];
      iExt istop  = istop_scr[irow];
      iExt nadd   = istop-istart;

      iExt j = ind_G-1;
      for ( iExt i = istart-1; i < istop; i++ ) {
         ja_G[j]   = ja_scr[i];
         coef_G[j] = coef_scr[i];
         j += 1; 
      }

      ind_G += nadd+1;
   }

   if (myid + 1 == nthread){ iat_G[my_lastrow] = ind_G;}
  
} // end parallel region

}

