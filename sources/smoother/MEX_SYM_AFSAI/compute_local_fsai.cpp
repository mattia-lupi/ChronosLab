
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstdint>   // int64_t, int32_t
#include <string>
#include <cassert>
#include <vector>     // to use: vector
#include <cmath>      // to use: abs,sqrt
#include <algorithm>  // to use: min
#include <lapacke.h>  // to use: dpotrf,dpotrs
#include <omp.h>      // omp_get_wtime
#include <iostream>
#include <iomanip>
using namespace std;

//----------------------------------------------------------------------------------------

// Regular integer
typedef int iReg;

// Extended integer
typedef long int iExt;

// Regular real
typedef float rReg;

// Extended real
typedef double rExt;

// Blas integer
typedef int iBlas;

// Lapack integer
typedef int iLapack;

//----------------------------------------------------------------------------------------

rExt cpt_ddot(iReg n, rExt *x1, rExt *x2){
   rExt ddot = 0.;
   for ( iReg i = 0; i < n; i++ ){
      ddot += x1[i] * x2[i];
   }
   return ddot;
}

//----------------------------------------------------------------------------------------

rExt cpt_dnrm2(iReg n, rExt *x1){
   rExt dnrm2 = cpt_ddot(n,x1,x1);
   dnrm2 = sqrt(dnrm2);
   return dnrm2;
}

//----------------------------------------------------------------------------------------

void copy_ja(iExt n, iReg *x1, iReg *x2){
   for ( iExt i = 0; i < n; i++ ){
      x2[i] = x1[i];
   }
}

//----------------------------------------------------------------------------------------

void copy_coef(iExt n, rExt *x1, rExt *x2){
   for ( iExt i = 0; i < n; i++ ){
      x2[i] = x1[i];
   }
}

//----------------------------------------------------------------------------------------

void ri_sortsplit(iReg n, iReg ncut, rExt * const R_vec, iReg * const I_vec){

iReg first = 1;
iReg last = n;

// Check ncut consistency
if (ncut < first || ncut > last) return;

while(true){

   iReg mid = first;
   rExt absval = abs( R_vec[mid-1] );

   for (iReg j = first+1; j < last+1; j ++){
      //cout << showpos;
      //cout << scientific;
      //cout << setprecision(17);
      //cout << " " << j << " " << abs(R_vec[j-1]) << " " << absval << endl;
      if (abs(R_vec[j-1]) > absval){
         mid += 1;
         // Exchange
         rExt tmp  = R_vec[mid-1];
         iReg itmp = I_vec[mid-1];
         R_vec[mid-1] = R_vec[j-1];
         I_vec[mid-1] = I_vec[j-1];
         R_vec[j-1]   = tmp;
         I_vec[j-1]   = itmp;
      }
   }

   // Exchange
   rExt tmp        = R_vec[mid-1];
   R_vec[mid-1]    = R_vec[first-1];
   R_vec[first-1]  = tmp;

   iReg itmp      = I_vec[mid-1];
   I_vec[mid-1]   = I_vec[first-1];
   I_vec[first-1] = itmp;

   // Exit test
   if (mid == ncut) return;
   if (mid >  ncut){
      last = mid-1;
   }else{
      first = mid+1;
   }

}

}

//----------------------------------------------------------------------------------------

void swapi(iReg &i1, iReg &i2){
iReg tmp = i1;
i1 = i2;
i2 = tmp;
}

//----------------------------------------------------------------------------------------

void iheapsort(iReg * const x1, iReg n){

for (iReg node = 2; node < n+1; node ++){
   iReg i = node;
   iReg j = i/2;
   while( x1[j-1] < x1[i-1] ){
      swapi(x1[j-1],x1[i-1]);
      i = j;
      j = i/2;
      if (i == 1) break;
   }
}

for (iReg i = n; i > 1; i --){
   swapi(x1[i-1],x1[0]);
   iReg k = i - 1;
   iReg ik = 1;
   iReg jk = 2;
   if (k >= 3){
      if (x1[2] > x1[1]) jk = 3;
   }
   bool cont_cycle = false;
   if (jk <= k){
      if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
   }
   while (cont_cycle){
      swapi(x1[jk-1],x1[ik-1]);
      ik = jk;
      jk = ik*2;
      if (jk+1 <= k){
         if (x1[jk] > x1[jk-1]) jk = jk+1;
      }
      cont_cycle = false;
      if (jk <= k){
         if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
      }
   }
}

}

//----------------------------------------------------------------------------------------

void gather_fullsys(bool &nulrhs, iReg irow, iReg mrow, iReg jendbloc, iReg nequ,
                    iExt nterm, iReg mmax, const iExt * const iat, const iReg * const ja,
                    const iReg * const vecinc, const rExt * const mat_A,
                    rExt * const full_A, rExt * const rhs){

// Load all the rows of full_A from the first to the last
for (iReg i = 1; i < mrow+1; i++){

   // Load the i-th column exploring the row-th row of mat_A
   iReg ii     = i;
   iReg row    = vecinc[i-1];
   iExt jj     = iat[row-1];
   iExt endrow = iat[row];

   while (ii <= mrow){

      // Make sure that ja(jj-1) >= vecinc(ii-1) (only upper part)
      while (ja[jj-1] < vecinc[ii-1]){
         jj += 1;
         if (jj == endrow) {
            // If the row end is reached move to the next column
            for (iReg k = ii; k < mrow+1; k++){
               full_A[(i-1)*mmax+k-1] = 0.;
            }
            goto exit_column_loop;
         }
      }

      if (vecinc[ii-1] == ja[jj-1]){
         // If ja(jj-1) == vecinc(ii-1), load the term mat_A(jj-1)
         full_A[(i-1)*mmax+ii-1] = mat_A[jj-1];
         ii += 1;
      }else{
         // If ja(jj-1) != vecinc(ii-1), the term is missing and set a 0
         full_A[(i-1)*mmax+ii-1] = 0.;
         ii += 1;
      }

   } // end column loop
   exit_column_loop: ;

} // end row loop


// Load the right-hand-side exploring the row/column irow of mat_A
// The sign of the right-hand-side is changed
iReg ii = 1;
iExt jj = iat[irow-1];
nulrhs = ja[jj-1] >= jendbloc;
if (nulrhs) return;
while (ja[jj-1] < jendbloc){
   if (vecinc[ii-1] > ja[jj-1]){
      jj += 1;
   }else if (vecinc[ii-1] == ja[jj-1]){
      rhs[ii-1] = - mat_A[jj-1];
      jj += 1;
      ii += 1;
      if (ii > mrow) return;
   }else{
      rhs[ii-1] = 0.;
      ii += 1;
      if (ii > mrow) return;
   }
}
for (iReg k = ii; k < mrow+1; k++){
   rhs[k-1] = 0.;
}

}

//----------------------------------------------------------------------------------------

void kap_grad(iReg irow, iReg nequ, iExt nterm, iReg mmax, iReg &mrow, iReg jendbloc,
              iReg lfil, const iExt * const iat, const iReg * const ja,
              const rExt * const mat_A, const rExt * const rhs, iReg * const IWN,
              iReg * const JWN, rExt * const WR){

iReg ind_WR = 0;

// Get the degree zero entries (A pattern)
iExt j = iat[irow-1];
iReg jjcol = ja[j-1];
while (jjcol < jendbloc){
   // If JWN is positive add this entry to the tentative pattern
   if (JWN[jjcol-1] >= 0) {
      ind_WR += 1;
      JWN[jjcol-1] = ind_WR;
      IWN[mrow+ind_WR-1] = ja[j-1];
      WR[ind_WR-1] = mat_A[j-1];
   }
   j += 1;
   jjcol = ja[j-1];
}

/////////////////////////////////////////////////////////////////////
//DEBUG PRINT
//cout << "WR-1" << endl;
//for (iReg i = 0; i < mrow; i ++){
//   cout << i+1 << " " << WR[i] << endl;
//}
/////////////////////////////////////////////////////////////////////

// Get the higher degree entries
for (iReg i = 1; i < mrow+1; i ++){
   iReg iirow = IWN[i-1];
   for (iExt j = iat[iirow-1]; j < iat[iirow]; j ++){
      jjcol = ja[j-1];
      if (jjcol < jendbloc){
         iReg ind = JWN[jjcol-1];
         if (ind == 0){
            // New fill-in entry
            ind_WR += 1;
            JWN[jjcol-1] = ind_WR;
            IWN[mrow+ind_WR-1] = ja[j-1];
            WR[ind_WR-1] = rhs[i-1]*mat_A[j-1];
	 }else if (ind > 0){
            // Already found fill-in entry
            WR[ind-1] += rhs[i-1]*mat_A[j-1];
         }
      }
   }
}

/////////////////////////////////////////////////////////////////////
//DEBUG PRINT
//cout << "WR-2" << endl;
//for (iReg i = 0; i < mrow; i ++){
//   cout << i+1 << " " << WR[i] << endl;
//}
/////////////////////////////////////////////////////////////////////

// Keep only the largest lfil entries
iReg ncut = min(lfil,ind_WR);
ri_sortsplit(ind_WR,ncut,WR,&IWN[mrow]);

/////////////////////////////////////////////////////////////////////
//DEBUG PRINT
//cout << "WR-3" << endl;
//for (iReg i = 0; i < mrow; i ++){
//   cout << i+1 << " " << WR[i] << endl;
//}
/////////////////////////////////////////////////////////////////////

// Mark retained JWN with a -1
for (iReg i = mrow+1; i < mrow+ncut+1; i ++){
   JWN[IWN[i-1]-1] = -1;
}

// Reset the other non-zero indicators
for (iReg i = mrow+ncut+1; i < mrow+ind_WR+1; i ++){
   JWN[IWN[i-1]-1] = 0;
}

// Update IWN dimension and sort it
mrow += ncut;
iheapsort(IWN,mrow);

}

//----------------------------------------------------------------------------------------

void cpt_afsai_coef(iReg chunk_size, iReg n_step, iReg step_size, rExt tau, rExt eps,
                    iReg shift, iReg nrows, iReg nequ, iExt nterm, iExt &nterm_G,
                    const iExt * const iat, const iReg * const ja, iExt * const istart_G,
                    iExt * const istop_G, iReg * const ja_G, const rExt * const coef_A,
                    rExt * const coef_G){

// Set parameters
iReg mrow_min = 5;

// warning message for the use of BLAS and LAPACK
if ( sizeof(iReg) != sizeof(iBlas) || sizeof(iReg) != sizeof(iLapack) ){
   throw "WARNING: check integer kind compatibility with BLAS and LAPACK";
}

// Compute maximum number of steps
iReg mmax = n_step * step_size;

// Allocate local work arrays
vector<iReg> IWN(nequ,0);
vector<iReg> JWN(nequ,0);
vector<rExt> WR(nequ,0);
vector<rExt> full_A(mmax*mmax,0);
vector<rExt> rhs(mmax+1,0);
vector<rExt> rhs_sav(mmax+1,0);

// Initialize number of entries
nterm_G = 0;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma omp for nowait schedule(dynamic,chunk_size)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Loop over the rows of the current processor
for( iReg irow = 1; irow < nrows+1; irow++){
//for( iReg irow = 1; irow < 59+1; irow++){

   iReg mrow = 0;
   iReg irow_glo = irow + shift;

   // Loop for the refinement of the row pattern
   iReg istep = 0;
   rExt DKap_old = 0.;
   bool Refine = n_step >= 1;
   while (Refine){

      istep += 1;

      // Apply the dropping strategy only with at least mrow_min terms
      if ((tau > 0.) && (mrow > mrow_min)){

         rhs[mrow] = 1.;
         rExt asstol = tau * cpt_dnrm2(mrow+1,rhs.data());
         // Drop the row
         iReg ind = 0;
         for (iReg i = 0; i < mrow; i ++){
            if (abs(rhs[i]) > asstol) {
               // Retain the term and add it to IWN
               rhs[ind] = rhs[i];
               rhs_sav[ind] = rhs_sav[i];
               IWN[ind] = IWN[i];
               ind += 1;
	    }else{
               // Drop this term and eliminate from JWN
               JWN[IWN[i]-1] = -1;
            }
         }
         // Update mrow
         mrow = ind;
      }
 
      /////////////////////////////////////////////////////////////////////
      //DEBUG PRINT
 /*   if( (irow == 59 && istep == 16) ) {
         cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
         cout << irow  << endl;
         cout << istep << endl;
         cout << "IWN" << endl;
         for (iReg i = 0; i < mrow; i ++){
            cout << i+1 << " " << IWN[i] << endl;
         }
         cout << "JWN" << endl;
         for (iReg i = 0; i < mrow; i ++){
            cout << i+1 << " " << JWN[i] << endl;
         }
         cout << "WR" << endl;
         for (iReg i = 0; i < mrow; i ++){
            cout << i+1 << " " << WR[i] << endl;
         }
      }
  */  /////////////////////////////////////////////////////////////////////

      // Compute the Kaporin gradient
      iReg mrow_old = mrow;
      kap_grad(irow_glo,nequ,nterm,mmax,mrow,irow_glo,step_size,iat,ja,coef_A,rhs.data(),
               IWN.data(),JWN.data(),WR.data());

      /////////////////////////////////////////////////////////////////////
      //DEBUG PRINT
  /*  if( (irow == 59 && istep == 16) ) {
         cout << "**************************************************" << endl;
         cout << irow  << endl;
         cout << istep << endl;
         cout << "IWN" << endl;
         for (iReg i = 0; i < mrow; i ++){
            cout << i+1 << " " << IWN[i] << endl;
         }
         cout << "JWN" << endl;
         for (iReg i = 0; i < mrow; i ++){
            cout << i+1 << " " << JWN[i] << endl;
         }
         cout << "WR" << endl;
         for (iReg i = 0; i < mrow; i ++){
            cout << i+1 << " " << WR[i] << endl;
         }
         break;
      }
  */  /////////////////////////////////////////////////////////////////////

      // Compute the G row if the pattern is not null
      if (mrow > mrow_old){

         // Gather the system coefficients
         bool nulrhs;
         gather_fullsys(nulrhs,irow_glo,mrow,irow_glo,nequ,nterm,mmax,iat,ja,IWN.data(),
                        coef_A,full_A.data(),rhs.data());

         if (nulrhs == false){

            // If the rhs is not null
            lapack_int info;
            //info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',mrow,full_A.data(),mmax);
            info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',mrow,full_A.data(),mmax);

	    // Backward and forward substitution
	    rhs_sav = rhs;
            //info = LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'U',mrow,1,full_A.data(),mmax,rhs.data(),1);
            info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',mrow,1,full_A.data(),mmax,rhs.data(),mrow);

         }

         // Compute the Kaporin number decrease
         rExt DKap_new = cpt_ddot(mrow,rhs.data(),rhs_sav.data());

         // Exit check
         if (istep == n_step){
            Refine = false;
         }else{
            Refine = (abs(DKap_new-DKap_old) >= eps*DKap_old) && (DKap_new != 0.);
            DKap_old = DKap_new;
         }

      }else{

         // If the pattern is empty the row is uncoupled
         Refine = false;

      }

   } // end refinement loop

   // Compute the scaling factor for this row
   iExt ind = iat[irow_glo-1];
   while (ja[ind-1] < irow_glo){
      ind += 1;
   }

   rExt scal_fac = coef_A[ind-1] - cpt_ddot(mrow,rhs_sav.data(),rhs.data());

   if (scal_fac < 0.){

      // Recompute this row in safer way
      bool nulrhs;
      gather_fullsys(nulrhs,irow_glo,mrow,irow_glo,nequ,nterm,mmax,iat,ja,IWN.data(),
                     coef_A,full_A.data(),rhs.data());
      lapack_int info;
      //info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',mrow,full_A.data(),mmax);
      info = LAPACKE_dpotrf(LAPACK_COL_MAJOR,'L',mrow,full_A.data(),mmax);
      rhs_sav = rhs;
      //info = LAPACKE_dpotrs(LAPACK_ROW_MAJOR,'U',mrow,1,full_A.data(),mmax,rhs.data(),1);
      info = LAPACKE_dpotrs(LAPACK_COL_MAJOR,'L',mrow,1,full_A.data(),mmax,rhs.data(),mrow);
      scal_fac = coef_A[ind-1] - cpt_ddot(mrow,rhs_sav.data(),rhs.data());
   }

   scal_fac = 1. / sqrt(scal_fac);

   // Get the beginning of G
   iExt ind_G  = istart_G[irow-1];
   iExt ind_G0 = ind_G;

   // Compute the terms found in G
   for (iReg i = 1; i < mrow+1; i++){
      coef_G[ind_G-1] = scal_fac*rhs[i-1];
      ja_G[ind_G-1] = IWN[i-1];
      ind_G += 1;
      // Reset the non-zero indicator
      JWN[IWN[i-1]-1] = 0;
   }

   // Load the diagonal entry
   coef_G[ind_G-1] = scal_fac;
   ja_G[ind_G-1]   = irow_glo;

   // Set the pointer to the end of the row
   istop_G[irow-1] = ind_G;

   // Count the entries of G
   nterm_G += (ind_G - ind_G0) + 1;

} // end row loop

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// end omp for
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

}

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

//----------------------------------------------------------------------------------------
// compute_local_fsai_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [nterm_G, iat_G, ja_G, coef_G] =
//       compute_local_fsai_wrap(nthread, n_step, step_size, tau, eps,
//                               nrows, nrows_M, nterm_M,
//                               iat_M, ja_M, coef_M)
//
// NOTE: Add the appropriate header include for compute_local_fsai(),
//       copy_ja(), and copy_coef() — they are not declared in the
//       original source and must come from a project-specific header.
//
// Build command:
//   See compile.m — ensure -R2018a is on its own line
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexPrintf() not declared in the pure C++ MEX API (mex.hpp does
//         not pull in mex.h). Replaced with mprint() helper routing
//         through getEngine()->feval(u"fprintf", createCharArray(...)).
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         Input arrays copied into std::vector<T>; .data() passed to kernel.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() requires std::is_arithmetic<T>.
//         All string arguments use factory.createCharArray() instead.
//
// [NEW-1] mxGetData() with void* cast replaced by TypedArray<T> accessors:
//           iat_M  → TypedArray<int64_t>  (was mxINT64_CLASS / iExt*)
//           ja_M   → TypedArray<int32_t>  (was mxINT32_CLASS / iReg*)
//           coef_M → TypedArray<double>
//
// [NEW-2] Output/work arrays were created as mxArray* and written to by
//         the kernel via raw pointers obtained from mxGetData().
//         In the C++ API TypedArray elements are proxies — raw pointers
//         cannot be obtained directly. Solution: use std::vector<T> as
//         the kernel-writable buffer, then copy into TypedArray after.
//           - nterm_G : single int64_t variable (scalar output of kernel)
//           - iat_G   : std::vector<int64_t>, length nrows+1
//           - ja_G    : std::vector<int32_t>, length nzmax_G (oversized work)
//           - coef_G  : std::vector<double>,  length nzmax_G (oversized work)
//         The mxDestroyArray(ja_W) / mxDestroyArray(coef_W) pattern
//         (oversized alloc → copy to resized → destroy original) is
//         replaced naturally: vectors are resized implicitly via the
//         TypedArray constructor that takes an exact count.
//
// [NEW-3] iExt = long int maps to int64_t on 64-bit Linux/macOS (LP64).
//         std::vector<iExt> is used for iat_M, iat_G so the kernel
//         receives the exact pointer type it expects. A static_assert
//         guards against builds where long int is not 64 bits.
//
// [NEW-4] Removed #include <iostream> / <iomanip> / using namespace std.
//         All diagnostic output goes through mprint().
//
// [NEW-5] std::size_t used for all array sizes and loop bounds.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

// [NEW-3] Guard against LP32/ILP64 platforms where long int != 64 bits
static_assert(sizeof(iExt) == sizeof(int64_t),
              "iExt (long int) must be 64 bits on this platform. "
              "Update the typedef or the TypedArray<int64_t> mapping.");

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    // [FIX-A] Route diagnostic output through MATLAB's fprintf
    void mprint(const std::string& msg)
    {
        getEngine()->feval(u"fprintf", 0,
            std::vector<Array>{ factory.createCharArray(msg) });
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Read input scalars
        // -----------------------------------------------------------------------
        mprint("- get input scalars\n");

        const iReg nthread   = static_cast<iReg>(TypedArray<double>(inputs[0])[0]);
        const iReg n_step    = static_cast<iReg>(TypedArray<double>(inputs[1])[0]);
        const iReg step_size = static_cast<iReg>(TypedArray<double>(inputs[2])[0]);
        const rExt tau       = static_cast<rExt>(TypedArray<double>(inputs[3])[0]);
        const rExt eps       = static_cast<rExt>(TypedArray<double>(inputs[4])[0]);
        const iReg nrows     = static_cast<iReg>(TypedArray<double>(inputs[5])[0]);
        const iReg nrows_M   = static_cast<iReg>(TypedArray<double>(inputs[6])[0]);
        const iExt nterm_M   = static_cast<iExt>(TypedArray<double>(inputs[7])[0]);

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-1] Copy input arrays into std::vector.
        //   iat_M  → int64_t  (mxINT64_CLASS — was iExt* via mxGetData void* cast)
        //   ja_M   → int32_t  (mxINT32_CLASS — was iReg* via mxGetData void* cast)
        //   coef_M → double
        // -----------------------------------------------------------------------
        mprint("- get input arrays\n");

        const TypedArray<int64_t> iat_M_arr  = inputs[8];
        const TypedArray<int32_t> ja_M_arr   = inputs[9];
        const TypedArray<double>  coef_M_arr = inputs[10];

        // [NEW-3] Use std::vector<iExt> (= long int) so kernel receives
        //         the exact pointer type — safe on LP64 (Linux/macOS 64-bit)
        std::vector<iExt> iat_M_vec (iat_M_arr.begin(),  iat_M_arr.end());
        std::vector<iReg> ja_M_vec  (ja_M_arr.begin(),   ja_M_arr.end());
        std::vector<rExt> coef_M_vec(coef_M_arr.begin(), coef_M_arr.end());

        // -----------------------------------------------------------------------
        // [NEW-2] Allocate kernel-writable output / work buffers as std::vector.
        //
        //   Original pattern:
        //     plhs[0] = mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL)  ← scalar
        //     plhs[1] = mxCreateNumericMatrix(1,nrows+1,mxINT64_CLASS,...) ← iat_G
        //     ja_W    = mxCreateNumericMatrix(1,nzmax_G,mxINT32_CLASS,...) ← work
        //     coef_W  = mxCreateNumericMatrix(1,nzmax_G,mxDOUBLE_CLASS,...) ← work
        //
        //   Modern equivalent: plain variables / vectors the kernel writes into.
        //   After the kernel, exact-size TypedArrays are constructed for output.
        //   No mxDestroyArray() needed — vectors self-destruct.
        // -----------------------------------------------------------------------
        mprint("- allocate work arrays\n");

        const iExt kmax    = 1 + static_cast<iExt>(n_step) * static_cast<iExt>(step_size);
        const iExt nzmax_G = static_cast<iExt>(nrows) * kmax;

        iExt              nterm_G_val = 0;                      // scalar output
        std::vector<iExt> iat_G_vec(static_cast<std::size_t>(nrows) + 1, 0);
        std::vector<iReg> ja_G_vec (static_cast<std::size_t>(nzmax_G),   0);
        std::vector<rExt> coef_G_vec(static_cast<std::size_t>(nzmax_G),  0.0);

        // -----------------------------------------------------------------------
        // Call the computational kernel
        // -----------------------------------------------------------------------
        mprint("- compute G terms\n");

        compute_local_fsai(nthread, n_step, step_size, tau, eps,
                           nrows, nrows_M, nterm_M,
                           iat_M_vec.data(),
                           ja_M_vec.data(),
                           coef_M_vec.data(),
                           &nterm_G_val,
                           iat_G_vec.data(),
                           ja_G_vec.data(),
                           coef_G_vec.data());

        // -----------------------------------------------------------------------
        // [NEW-2] Pack results into exact-size MATLAB output arrays.
        //
        //   Original: oversized ja_W / coef_W copied to resized ja_R / coef_R
        //             then mxDestroyArray(ja_W) / mxDestroyArray(coef_W).
        //   Modern:   just copy the first nterm_G_val elements from the vectors
        //             directly into fresh TypedArrays of the exact size.
        //             The oversized vectors self-destruct at end of scope.
        // -----------------------------------------------------------------------
        mprint("- resize output arrays\n");

        const std::size_t sz_rows1 = static_cast<std::size_t>(nrows) + 1;
        const std::size_t sz_nt    = static_cast<std::size_t>(nterm_G_val);

        // --- nterm_G : int64 scalar -------------------------------------------
        TypedArray<int64_t> out_nterm_G =
            factory.createArray<int64_t>({1, 1});
        out_nterm_G[0] = static_cast<int64_t>(nterm_G_val);

        // --- iat_G : int64, length nrows+1 ------------------------------------
        TypedArray<int64_t> out_iat_G =
            factory.createArray<int64_t>({1, sz_rows1});
        {
            auto it = out_iat_G.begin();
            for (std::size_t i = 0; i < sz_rows1; ++i, ++it)
                *it = static_cast<int64_t>(iat_G_vec[i]);
        }

        // --- ja_G : int32, length nterm_G_val (trimmed from oversized work) ---
        TypedArray<int32_t> out_ja_G =
            factory.createArray<int32_t>({1, sz_nt});
        {
            auto it = out_ja_G.begin();
            for (std::size_t i = 0; i < sz_nt; ++i, ++it)
                *it = static_cast<int32_t>(ja_G_vec[i]);
        }

        // --- coef_G : double, length nterm_G_val ------------------------------
        TypedArray<double> out_coef_G =
            factory.createArray<double>({1, sz_nt});
        std::copy(coef_G_vec.begin(),
                  coef_G_vec.begin() + static_cast<std::ptrdiff_t>(sz_nt),
                  out_coef_G.begin());

        // Work vectors (ja_G_vec, coef_G_vec) destruct here automatically —
        // equivalent to mxDestroyArray(ja_W) / mxDestroyArray(coef_W)

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(out_nterm_G);
        outputs[1] = std::move(out_iat_G);
        outputs[2] = std::move(out_ja_G);
        outputs[3] = std::move(out_coef_G);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 11)
            throwError("LocalFSAI:badInputCount",
                       "Expected 11 input arguments, got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() != 4)
            throwError("LocalFSAI:badOutputCount",
                       "Expected 4 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        // Scalar double inputs: nthread(0)..nterm_M(7)
        for (std::size_t i = 0; i < 8; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE ||
                inputs[i].getNumberOfElements() != 1)
                throwError("LocalFSAI:badScalar",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double scalar.");

        // int64 array input: iat_M(8)
        if (inputs[8].getType() != ArrayType::INT64)
            throwError("LocalFSAI:badArray",
                       "Input argument 9 (iat_M) must be an int64 array.");

        // int32 array input: ja_M(9)
        if (inputs[9].getType() != ArrayType::INT32)
            throwError("LocalFSAI:badArray",
                       "Input argument 10 (ja_M) must be an int32 array.");

        // double array input: coef_M(10)
        if (inputs[10].getType() != ArrayType::DOUBLE)
            throwError("LocalFSAI:badArray",
                       "Input argument 11 (coef_M) must be a double array.");
    }

    // [FIX-E] createCharArray for strings; routes through MATLAB error()
    void throwError(const std::string& id, const std::string& msg)
    {
        getEngine()->feval(u"error", 0,
            std::vector<Array>{
                factory.createCharArray(id),
                factory.createCharArray(msg)
            });
    }
};

//----------------------------------------------------------------------------------------
