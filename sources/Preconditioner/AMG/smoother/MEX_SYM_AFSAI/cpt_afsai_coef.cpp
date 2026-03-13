
#include "cpt_afsai_coef.h"

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


