#include "omp.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace std;

#include "parm_EMIN.h"
#include "spmatv_blk.h"
#include "apply_perm.h"
#include "Orth_Q.h"
#include "ddot_par.h"
#include "dnrm2_par.h"

/*****************************************************************************************
 *
 * Preconditioned Conjugate Gradient for Energy Minimization.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 *
*****************************************************************************************/
int DEFL_PCG(const int np, const int prec_type, const int nn, const int nn_C,
             const int nn_K, const int ntv, const int *perm, const int *iperm,
             const double Tr_A, const int *iat_K, const int *ja_K, const double *coef_K,
             const int *it_U, const int *jcol_U, const double *coef_U,
             const double *D_inv, const int *iat_patt, const int *ja_patt,
             const int *iat_Tpatt, const int *ja_Tpatt, const double *mat_Q,
             const double *vec_P0, const double *vec_f, const int itmax, int &iter,
             double &relres, double *vec_DP){

   // Init error code
   int ierr = 0;

   // Allocate scratches
   double *rhs = (double*) malloc( nn_K*sizeof(double) );
   double *vscr = (double*) malloc( nn_K*sizeof(double) );
   double *wscr = (double*) malloc( nn_K*sizeof(double) );
   double *res = (double*) malloc( nn_K*sizeof(double) );
   double *zvec = (double*) malloc( nn_K*sizeof(double) );
   double *pvec = (double*) malloc( nn_K*sizeof(double) );
   double *QKpvec = (double*) malloc( nn_K*sizeof(double) );
   double *ridv = (double*) malloc( np*sizeof(double) );
   double *v_ntv = (double*) malloc( (np*ntv)*sizeof(double) );
   double *tvec = (double*) malloc( nn_K*sizeof(double) );
   if (rhs == nullptr || vscr == nullptr || wscr == nullptr || res == nullptr ||
       zvec == nullptr || pvec == nullptr || QKpvec == nullptr ||
       ridv == nullptr || v_ntv == nullptr || tvec == nullptr) return ierr = 1;

   // Init vec_DP to zero
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++)
   {
     vec_DP[i] = 0.0;
     tvec[i] = 0.0;
   }

   // Compute rhs
   // 1 - Permute vec_P0 from row-major to col-major
   apply_perm(np,nn_K,iperm,vec_P0,vscr);
   // 2 - Perform wscr <-- K*vscr
   spmatv_blk(np,nn_C,iat_Tpatt,iat_K,ja_K,coef_K,vscr,wscr);
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //#if COMP_ENRG
   double DE_old = 1.0;
   // Allocate back-up vector for entire rhs
   double *rhs_sav = (double*) malloc( nn_K*sizeof(double) );
   if (rhs_sav == nullptr) return ierr = 1;
   double init_energy = 0;
   #if COMP_ENRG
   // Compute initial energy
   init_energy = Tr_A + ddot_par(np,nn_K,vscr,wscr,ridv) -
                        2.0*ddot_par(np,nn_K,vscr,vec_f,ridv);
   cout << setprecision(10) << scientific;
   cout << "INIT ENERGY:  " << init_energy << endl;
   cout << "TRACCIA: " << Tr_A << endl;
   cout << "P0KP0: " << ddot_par(np,nn_K,vscr,wscr,ridv) << endl;
   cout << "-2*f'*P0: " << -2.0*ddot_par(np,nn_K,vscr,vec_f,ridv) << endl;
   #endif
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // 3 - Subtract vec_f to wscr: wscr <-- wscr - vec_f
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++){
      wscr[i] -= vec_f[i];
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //#if COMP_ENRG
      // Save a copy of the rhs in col-major format
      rhs_sav[i] = wscr[i];
      //#endif
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   }
   // 4 - Permute wscr from col-major to row-major
   apply_perm(np,nn_K,perm,wscr,vscr);
   // 5 - Multiply by Q: rhs <-- (I - Q*Q')*vscr
   Orth_Q(np,nn,nn_K,ntv,iat_patt,ja_patt,mat_Q,vscr,v_ntv,rhs);

   // Init residual
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++) res[i] = rhs[i];
   double bnorm = dnrm2_par(np,nn_K,rhs,ridv);
   //////////////////////////////////////////
   //cout << "BNORM: " << bnorm << endl;
   //cout << "ITMAX: " << itmax << endl;
   //////////////////////////////////////////

   // Init PCG
   iter = 0;
   double resiter = 1.0;
   bool exit_test = (itmax <= 0);

   // PCG loop
   double alpha, beta, gamma, gamma_old;
   while (!exit_test){

      // Increase iter count
      iter++;

      // Compute zvec
      // 1 - Permute res from row-major to col-major
      apply_perm(np,nn_K,iperm,res,vscr);
      // 2 - Apply preconditioner to wscr: wscr <-- M_inv*vscr
      if (prec_type == DIAG){
         #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nn_K; i++) wscr[i] = D_inv[i]*vscr[i];
      }
      // 3 - Permute wscr from col-major to row-major
      apply_perm(np,nn_K,perm,wscr,vscr);
      // 4 - Multiply by Q: zvec <-- (I - Q*QT)*vscr
      Orth_Q(np,nn,nn_K,ntv,iat_patt,ja_patt,mat_Q,vscr,v_ntv,zvec);

      // Compute gamma <-- resT*zvec
      gamma = ddot_par(np,nn_K,res,zvec,ridv);
      //////////////////////////////////
      //cout << "ITER " << iter << endl;
      //cout << "gamma " << gamma << endl;
      //////////////////////////////////

      // Compute pvec
      if (iter == 1){
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) pvec[i] = zvec[i];
      } else {
         beta = gamma / gamma_old;
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) pvec[i] = zvec[i] + beta*pvec[i];
      }
      //////////////////////////////////
      //if (iter > 1) cout << "beta " << beta << endl;
      //cout << "NORM(P) " << dnrm2_par(np,nn_K,pvec,ridv) << endl;
      //////////////////////////////////

      // Save gamma
      gamma_old = gamma;

      // Premultiply the search direction by the operator
      // 1 - Permute pvec from row-major to col-major
      apply_perm(np,nn_K,iperm,pvec,vscr);
      // 2 - Multiply by K: wscr <-- K*vscr
      spmatv_blk(np,nn_C,iat_Tpatt,iat_K,ja_K,coef_K,vscr,wscr);
      // 3 - Permute wscr from col-major to row-major
      apply_perm(np,nn_K,perm,wscr,vscr);
      // 4 - Multiply by Q: QKpvec <-- (I - Q*QT)*vscr
      Orth_Q(np,nn,nn_K,ntv,iat_patt,ja_patt,mat_Q,vscr,v_ntv,QKpvec);

      // Compute alpha
      alpha = ddot_par(np,nn_K,QKpvec,pvec,ridv);
      alpha = gamma / alpha;
      //////////////////////////////////
      //cout << "alpha " << alpha << endl;
      //////////////////////////////////

      // Update Prolongation and residual
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++){
         vec_DP[i] += alpha*pvec[i];
         res[i] -= alpha*QKpvec[i];
         tvec[i] += alpha*wscr[i];
      }

      ////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //#if COMP_ENRG
      //   apply_perm(np,nn_K,iperm,vec_DP,vscr);
      //   double rhs_DP = ddot_par(np,nn_K,vscr,rhs_sav,ridv);
      //   spmatv_blk(np,nn_C,iat_Tpatt,iat_K,ja_K,coef_K,vscr,wscr);
      //   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //   //double DWAffDW = ddot_par(np,nn_K,vscr,wscr,ridv);
      //   //cout << endl << endl;
      //   //cout << "Afc_DP  " << Afc_DP << endl;
      //   //cout << "DWAffDW  " << DWAffDW << endl;
      //   //cout << endl << endl;
      //   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //   double DE = ddot_par(np,nn_K,vscr,wscr,ridv) - 2.0*rhs_DP;
      //   cout << fixed << setw(6) << iter << " ";
      //   cout << setprecision(6) << scientific;
      //   cout << resiter << "   " << init_energy << "   " << DE << "   " <<
      //           fabs(DE-DE_old) <<  "   " << init_energy+DE << endl;
      //   DE_old = DE;
      //#else
      //   cout << setprecision(6) << scientific;
      //   cout << iter << " " << resiter << endl;
      //#endif

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Compute energy and energy variation
      apply_perm(np,nn_K,iperm,vec_DP,vscr);
      double const rhs_DP = ddot_par(np,nn_K,vscr,rhs_sav,ridv);
      double const DE = ddot_par(np,nn_K,vscr,tvec,ridv) - 2.0*rhs_DP;
      double const dDE = std::abs((DE-DE_old)/DE_old);
      cout << fixed << setw(6) << iter << " ";
      cout << setprecision(6) << scientific;
      cout << resiter << "   " << init_energy << "   " << DE << "   " <<
              dDE <<  "   " << init_energy+DE << endl;
      DE_old = DE;

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Check convergence
      resiter = dnrm2_par(np,nn_K,res,ridv) / bnorm;
      exit_test = (resiter < EMIN_CGTOL) || (iter == itmax) || (dDE < ENER_TOL);

   }
   // Store final relative residual
   relres = resiter;

   // Deallocate local scratches
   free(v_ntv);

   // Deallocate scratches
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //#if COMP_ENRG
   free(rhs_sav);
   //#endif
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   free(rhs);
   free(res);
   free(vscr);
   free(wscr);
   free(zvec);
   free(pvec);
   free(QKpvec);
   free(ridv);
   free(tvec);

   return ierr;

}
