#include "omp.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "parm_EMIN.h"
#include "mult_K_matfree.h"
#include "apply_perm.h"
#include "Orth_Q.h"
#include "ddot_par.h"
#include "dnrm2_par.h"
#include "KP_spmat.h"
#include "LinvP_spmat.h"
#include "UinvDP_spmat.h"

/*****************************************************************************************
 *
 * Preconditioned Conjugate Gradient for Energy Minimization, matrix free version.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 *
*****************************************************************************************/
int DEFL_PCG_matfree(const int np, const int prec_type, const int sol_type, const int nn,
                     const int nn_C, const int nn_K, const int ntv, const int *perm,
                     const int *iperm, const double *D_inv, const int *iat_A,
                     const int *ja_A, const double *coef_A, const double Tr_A,
                     const int *iat_patt, const int *ja_patt, const int *iat_Tpatt,
                     const int *ja_Tpatt, const double *mat_Q, const double *vec_P0,
                     const double *vec_f, const int itmax, const double energy_tol,
                     int &iter, double *vec_DP){

   // Init error code
   int ierr = 0;

   // Variables for energy computation
   double init_energy = 0.0;
   double DE = 0, DEk, DE0;

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
   //double *tvec = (double*) malloc( nn_K*sizeof(double) );
   if (rhs == nullptr || vscr == nullptr || wscr == nullptr || res == nullptr ||
       zvec == nullptr || pvec == nullptr || QKpvec == nullptr ||
       ridv == nullptr || v_ntv == nullptr) return ierr = 1;
       //ridv == nullptr || v_ntv == nullptr || tvec == nullptr) return ierr = 1;
   int *WNALL = nullptr;
   int *WNALLA = nullptr;
   if (sol_type == SPMAT){
      WNALL = (int*) malloc( (np*nn_C)*sizeof(int) );
      if (WNALL == nullptr) return ierr = 1;
      WNALLA = (int*) malloc( (np*nn)*sizeof(int) );
      if (WNALLA == nullptr) return ierr = 1;
   }

   // Init vec_DP to zero
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++)
   {
     vec_DP[i] = 0.0;
   }

   // Compute rhs
   if (sol_type == MATFREE){
      // 1 - Permute vec_P0 from row-major to col-major
      apply_perm(np,nn_K,iperm,vec_P0,vscr);
      // 2 - Perform wscr <-- K*vscr
      mult_K_matfree(np,nn_C,iat_A,ja_A,coef_A,iat_Tpatt,ja_Tpatt,vscr,wscr);
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #if COMP_ENRG
      // Compute initial energy
      init_energy = Tr_A + ddot_par(np,nn_K,vscr,wscr,ridv) -
                           2.0*ddot_par(np,nn_K,vscr,vec_f,ridv);
      std::cout << std::setprecision(6) << std::scientific;
      std::cout << "Initial Energy:  " << init_energy << std::endl;
      std::cout << "Trace of A:      " << Tr_A << std::endl;
      #endif
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // 3 - Subtract vec_f to wscr: wscr <-- wscr - vec_f
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++){
         wscr[i] -= vec_f[i];
      }
      // 4 - Permute wscr from col-major to row-major
      apply_perm(np,nn_K,perm,wscr,vscr);
   } else {
      // 1 - No need for permutation
      // 2 - Perform vscr <-- K*vec_P0
      KP_spmat(np,nn,iat_A,ja_A,coef_A,nn_C,iat_patt,ja_patt,vec_P0,vscr,WNALL);

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #if COMP_ENRG
      // Permute vec_f from col-major to row-major
      apply_perm(np,nn_K,perm,vec_f,wscr);
      // Compute initial energy
      init_energy = Tr_A + ddot_par(np,nn_K,vec_P0,vscr,ridv) -
                           2.0*ddot_par(np,nn_K,vec_P0,wscr,ridv);
      std::cout << std::setprecision(6) << std::scientific;
      std::cout << "Initial Energy:  " << init_energy << std::endl;
      std::cout << "Trace of A:      " << Tr_A << std::endl;
      #endif
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // 3 - Subtract vec_f to vscr: vscr <-- vscr - vec_f
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++){
         vscr[i] -= vec_f[perm[i]];//CHECK
      }
      // 4 - No need for permutation
   }
   // 5 - Multiply by Q: rhs <-- (I - Q*Q')*vscr
   Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,rhs);

   // Init residual
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++) res[i] = rhs[i];

   // Init PCG
   iter = 0;
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
      } else if (prec_type == SGS){
        // symmetric Gauss-Seidel U^-1 * D * L^-1 * vscr
        // works in place
        LinvP_spmat(np,nn_C,iat_Tpatt,ja_Tpatt,vscr,nn,iat_A,ja_A,coef_A,WNALLA);
        UinvDP_spmat(np,nn_C,iat_Tpatt,ja_Tpatt,vscr,nn,iat_A,ja_A,coef_A,WNALLA);
        // TODO this can be removed ...
        #pragma omp parallel for num_threads(np)
          for (int i = 0; i < nn_K; i++) wscr[i] = vscr[i];
      }
      // 3 - Permute wscr from col-major to row-major
      apply_perm(np,nn_K,perm,wscr,vscr);
      // 4 - Multiply by Q: zvec <-- (I - Q*QT)*vscr
      Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,zvec);

      // Compute gamma <-- resT*zvec
      gamma = ddot_par(np,nn_K,res,zvec,ridv);

      // Compute pvec
      if (iter == 1){
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) pvec[i] = zvec[i];
      } else {
         beta = gamma / gamma_old;
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) pvec[i] = zvec[i] + beta*pvec[i];
      }

      // Save gamma
      gamma_old = gamma;

      // Premultiply the search direction by the operator
      if (sol_type == MATFREE){
         // 1 - Permute pvec from row-major to col-major
         apply_perm(np,nn_K,iperm,pvec,vscr);
         // 2 - Multiply by K: wscr <-- K*vscr
         mult_K_matfree(np,nn_C,iat_A,ja_A,coef_A,iat_Tpatt,ja_Tpatt,vscr,wscr);
         // 3 - Permute wscr from col-major to row-major
         apply_perm(np,nn_K,perm,wscr,vscr);
      } else {
         // 1 - No need for permutation
         // 2 - Multiply by K: vscr <-- K*pvec
         KP_spmat(np,nn,iat_A,ja_A,coef_A,nn_C,iat_patt,ja_patt,pvec,vscr,WNALL);
         // 3 - No need for permutation
      }
      // 4 - Multiply by Q: QKpvec <-- (I - Q*QT)*vscr
      Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,QKpvec);

      // Compute alpha
      alpha = ddot_par(np,nn_K,QKpvec,pvec,ridv);
      // Compute energy reduction
      DEk = gamma*gamma / alpha;
      alpha = gamma / alpha;

      // Update Prolongation and residual
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++){
         vec_DP[i] += alpha*pvec[i];
         res[i] -= alpha*QKpvec[i];
      }

      if (iter==1){
         DE0 = DEk;
         printf("%4s %15s %15s\n","iter","Energy","DE");
      }
      double const dDE = DEk/DE0;
      DE -= DEk;
      printf("%4d %15.6e %15.6e\n",iter,init_energy+DE,dDE);

      // Check convergence
      exit_test = (iter == itmax) || (dDE < energy_tol);

   }

   // Deallocate local scratches
   free(v_ntv);

   // Deallocate scratches
   free(rhs);
   free(res);
   free(vscr);
   free(wscr);
   free(zvec);
   free(pvec);
   free(QKpvec);
   free(ridv);
   free(WNALL);
   free(WNALLA);

   return ierr;

}
