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
 * Preconditioned Generalized Minimal Residual (GMRES) for pseudo Energy Minimization,
 * matrix free version.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 *
*****************************************************************************************/
int DEFL_GMRES_matfree(const int np, const int prec_type, const int sol_type, const int nn,
                       const int nn_C, const int nn_K, const int ntv, const int *perm,
                       const int *iperm, const double *D_inv, const int *iat_A,
                       const int *ja_A, const double *coef_A, const int *iat_patt,
                       const int *ja_patt, const int *iat_Tpatt, const int *ja_Tpatt,
                       const double *mat_Q, const double *vec_P0, const double *vec_f,
                       const int itmax, const double energy_tol, int &iter, double *vec_DP,
                       bool verb){

   // Init error code
   int ierr = 0;

   // Variables for energy computation
   double init_energy = 0.0;

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
   for (int i = 0; i < nn_K; i++){
      res[i] = rhs[i];
   }

   // Restart size
   int const m = 50;

   // Krylov basis: (m+1) vectors of size nn_K
   double *V = (double*) aligned_alloc(64, (m+1)*nn_K*sizeof(double));
   
   // Hessenberg matrix: (m+1) x m
   double *H = (double*) aligned_alloc(64, (m+1)*m*sizeof(double));
   
   // Givens + RHS
   double *cs = (double*) malloc(m*sizeof(double));
   double *sn = (double*) malloc(m*sizeof(double));
   double *g  = (double*) malloc((m+1)*sizeof(double));
   
   // Access macros
   #define VEC(j) (&V[(j)*nn_K])
   #define HESS(i,j) H[(j)*(m+1) + (i)]

   // Init GMRES
   iter = 0;
   bool exit_test = (itmax <= 0);

   // Initial residual already in res
   double beta = dnrm2_par(np,nn_K,res,ridv);

   init_energy = beta * beta;
   double Eold, DE0, DEk;

   #if COMP_ENRG
   if(verb){
      Orth_Q(np,nn,ntv,iat_patt,mat_Q,vec_f,v_ntv,vscr);
      double const E0 = dnrm2_par(np,nn_K,vscr,ridv);
      std::cout << std::setprecision(6) << std::scientific;
      std::cout << "Initial Energy:  " << init_energy << std::endl;
   }
   #endif

   // v0 = r / ||r||
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++) {
      VEC(0)[i] = res[i] / beta;
   };

   // initialize g
   for (int i = 0; i < m+1; i++) {
      g[i] = 0.0;
   }
   g[0] = beta;

   while (!exit_test){

      int j;
      for (j = 0; j < m && !exit_test; j++){
         iter++;

         // z = (I - Q Q^t) M^-1 v_j
         apply_perm(np,nn_K,iperm,VEC(j),vscr);

         if (prec_type == DIAG){
            #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nn_K; i++) {
               wscr[i] = D_inv[i]*vscr[i];
            }
         } else if (prec_type == SGS){
            LinvP_spmat(np,nn_C,iat_Tpatt,ja_Tpatt,vscr,nn,iat_A,ja_A,coef_A,WNALLA);
            UinvDP_spmat(np,nn_C,iat_Tpatt,ja_Tpatt,vscr,nn,iat_A,ja_A,coef_A,WNALLA);

            #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nn_K; i++) {
               wscr[i] = vscr[i];
            }
         }

         apply_perm(np,nn_K,perm,wscr,vscr);

         Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,zvec);

         // w = (I - Q Q^t) K z
         if (sol_type == MATFREE){
            apply_perm(np,nn_K,iperm,zvec,vscr);
            mult_K_matfree(np,nn_C,iat_A,ja_A,coef_A,iat_Tpatt,ja_Tpatt,vscr,wscr);
            apply_perm(np,nn_K,perm,wscr,vscr);
         } else {
            KP_spmat(np,nn,iat_A,ja_A,coef_A,nn_C,iat_patt,ja_patt,zvec,vscr,WNALL);
         }

         Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,wscr);

         // Arnoldi
         for (int i = 0; i <= j; i++){

            double const hij = ddot_par(np,nn_K,wscr,VEC(i),ridv);
            HESS(i,j) = hij;

            #pragma omp parallel for num_threads(np)
            for (int k = 0; k < nn_K; k++) {
               wscr[k] -= hij * VEC(i)[k];
            }
         }

         double const hnext = dnrm2_par(np,nn_K,wscr,ridv);
         HESS(j+1,j) = hnext;

         // v_{j+1}
         #pragma omp parallel for num_threads(np)
         for (int k = 0; k < nn_K; k++) {
            VEC(j+1)[k] = wscr[k] / hnext;
         }

         // Apply previous Givens
         for (int i = 0; i < j; i++){
            double const temp = cs[i]*HESS(i,j) + sn[i]*HESS(i+1,j);
            HESS(i+1,j) = -sn[i]*HESS(i,j) + cs[i]*HESS(i+1,j);
            HESS(i,j)   = temp;
         }

         // New Givens
         double const r = std::sqrt(HESS(j,j)*HESS(j,j) + HESS(j+1,j)*HESS(j+1,j));

         cs[j] = HESS(j,j) / r;
         sn[j] = HESS(j+1,j) / r;

         HESS(j,j)   = r;
         HESS(j+1,j) = 0.0;

         // update g
         double const temp = cs[j]*g[j];
         g[j+1] = -sn[j]*g[j];
         g[j]   = temp;

         double const resid = std::abs(g[j+1]);

         double const E = resid * resid;
         DEk = E - Eold;

         if (iter == 1) {
            DE0 = DEk;
            if(verb){
               std::cout << std::setw(4)  << "iter"
                         << std::setw(15) << "Energy"
                         << std::setw(15) << "DE"
                         << std::endl;
            }
         }

         double const dDE = DEk/DE0;

         if(verb){
            std::cout << std::setw(4)  << iter
                      << std::setw(15) << std::scientific << std::setprecision(6) << E
                      << std::setw(15) << std::scientific << std::setprecision(6) << dDE
                      << std::endl;
         }

         // Check convergence
         exit_test = (iter == itmax) || (dDE < energy_tol);
      }

      int const inner_iter = j;

      // Backsolve H y = g
      double y[inner_iter];

      for (int i = inner_iter-1; i >= 0; i--){
         y[i] = g[i];
         for (int k = i+1; k < inner_iter; k++) {
            y[i] -= HESS(i,k)*y[k];
         }
         y[i] /= HESS(i,i);
      }

      // Update solution
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++) {
         vscr[i] = 0.0;
      }

      for (int j = 0; j < inner_iter; j++){
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) {
            vscr[i] += y[j] * VEC(j)[i];
         }
      }

      // z = (I - Q Q^t) M^-1 v_j
      apply_perm(np,nn_K,iperm,vscr,wscr);

      if (prec_type == DIAG){
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) {
            wscr[i] = D_inv[i]*wscr[i];
         }
      } else if (prec_type == SGS){
         LinvP_spmat(np,nn_C,iat_Tpatt,ja_Tpatt,wscr,nn,iat_A,ja_A,coef_A,WNALLA);
         UinvDP_spmat(np,nn_C,iat_Tpatt,ja_Tpatt,wscr,nn,iat_A,ja_A,coef_A,WNALLA);
      }

      apply_perm(np,nn_K,perm,wscr,vscr);

      Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,wscr);

      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++) {
         vec_DP[i] += wscr[i];
      }

      // restart: recompute residual + beta
      // Compute K * x
      if (sol_type == MATFREE){
         apply_perm(np,nn_K,iperm,vec_DP,vscr);
         mult_K_matfree(np,nn_C,iat_A,ja_A,coef_A,iat_Tpatt,ja_Tpatt,vscr,wscr);
         apply_perm(np,nn_K,perm,wscr,vscr);
      } else {
         KP_spmat(np,nn,iat_A,ja_A,coef_A,nn_C,iat_patt,ja_patt,vec_DP,vscr,WNALL);
      }

      // Apply deflation: (I - Q Q^t)
      Orth_Q(np,nn,ntv,iat_patt,mat_Q,vscr,v_ntv,wscr);

      // r = b - deflated(Kx)
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++) {
         res[i] = rhs[i] - wscr[i];
      }

      beta = dnrm2_par(np,nn_K,res,ridv);

      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nn_K; i++) {
         VEC(0)[i] = res[i] / beta;
      }

      for (int i = 0; i < m+1; i++) {
         g[i] = 0.0;
      }
      g[0] = beta;
   }

   free(V);
   free(H);
   free(cs);
   free(sn);
   free(g);

   #if COMP_ENRG
   double const final_energy = beta * beta;
   if(verb){
      std::cout << std::setprecision(6) << std::scientific;
      std::cout << "Final Energy:    " << final_energy << std::endl;
   }
   #endif

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
