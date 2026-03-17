#include "omp.h"
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "parm_EMIN.h"
#include "spmatv_blk.h"
#include "apply_perm.h"
#include "Z_mult.h"
#include "ZT_mult.h"
#include "ddot_par.h"
#include "dnrm2_par.h"

/*****************************************************************************************
 *
 * Preconditioned Conjugate Gradient for Energy Minimization.
 * vec_P contain initial solution that is corrected with Z*sol
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 *
*****************************************************************************************/
int NLSP_PCG(const int np, const int prec_type, const int nn, const int nn_C,
             const int nn_K, const int ntv, const int *perm, const int *iperm,
             const double Tr_A, const int *iat_K, const int *ja_K, const double *coef_K,
             const double *D_inv, const int *iat_patt, const int *iat_Tpatt,
             const int *pt_Z, const int *pt_col_Z, const double *mat_Z, const double *vec_f,
             const int itmax, int &iter, double &relres, double *vec_P){

   // Init error code
   int ierr = 0;

   // Get size of nullspace system
   int nequ = pt_col_Z[nn];

   // Allocate scratches @@@@@@@@@MOLTI VETTORI POSSONO ESSERE RIDOTTI
   double *sol = (double*) malloc( nequ*sizeof(double) );
   double *rhs = (double*) malloc( nequ*sizeof(double) );
   double *vscr = (double*) malloc( nn_K*sizeof(double) );
   double *wscr = (double*) malloc( nn_K*sizeof(double) );
   double *res = (double*) malloc( nequ*sizeof(double) );
   double *pres = (double*) malloc( nequ*sizeof(double) );
   double *pvec = (double*) malloc( nequ*sizeof(double) );
   double *Kxp = (double*) malloc( nequ*sizeof(double) );
   double *ridv = (double*) malloc( np*sizeof(double) );
   if (rhs == nullptr || vscr == nullptr || wscr == nullptr || res == nullptr ||
       pres == nullptr || pvec == nullptr || Kxp == nullptr || ridv == nullptr)
       return ierr = 1;

   // Compute rhs
   // 1 - Permute vec_P from row-major to col-major
   apply_perm(np,nn_K,iperm,vec_P,vscr);
   // 2 - Perform vscr <-- K*vec_P
   spmatv_blk(np,nn_C,iat_Tpatt,iat_K,ja_K,coef_K,vscr,wscr);
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   #if COMP_ENRG
   double init_energy = Tr_A + ddot_par(np,nn_K,vscr,wscr,ridv) -
                               2.0*ddot_par(np,nn_K,vscr,vec_f,ridv);
   std::cout << std::setprecision(10) << std::scientific;
   std::cout << "INIT ENERGY:  " << init_energy << std::endl;
   std::cout << "TRACCIA: " << Tr_A << std::endl;
   std::cout << "P0KP0: " << ddot_par(np,nn_K,vscr,wscr,ridv) << std::endl;
   std::cout << "-2*f'*P0: " << -2.0*ddot_par(np,nn_K,vscr,vec_f,ridv) << std::endl;
   #endif
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // 3 - Compute wscr <-- vec_f - wscr = vec_f - K*vec_P
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++) wscr[i] = vec_f[i] - wscr[i];
   // 4 - Permute wscr from col-major to row-major
   apply_perm(np,nn_K,perm,wscr,vscr);
   // 5 - Compute rhs <-- ZT*vscr
   ZT_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,vscr,rhs);

   // Init residual and solution
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nequ; i++){
      res[i] = rhs[i];
      sol[i] = 0.0;
   }
   double bnorm = dnrm2_par(np,nequ,rhs,ridv);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //FILE *rf = fopen("res0","w");
   //for (int i = 0; i < nequ; i++) fprintf(rf,"%20.11e\n",res[i]);
   //fflush(rf);
   //fclose(rf);
   //cout << "BNORM: " << bnorm << endl;
   //cout << "ITMAX: " << itmax << endl;
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Init PCG
   iter = 0;
   double resiter = 1.0;
   bool exit_test = (itmax <= 0);

   // PCG loop
   double alpha, beta, ptKp;
   while (!exit_test){

      // Increase iter count
      iter++;

      // Compute pres
      // 1 - Perform vscr <-- Z*res
      Z_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,res,vscr);
      // 2 - Permute vscr from row-major to col-major
      apply_perm(np,nn_K,iperm,vscr,wscr);
      // 3 - Apply preconditioner to wscr: vscr <-- M_inv*wscr
      if (prec_type == DIAG){
         #pragma omp parallel for num_threads(np)
            for (int i = 0; i < nn_K; i++) vscr[i] = D_inv[i]*wscr[i];
      }
      // 4 - Permute vscr from col-major to row-major
      apply_perm(np,nn_K,perm,vscr,wscr);
      // 5 - Perform pvec <-- ZT*wscr
      ZT_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,wscr,pres);

      // Compute beta and pvec
      if (iter == 1){
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nequ; i++) pvec[i] = pres[i];
      } else {
         beta = - ddot_par(np,nequ,pres,Kxp,ridv) / ptKp;
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nequ; i++) pvec[i] = pres[i] + beta*pvec[i];
      }

      // Multiply pvec by ZT*K*Z
      // 1 - Perform vscr <-- Z*pvec
      Z_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,pvec,vscr);
      // 2 - Permute vscr from row-major to col-major
      apply_perm(np,nn_K,iperm,vscr,wscr);
      // 3 - Multiply by K: vscr <-- K*wscr
      spmatv_blk(np,nn_C,iat_Tpatt,iat_K,ja_K,coef_K,wscr,vscr);
      // 4 - Permute vscr from col-major to row-major
      apply_perm(np,nn_K,perm,vscr,wscr);
      // 5 - Perform Kxp <-- ZT*wscr
      ZT_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,wscr,Kxp);

      // Compute alpha
      ptKp = ddot_par(np,nequ,pvec,Kxp,ridv);
      alpha = ddot_par(np,nequ,pvec,res,ridv) / ptKp;

      // Update solution and residual
      #pragma omp parallel for num_threads(np)
      for (int i = 0; i < nequ; i++){
         sol[i] += alpha*pvec[i];
         res[i] -= alpha*Kxp[i];
      }

      // Check convergence
      resiter = dnrm2_par(np,nequ,res,ridv) / bnorm;
      exit_test = (resiter < EMIN_CGTOL) || (iter == itmax);
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      #if COMP_ENRG
         Z_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,sol,vscr);
         apply_perm(np,nn_K,iperm,vscr,wscr);
         spmatv_blk(np,nn_C,iat_Tpatt,iat_K,ja_K,coef_K,wscr,vscr);
         apply_perm(np,nn_K,perm,vscr,wscr);
         ZT_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,wscr,vscr);
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
         std::cout << std::endl << std::endl;
         std::cout << "P Aff P " << ddot_par(np,nequ,vscr,sol,ridv) << std::endl;
         std::cout << "P rhs " << ddot_par(np,nequ,rhs,sol,ridv) << std::endl;
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nequ; i++)
            vscr[i] -= 2.0*rhs[i];
         double DE = ddot_par(np,nequ,vscr,sol,ridv);
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
         std::cout << "DE   " << DE << std::endl << std::endl;
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
         std::cout << std::fixed << std::setw(6) << iter << " ";
         std::cout << std::setprecision(6) << std::scientific;
         std::cout << resiter << "   " << init_energy << "   " << DE << "   " << init_energy+DE << std::endl;
      #else
         std::cout << std::setprecision(6) << std::scientific;
         std::cout << iter << " " << resiter << std::endl;
      #endif
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   }
   // Store final relative residual
   relres = resiter;
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //FILE *of = nullptr;
   //of = fopen("Pvec0","w");
   //for (int i = 0; i < nn_K; i++) fprintf(of,"%20.11e\n",vec_P[i]);
   //fflush(of);
   //fclose(of);
   //FILE *vf = fopen("v_sol","w");
   //for (int i = 0; i < nequ; i++) fprintf(vf,"%20.11e\n",sol[i]);
   //fflush(vf);
   //fclose(vf);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Compute vec_P = vec_P + Z*sol
   // 1 - Compute vscr <-- Z*sol
   Z_mult(np,nn,ntv,iat_patt,pt_Z,pt_col_Z,mat_Z,sol,vscr);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //of = fopen("DP","w");
   //for (int i = 0; i < nn_K; i++) fprintf(of,"%20.11e\n",vscr[i]);
   //fflush(of);
   //fclose(of);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   // 2 - Compute vec_P <-- vec_P + vscr
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_K; i++) vec_P[i] += vscr[i];

   // Deallocate scratches
   free(sol);
   free(rhs);
   free(vscr);
   free(wscr);
   free(res);
   free(pres);
   free(pvec);
   free(Kxp);
   free(ridv);

   return ierr;

}
