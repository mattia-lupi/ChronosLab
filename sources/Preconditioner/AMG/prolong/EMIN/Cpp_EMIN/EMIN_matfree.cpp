#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>

#if defined PRINT
#define dump true
#else
#define dump false
#endif

#include "parm_EMIN.h"
#include "wrCSRmat.h"
#include "Transp_Patt.h"
#include "load_Jacobi.h"
#include "copy_Prol.h"
#include "gather_B_QR.h"
#include "gather_B_dump.h"
#include "gather_f.h"
#include "cpt_Trace_Acc.h"
#include "print_Q.h"
#include "Prol_add_Cnodes.h"
#include "DEFL_PCG_matfree.h"

/*****************************************************************************************
 *
 * Function that starting from an input prolongation and a given non-zero pattern,
 * improves the prolongation through an energy minimization approach while maintaining the
 * near kernel included in the prolongation range. This routine works in matfree mode.
 *
 * Input:
 *
 * np:                              number of openMP threads.
 * itmax:                           number of EnerMinCG iterations.
 * condmax:                         max conditioning allowed for a B block.
 * prec_type:                       preconditioner for energy minimization
 * sol_type:                        solution method for energy minimization
 * nn:                              size of the system matrix.
 * nn_C:                            number of coarse nodes.
 * ntv:                             number of test vectors.
 * nt_patt:                         number of non-zeroes of the non-zero pattern.
 * fcnode:                          F/C indicator.
 * iat_A, ja_A, coef_A:             system matrix.
 * iat_Pin, ja_Pin, coef_Pin:       input prolongation.
 * iat_patt, ja_patt:               prescribed prolongation pattern.
 * TV:                              nn x ntv test space.
 *
 * Output:
 *
 * iat_Pout, ja_Pout, coef_Pout:    output prolongation
 * info:                            time and nonzero info
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> error in gather_K
 * 3 ---> error in preconditioner computation
 * 4 ---> error in gather_B
 * 5 ---> error in PCG
 * 6 ---> allocation for final result
 *
 * NOTE: It is assumed that the input pattern does not contain row entries relative
 * to coarse nodes.
 *
 *****************************************************************************************/

int EMIN_matfree(const int np, const int itmax, const double en_tol, const double condmax,
                 const int prec_type, const int sol_type, const int nn, const int nn_C,
                 const int ntv, const int nt_patt, const int *fcnode, const int *iat_A,
                 const int *ja_A, const double *coef_A, const int *iat_Pin, const int *ja_Pin,
                 const double *coef_Pin, const int *iat_patt, const int *ja_patt,
                 const double *const *TV, int *&iat_Pout, int *&ja_Pout, double *&coef_Pout,
                 double *info, bool verb)
{

   // Init error code
   int ierr = 0;

   // --- Local variables for timing -----------------------------------------------------
   std::chrono::time_point<std::chrono::system_clock> start, end;
   std::chrono::duration<double> elaps_sec;

   // ------ Transpose the pattern ----------------------------------------------------

   int *iat_Tpatt = (int*) malloc( (nn_C+1)*sizeof(int) );
   int *ja_Tpatt  = (int*) malloc( nt_patt*sizeof(int) );
   int *perm  = (int*) malloc( nt_patt*sizeof(int) );
   int *iperm  = (int*) malloc( nt_patt*sizeof(int) );
   if (iat_Tpatt == nullptr || ja_Tpatt == nullptr || perm == nullptr ||
       iperm == nullptr)
      return ierr = 1;
   ierr = Transp_Patt(np,nn,nn_C,nt_patt,iat_patt,ja_patt,iat_Tpatt,ja_Tpatt,perm,
                      iperm);
   if (ierr != 0) return ierr = 1;

   // ----- Compute the preconditioner ---------------------------------------------------

   int nnz_J = 0;
   double *D_inv = nullptr;
   double time_prec_K = 0.0;
   if (prec_type == DIAG){
      // Compute a Jacobi preconditioner for K
      start = std::chrono::system_clock::now();
      D_inv = (double*) malloc( nt_patt*sizeof(double) );
      double *scr = (double*) malloc( nn*sizeof(double) );
      if (D_inv == nullptr || scr == nullptr) return ierr = 3;
      load_Jacobi(np,nn,nt_patt,ja_Tpatt,iat_A,ja_A,coef_A,scr,D_inv);
      free(scr);
      nnz_J = nt_patt;
      end = std::chrono::system_clock::now();
      elaps_sec = end - start;
      time_prec_K = elaps_sec.count();
      if (DUMP_PREC){
         FILE *dfile = fopen("D_mat","w");
         for (int k = 0; k < nnz_J; k++) fprintf(dfile,"%20.11e\n",D_inv[k]);
         fflush(dfile);
         fclose(dfile);
      }
   }

   // ----- Copy the initial prolongation into the one with extended pattern
   if (dump) std::cout << "---- COPY PROL ----" << std::endl << std::endl;
   double *coef_P0 = (double*) calloc( nt_patt , sizeof(double) );
   if (coef_P0 == nullptr) return ierr = 1;
   copy_Prol(np,nn,fcnode,iat_Pin,ja_Pin,coef_Pin,iat_patt,ja_patt,coef_P0); 
   if (DUMP_PREC){
      FILE *origPfile = fopen("origProl.csr","w");
      wrCSRmat(origPfile,false,nn,iat_Pin,ja_Pin,coef_Pin);
      fclose(origPfile);
      FILE *extPfile = fopen("extProl.csr","w");
      wrCSRmat(extPfile,false,nn,iat_patt,ja_patt,coef_P0);
      fclose(extPfile);
   }

   // ----- Assemble and decompose with QR the constraint part B ----------------------
   if (dump) std::cout << "---- gather_B_Qr ----" << std::endl << std::endl;
   start = std::chrono::system_clock::now();
   double *mat_Q = nullptr;
   double *vec_f = nullptr;
   if (DUMP_PREC){
      double *mat_B = nullptr;
      ierr = gather_B_dump(np,nn,nn_C,ntv,fcnode,iat_patt,ja_patt,TV,mat_B);
      FILE *Bfile = fopen("mat_B","w");
      print_Q(Bfile,nn,ntv,iat_patt,mat_B);
      fflush(Bfile);
      fclose(Bfile);
      free(mat_B);
   }
   ierr = gather_B_QR(np,condmax,nn,nn_C,ntv,fcnode,iat_patt,ja_patt,TV,mat_Q,coef_P0);
   if (ierr != 0) return ierr = 4;
   int nnz_Q = ntv*iat_patt[nn];

   int nrows_Q = iat_patt[nn];
   vec_f = (double*) malloc( nrows_Q*sizeof(double) );
   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   if (vec_f == nullptr || c2glo == nullptr) return ierr = 1;
    // Create mapping from coarse node numbering to global (original) numbering
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn; i++){
      int k = fcnode[i];
      if (k >= 0) c2glo[k] = i;
   }
   // Gather vec_f
   #pragma omp parallel for num_threads(np)
   for (int icol = 0; icol < nn_C; icol++){
       int istart = iat_Tpatt[icol];
       int iend = iat_Tpatt[icol+1];
       int irow = c2glo[icol];
       int istart_A = iat_A[irow];
       int iend_A = iat_A[irow+1];
       gather_f(iend-istart,&(ja_Tpatt[istart]),iend_A-istart_A,&(ja_A[istart_A]),
                &(coef_A[istart_A]),&(vec_f[istart]));
   }
   end = std::chrono::system_clock::now();
   elaps_sec = end - start;
   free(c2glo);
   double time_gath_B = elaps_sec.count();
   if (DUMP_PREC){
      FILE *corrPfile = fopen("corrProl.csr","w");
      wrCSRmat(corrPfile,false,nn,iat_patt,ja_patt,coef_P0);
      fclose(corrPfile);
      FILE *Qfile = fopen("mat_Q2","w");
      print_Q(Qfile,nn,ntv,iat_patt,mat_Q);
      fflush(Qfile);
      fclose(Qfile);
      FILE *ffile = fopen("vec_f","w");
      for (int i = 0; i < nt_patt; i++) fprintf(ffile,"%15.6e\n",vec_f[i]);
      fflush(ffile);
      fclose(ffile);
   }
   double Tr_A = 0;
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   #if COMP_ENRG
   // Compute Trace of A
   Tr_A = cpt_Trace_Acc(np,nn,fcnode,iat_A,ja_A,coef_A);
   #endif
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   int iter;

   // Compute prolongation correction with PCG
   if (dump) std::cout << "---- DEFL_PCG_matfree ----" << std::endl << std::endl;
   start = std::chrono::system_clock::now();
   // Allocate prolongation correction
   double *DP = (double*) malloc( nt_patt*sizeof(double) );
   if (DP == nullptr) return ierr = 1;
   ierr = DEFL_PCG_matfree(np,prec_type,sol_type,nn,nn_C,nt_patt,ntv,perm,iperm,D_inv,
                           iat_A,ja_A,coef_A,Tr_A,iat_patt,ja_patt,iat_Tpatt,ja_Tpatt,
                           mat_Q,coef_P0,vec_f,itmax,en_tol,iter,DP,verb);
   if (ierr != 0) return ierr = 5;
   end = std::chrono::system_clock::now();
   elaps_sec = end - start;
   double time_PCG = elaps_sec.count();
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if(verb){
      std::cout << "PCG TIME "<< time_PCG << std::endl;
   }
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Update prolongation with DP
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nt_patt; i++) coef_P0[i] -= DP[i];
   
   // Free scracthes
   free(DP);
   free(vec_f);
   free(mat_Q);
   free(D_inv);
   free(perm);
   free(iperm);
   free(iat_Tpatt);
   free(ja_Tpatt);

   // Extend prolongation including coarse nodes
   ierr = Prol_add_Cnodes(np,nn,nn_C,fcnode,iat_patt,ja_patt,coef_P0,iat_Pout,ja_Pout,
                          coef_Pout);
   if (ierr != 0) ierr = 6;

   // Free also coef_P0
   free(coef_P0);

   // Store info
   info[0]  = 0.0;
   info[1]  = time_prec_K;
   info[2]  = time_gath_B;
   info[3]  = time_PCG;
   info[6]  = 0.0;
   info[7]  = static_cast<double>(iter);
   info[8]  = 0.0;
   info[9]  = 0.0;
   info[10] = static_cast<double>(nnz_J);
   info[11] = static_cast<double>(nnz_Q);

   return ierr;

}
