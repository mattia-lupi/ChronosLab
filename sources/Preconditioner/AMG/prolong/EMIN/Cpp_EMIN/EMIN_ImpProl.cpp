#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>
using namespace std;

#if defined PRINT
#define dump true
#else
#define dump false
#endif

#include "parm_EMIN.h"
#include "wrCSRmat.h"
#include "Transp_Patt.h"
#include "EMIN_matfree.h"
#include "gather_K.h"
#include "copy_Prol.h"
#include "gather_B_QR.h"
#include "gather_B_Z.h"
#include "gather_f.h"
#include "cpt_Trace_Acc.h"
#include "print_Q.h"
#include "print_Z.h"
#include "DEFL_PCG.h"
#include "NLSP_PCG.h"
#include "Prol_add_Cnodes.h"

/*****************************************************************************************
 *
 * Function that starting from an input prolongation and a given non-zero pattern,
 * improves the prolongation through an energy minimization approach while maintaining the
 * near kernel included in the prolongation range.
 *
 * Input:
 *
 * np:                              number of openMP threads.
 * itmax:                           number of EnerMinCG iterations.
 * en_tol:                          exit tolerance controlling energy decrease.
 * condmax:                         max conditioning allowed for a B block.
 * prec_type:                       preconditioner for energy minimization
 * sol_type:                        solver for energy minimization
 * nn:                              size of the system matrix.
 * nn_C:                            number of coarse nodes.
 * ntv:                             number of test vectors.
 * nt_A:                            number of non-zeroes of the system matrix.
 * nt_P:                            number of non-zeroes of the input prolongation.
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
 * info:                            array with timings and information on EMIN process
 *                                  info[0]  --> time_gath_K
 *                                  info[1]  --> time_prec_K
 *                                  info[2]  --> time_gath_B
 *                                  info[3]  --> time_PCG
 *                                  info[4]  --> time_overhead
 *                                  info[5]  --> time_glob
 *                                  info[6]  --> avg_lfil
 *                                  info[7]  --> iter
 *                                  info[8]  --> relres
 *                                  info[9]  --> nnz_K
 *                                  info[10] --> nnz_PK
 *                                  info[11] --> nnz_ZQ
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

int EMIN_ImpProl(const int np, const int itmax, const double en_tol, const double condmax,
                 const int prec_type, const int sol_type, const int nn, const int nn_C,
                 const int ntv, const int nt_A, const int nt_P, const int nt_patt,
                 const int *fcnode, const int *iat_A, const int *ja_A, const double *coef_A,
                 const int *iat_Pin, const int *ja_Pin, const double *coef_Pin,
                 const int *iat_patt, const int *ja_patt, const double *const *TV,
                 int *&iat_Pout, int *&ja_Pout, double *&coef_Pout, double *info)
{

   // Init error code
   int ierr = 0;

   // Variables to store info
   int iter = 0;
   int nnz_K = 0, nnz_PK = 0, nnz_ZQ = 0;
   double avg_lfil = 0.0;
   double relres = 0.0;

   // --- Local variables for timing -----------------------------------------------------
   double  time_gath_K = 0.0,time_prec_K = 0.0, time_gath_B = 0.0, time_PCG = 0.0;
   chrono::time_point<std::chrono::system_clock> start, end, glob_start, glob_end;
   chrono::duration<double> elaps_sec;

   //---GLOBAL START-----------------------------
   glob_start = chrono::system_clock::now();
   //--------------------------------------------
   
   if (sol_type == MATFREE || sol_type == SPMAT){

      //**********************************
      //** SOLUTION IN MATRIX FREE MODE **
      //**********************************

      ierr = EMIN_matfree(np,itmax,en_tol,condmax,prec_type,sol_type,nn,nn_C,
                          ntv,nt_A,nt_P,nt_patt,fcnode,iat_A,ja_A,coef_A,iat_Pin,
                          ja_Pin,coef_Pin,iat_patt,ja_patt,TV,iat_Pout,ja_Pout,
                          coef_Pout,info);

   } else {

      //*****************************
      //** SOLUTION STORING MATRIX **
      //*****************************

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

      // ----- Gather the K matrix -------------------------------------------------------

      if (dump) cout << "---- GATHER K ----" << endl << endl;
      start = chrono::system_clock::now();
      int nn_K;
      int *iat_K = nullptr, *ja_K = nullptr;
      double *coef_K = nullptr;
      int nnmax_blk, ntmax_blk;
      ierr = gather_K(np,nn,nn_C,iat_A,ja_A,coef_A,iat_Tpatt,ja_Tpatt,nnmax_blk,ntmax_blk,
                      nn_K,iat_K,ja_K,coef_K);
      if (ierr != 0) return ierr = 2;
      end = chrono::system_clock::now();
      elaps_sec = end - start;
      time_gath_K = elaps_sec.count();
      if (DUMP_PREC){
         FILE *kfile = fopen("MATK","w");
         wrCSRmat(kfile,false,nn_K,iat_K,ja_K,coef_K);
         fflush(kfile);
         fclose(kfile);

         FILE *pfile = fopen("perm","w");
         for( int i=0; i<nt_patt; ++i )
         {
           fprintf(pfile, "%6i\n", perm[i]+1);
         }
         fflush(pfile);
         fclose(pfile);
      }
      nnz_K = iat_K[nn_K];

      // ----- Compute the preconditioner for K ------------------------------------------

      start = chrono::system_clock::now();
      int *it_U = nullptr, *jcol_U = nullptr;
      double *coef_U = nullptr, *D_inv = nullptr;
      if (prec_type == DIAG){
         // Compute a Jacobi preconditioner for K
         D_inv = (double*) malloc( nn_K*sizeof(double) );
         if (D_inv == nullptr) return ierr = 3;
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nn_K; i++) D_inv[i] = 1.0 / coef_K[iat_K[i]];
         nnz_PK = nn_K;
      } else {
         // Preconditioner not available
         cout << "Preconditioner not available" << endl;
         return ierr = 3;
      }
      end = chrono::system_clock::now();
      elaps_sec = end - start;
      time_prec_K = elaps_sec.count();
      if (DUMP_PREC){
         FILE *dfile = fopen("D_mat","w");
         for (int k = 0; k < nn_K; k++) fprintf(dfile,"%20.11e\n",D_inv[k]);
         fflush(dfile);
         fclose(dfile);
      }

      // ----- Copy the initial prolongation into the one with extended pattern
      if (dump) cout << "---- COPY PROL ----" << endl << endl;
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
      if (dump) cout << "---- gather_B_Qr ----" << endl << endl;
      start = chrono::system_clock::now();
      int *pt_Z = nullptr;
      int *pt_col_Z = nullptr;
      double *mat_Q = nullptr;
      double *mat_Z = nullptr;
      double *vec_f = nullptr;
      if (sol_type == DEFL_CG){
         ierr = gather_B_QR(np,condmax,nn,nn_C,ntv,fcnode,iat_patt,ja_patt,TV,mat_Q,coef_P0);
         if (ierr != 0) return ierr = 4;
         nnz_ZQ = ntv*iat_patt[nn];

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
      } else if (sol_type == NULL_SPACE){
         ierr = gather_B_Z(np,nn,nn_C,ntv,fcnode,iat_patt,ja_patt,iat_Tpatt,ja_Tpatt,
                           iat_A,ja_A,coef_A,TV,pt_Z,pt_col_Z,mat_Z,vec_f,coef_P0);
         if (ierr != 0) return ierr = 4;
         nnz_ZQ = pt_Z[nn];
      }
      end = chrono::system_clock::now();
      elaps_sec = end - start;
      time_gath_B = elaps_sec.count();
      if (DUMP_PREC){
         FILE *corrPfile = fopen("corrProl.csr","w");
         wrCSRmat(corrPfile,false,nn,iat_patt,ja_patt,coef_P0);
         fclose(corrPfile);
         if (sol_type == DEFL_CG){
            FILE *Qfile = fopen("mat_Q2","w");
            print_Q(Qfile,nn,ntv,iat_patt,mat_Q);
            fflush(Qfile);
            fclose(Qfile);
         } else if (sol_type == NULL_SPACE){
            FILE *Zfile = fopen("mat_Z","w");
            print_Z(Zfile,nn,ntv,iat_patt,mat_Z);
            fflush(Zfile);
            fclose(Zfile);
         }
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

      // Choose solution strategy
      if (sol_type == DEFL_CG){

         // Compute prolongation correction with PCG
         if (dump) cout << "---- DEFL_PCG ----" << endl << endl;
         start = chrono::system_clock::now();
         // Allocate prolongation correction
         double *DP = (double*) malloc( nt_patt*sizeof(double) );
         if (DP == nullptr) return ierr = 1;
         ierr = DEFL_PCG(np,prec_type,nn,nn_C,nt_patt,ntv,perm,iperm,Tr_A,iat_K,ja_K,
                         coef_K,it_U,jcol_U,coef_U,D_inv,iat_patt,ja_patt,iat_Tpatt,
                         ja_Tpatt,mat_Q,coef_P0,vec_f,itmax,iter,relres,DP);
         if (ierr != 0) return ierr = 5;
         end = chrono::system_clock::now();
         elaps_sec = end - start;
         time_PCG = elaps_sec.count();
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         cout << "PCG TIME "<< time_PCG << endl;
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         // Update prolongation with DP
         #pragma omp parallel for num_threads(np)
         for (int i = 0; i < nt_patt; i++) coef_P0[i] -= DP[i];
         
         // Free memory
         free(DP);

      } else if (sol_type == NULL_SPACE){

         // Compute prolongation with the null-space method
         if (dump) cout << "---- NULLSPACE_PCG ----" << endl << endl;
         start = chrono::system_clock::now();
         ierr = NLSP_PCG(np,prec_type,nn,nn_C,nt_patt,ntv,perm,iperm,Tr_A,iat_K,ja_K,
                         coef_K,it_U,jcol_U,coef_U,D_inv,iat_patt,ja_patt,iat_Tpatt,
                         ja_Tpatt,pt_Z,pt_col_Z,mat_Z,vec_f,itmax,iter,relres,coef_P0);
         if (ierr != 0) return ierr = 5;
         end = chrono::system_clock::now();
         elaps_sec = end - start;
         time_PCG = elaps_sec.count();

      }

      // Free scracthes
      free(pt_Z);
      free(pt_col_Z);
      free(mat_Z);
      free(vec_f);
      free(mat_Q);
      free(it_U);
      free(jcol_U);
      free(coef_U);
      free(D_inv);
      free(iat_K);
      free(ja_K);
      free(coef_K);
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
      info[0]  = time_gath_K;
      info[1]  = time_prec_K;
      info[2]  = time_gath_B;
      info[3]  = time_PCG;
      info[6]  = avg_lfil;
      info[7]  = static_cast<double>(iter);
      info[8]  = relres;
      info[9]  = static_cast<double>(nnz_K);
      info[10] = static_cast<double>(nnz_PK);
      info[11] = static_cast<double>(nnz_ZQ);

   }

   //---GLOBAL END---------------------------------------------------------------
   glob_end = chrono::system_clock::now();
   elaps_sec = glob_end - glob_start;
   double time_glob =  elaps_sec.count();
   double time_overhead = time_glob-time_gath_K-time_prec_K-time_gath_B-time_PCG;
   //----------------------------------------------------------------------------

   // Store info
   info[4]  = time_overhead;
   info[5]  = time_glob;

   return ierr;

}
