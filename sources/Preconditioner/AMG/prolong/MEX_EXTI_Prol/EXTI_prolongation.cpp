/*****************************************************************************************
 *
 * This function computes the Extended + I interpolation that is very effective on large
 * scale Poisson problems in conjuction with PMIS coarsening (see Hans De Sterck et al.,
 * Distance-two interpolation for parallel algebraic multigrid, Numerical Linear Algebra
 * with Applications 15, pp. 115–139, 2007).
 *
 * Input parameters:
 *
 * level:     AMG level indicator
 * nthreads:  number of openMP threads used
 * vecstart:  pointer to the beginning of the region owned by each thread
 * nn_A:      number of rows of the A matrix
 * nt_A:      number of entries of the A matrix
 * iat_A:     array of pointers of the A matrix
 * ja_A:      array of column indices of the A matrix
 * coef_A:    array of coefficients of the A matrix
 * coef_S:    array of coefficients of the SoC matrix (it has the same pattern as A)
 * iat_C:     array of pointers of the compressed A matrix
 * ja_C:      array of column indices of the compressed A matrix
 * coef_C:    array of coefficients of the compressed A matrix
 * fcnodes:   F/C indicator:    fcnodes <  0 ==> FINE node
 *                              fcnodes >= 0 ==> COARSE node with number fcnodes
 * nr_I:      Number of rows of the interpolation matrix
 * nc_I:      Number of columns of the interpolation matrix
 *
 * Output parameters:
 *
 * nt_I:      Number of entries of the interpolation matrix
 * iat_I:     array of pointers of the interpolation matrix
 * ja_I:      array of column indices of the interpolation matrix
 * coef_I:    array of coefficients of the interpolation matrix
 *
 * Return parameter:
 *
 * ierr:      Error code: TBD
 *
*****************************************************************************************/

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#include <chrono>
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#include <stdlib.h>  // to use: exit
#include <iostream>
#include <new>
#include <vector>
#include "omp.h"
using namespace std;

#include "MIS.h"
#include "ProlStripe_EXTI.h"
#include "ONELP_ProlStripe_EXTI.h"
#include "TWOLP_ProlStripe_EXTI.h"

int EXTI_prolongation(const int level, const int nthreads,
                      const int nn_A, const int nt_A,
                      const int *const iat_A, const int *const ja_A,
                      const double *const coef_A, const int *const coef_S,
                      const int *const fcnodes,
                      const int nr_I, const int nc_I, int &nt_I, int *&iat_I,
                      int *&ja_I, double *&coef_I){

   // Compute thread row partitions in C++
   std::vector<int> vecstart(nthreads + 1, 0);
   if (nthreads > 0 && nn_A > 0) {
      int bsize = nn_A / nthreads;
      int resto = nn_A % nthreads;
      for (int i = 1; i <= resto; ++i) vecstart[i] = vecstart[i - 1] + bsize + 1;
      for (int i = resto + 1; i <= nthreads; ++i) vecstart[i] = vecstart[i - 1] + bsize;
   }

   // Internal scratch variables
   int *ridv_i;
   int *ridvn_i;

   // Init error code
   int ierr = 0;

   // Allocate some scratch vectors
   ridv_i = new int [nthreads+1]();
   ridvn_i = new int [nthreads+1]();
   if (ridv_i == nullptr || ridvn_i == nullptr) {
      delete [] ridv_i;
      delete [] ridvn_i;
      return 1;
   }

   int nn_I;
   // Compute the coefficients of the Prolongation matrix
   //-------------------------------------------------------------------------------------
   #pragma omp parallel num_threads(nthreads)
   //-------------------------------------------------------------------------------------
   {
      std::chrono::time_point<std::chrono::high_resolution_clock> START_OLD,END_OLD;
      std::chrono::time_point<std::chrono::high_resolution_clock> START_ONE,END_ONE;
      std::chrono::time_point<std::chrono::high_resolution_clock> START_TWO,END_TWO;
      std::chrono::duration<double> time_OLD,time_ONE,time_TWO;
      // Declaration of some private variables
      int istart_I, shift, istart_scr, iend_scr;

      // Init local error code
      int ierr_L = 0;

      // Get the processor ID
      int myid = omp_get_thread_num();

      // Get first and last rows of the current thread
      int firstrow = vecstart[myid];
      int lastrow  = vecstart[myid + 1];
      int nn_loc = lastrow - firstrow;
      int nt_I_loc;
      int nn_I_loc;
      
      // Allocate local scratches
      int *iat_scr = new (std::nothrow) int [nn_loc+1]();
      std::vector<int> ja_scr;
      std::vector<double> coef_scr;
      if (nn_loc > 0) {
         ja_scr.reserve(static_cast<size_t>(16) * nn_loc);
         coef_scr.reserve(static_cast<size_t>(16) * nn_loc);
      }
      if (iat_scr == nullptr) {
         ierr_L = 1;
      }

      // Compute the local stripe of the prolongation
      if (ierr_L == 0) {
         ierr_L = ONELP_ProlStripe_EXTI(firstrow,lastrow,nn_loc,nn_A,nt_A,
                                        iat_A,ja_A,coef_A,coef_S,
                                        fcnodes,nn_I_loc,nt_I_loc,iat_scr,ja_scr,
                                        coef_scr);
      } else {
         nn_I_loc = 0;
         nt_I_loc = 0;
      }

      #pragma omp atomic update
      ierr += ierr_L;

      // Store local number of entries
      ridv_i[myid+1] = nt_I_loc;
      ridvn_i[myid+1] = nn_I_loc;
      #pragma omp barrier

      if (ierr == 0) {
         // Allocate the prolongation matrix
         #pragma omp single
         {
            // Count the number of non-zeroes
            ridv_i[0] = 0;
            for (int ip = 0; ip < nthreads; ip++){
                ridv_i[ip+1] = ridv_i[ip] + ridv_i[ip+1];
            }
            nt_I = ridv_i[nthreads];

            // Count the number of prolongation rows
            ridvn_i[0] = 0;
            for (int ip = 0; ip < nthreads; ip++){
                ridvn_i[ip+1] = ridvn_i[ip] + ridvn_i[ip+1];
            }
            nn_I = ridvn_i[nthreads];
            if (nn_I != nr_I) ierr += 1;

            iat_I = new (std::nothrow) int [nn_I+1](); if (iat_I == nullptr) ierr += 1;
            ja_I = new (std::nothrow) int [nt_I](); if (ja_I == nullptr) ierr += 1;
            coef_I = new (std::nothrow) double [nt_I](); if (coef_I == nullptr) ierr += 1;
         }
         #pragma omp barrier

         if (ierr == 0) {
            // Store iat_scr, ja_scr and coef_scr in the prolongation matrix
            istart_I = ridv_i[myid];
            shift = ridvn_i[myid];
            iend_scr = iat_scr[0];
            for (int irow = 0; irow < nn_I_loc; irow++){
               iat_I[shift+irow] = istart_I;
               istart_scr = iend_scr;
               iend_scr = iat_scr[irow+1];
               int len_I = iend_scr - istart_scr;
               for (int k = 0; k < len_I; k++){
                  ja_I[istart_I+k] = ja_scr[istart_scr+k];
                  coef_I[istart_I+k] = coef_scr[istart_scr+k];
               }
               istart_I += len_I;
            }
            if (myid == nthreads-1) iat_I[nn_I] = ridv_i[nthreads];
         }
      }

      // Deallocate local scratches safely on all paths
      delete [] iat_scr;

   }

   delete [] ridv_i;
   delete [] ridvn_i;

   return ierr;
}
