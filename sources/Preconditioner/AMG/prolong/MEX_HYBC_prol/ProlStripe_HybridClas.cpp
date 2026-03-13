/*****************************************************************************************
 *
 * Inner part of "HybridClas_prolongation" that is used to distribute work among threads.
 * See "HybridClas_prolongation" for parameters.
 *
*****************************************************************************************/

#include <stdlib.h>   // to use: NULL
#include <cmath>      // to use: abs
#include <algorithm>  // to use: fill_n,min
#include <limits>     // to use: numeric_limits<double>epsilon()
using namespace std;

#include "ir_heapsort.h"
#include "find_bad_FS.h"
#include "find_dist2CC.h"

const double ONE  = 1.0;
const double ZERO = 0.0;
const double EPS = numeric_limits<double>::epsilon();

int ProlStripe_HybridClas(const int firstrow, const int lastrow,
                          const int nn_loc, const int nn_A, const int nt_A,
                          const int ntmax_P, const int *const iat_A, const int *const ja_A,
                          const double *const coef_A, const int *const coef_S,
                          const int *const fcnodes, int &nt_P, int *iat_P, int *ja_P,
                          double *coef_P){

   // Allocate some scratch vectors
   double avg_nnz = 1.2*static_cast<double>(nt_A) / static_cast<double>(lastrow-firstrow);
   int size_scr = static_cast<int>(avg_nnz)*static_cast<int>(avg_nnz);
   int *WI         = new int [nn_A](); if (WI == NULL) return 1;
   fill_n(WI,nn_A,0);
   int *ja_CC      = new int [size_scr](); if (ja_CC == NULL) return 1;
   int *ja_FS      = new int [size_scr](); if (ja_FS == NULL) return 1;
   int *pos_kj     = new int [size_scr](); if (pos_kj == NULL) return 1;
   double *coef_CC = new double [size_scr](); if (coef_CC == NULL) return 1;
   double *coef_FS = new double [size_scr](); if (coef_FS == NULL) return 1;
   double *a_kj    = new double [size_scr](); if (a_kj == NULL) return 1;

   int *ind_bad_FS = new int [size_scr](); if (ind_bad_FS == NULL) return 1;
   int *ja_CC_2    = new int [size_scr](); if (ja_CC_2 == NULL) return 1;
   double *coef_CC_2 = new double [size_scr](); if (coef_CC_2 == NULL) return 1;
   int *deg_CC_2   = new int [nn_A](); if (deg_CC_2 == NULL) return 1;
   fill_n(WI,nn_A,0);

   // Init pointer to the prolongation matrix
   int ind_P = 0;
   iat_P[0] = ind_P;

   // Loop over all the nodes
   int shift = firstrow - 1;
   for (int inod = firstrow; inod < lastrow ; inod ++){ // Node_loop

      int inod_coarse = fcnodes[inod];

      // Check whether inod is FINE OR COARSE
      if (inod_coarse >= 0){

         // It is a coarse node
         ja_P[ind_P] = inod_coarse;
         coef_P[ind_P] = ONE;
         ind_P++;

      } else {

         // Explore and classify distance 1 neighbours of inod
         int n_FS = 0;
         int n_CC = 0;
         double a_ii;
         double denom = ZERO;
         for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){ // Loop on dist_1 neigh
            int jcol = ja_A[j];
            if (jcol == inod){
               // Add the diagonal entry to denominator
               a_ii = coef_A[j];
               denom += a_ii;
            } else {
               if (coef_S[j] > 0){
                  // Strong connection
                  int jnod = fcnodes[jcol];
                  if (jnod >= 0){
                     // This neighbour is COARSE
                     // This node was not in the interpolatory set
                     ja_CC[n_CC] = jcol;
                     coef_CC[n_CC] = coef_A[j];
                     WI[jcol] = n_CC+1;
                     n_CC++;
                  } else {
                     // This neighbour is FINE --> store coefficient and index in an
                     // auxiliary list
                     ja_FS[n_FS] = jcol;
                     coef_FS[n_FS] = coef_A[j];
                     n_FS++;
                     // Explore its strong neighbours to expand the interpolatory set
                  }
               } else {
                  // Weak connection increase the denominator with the matrix entry
                  // if it is not yet in the interpolatory set
                  denom += coef_A[j];
               }
            }
         } // End of loop on dist_1 neigh

         // Find bad fine nodes
         int n_bad_FS;
         find_bad_FS(n_FS,ja_FS,WI,iat_A,ja_A,n_bad_FS,ind_bad_FS);
         
         // Find distance two coarse node
         while (n_bad_FS>0){
         // Find distance two coarse node
           int n_CC_2;
           find_dist2CC(n_bad_FS,ind_bad_FS,ja_FS,iat_A,ja_A,coef_A,coef_S,fcnodes, 
                        n_CC_2, ja_CC_2, coef_CC_2,deg_CC_2);
            // Add node with maximum degree
            int ind_max=0;
            int max_deg=0;
            for(int l=0; l<n_CC_2; l++){
               int lcol=ja_CC_2[l];
               if (deg_CC_2[lcol]>max_deg){
                  ind_max=l;
                  max_deg=deg_CC_2[lcol];
               } 
            }
            ja_CC[n_CC]=ja_CC_2[ind_max];
            coef_CC[n_CC]=0.0;
            WI[ja_CC_2[ind_max]]=n_CC+1;   
            n_CC++;

            // Reset the deg_CC_2
            for (int k=0; k<n_CC_2; k++){
               int kcol=ja_CC_2[k];
               deg_CC_2[kcol]=0;
            }   

            find_bad_FS(n_FS,ja_FS,WI,iat_A,ja_A,n_bad_FS,ind_bad_FS);
            
         }

         // Init weights of interpolatory points and sparse reset of WI relative to
         // weak COARSE neighbours
         int n_int = 0;
         for (int i = 0; i < n_CC; i++){
            int jcol = ja_CC[i];
            // This is an interpolatory point
            ja_P[ind_P+n_int] = jcol;
            coef_P[ind_P+n_int] = coef_CC[i];
            WI[jcol] = n_int+1;
            n_int++;
         }
        
         // Correct the interpolation weight to take into account strong FINE neighs
         // Loop over strong FINE neighbors
         for (int k = 0; k < n_FS; k++){
            int knod = ja_FS[k];
            double a_ik = coef_FS[k];
            // Init the  sum 
            double ext_sum = ZERO;
            // Explore closest neighbors of knod considering only nodes in the
            // interpolatory set
            int ind = 0;
            for (int l = iat_A[knod]; l < iat_A[knod+1]; l++){
               int lcol = ja_A[l];
               if ( WI[lcol] > 0){
                  // This is a node in the interpolatory set
                  double a_kl = min(ZERO,coef_A[l]);
                  ext_sum += a_kl;
                  a_kj[ind] = a_kl;
                  pos_kj[ind] = WI[lcol] - 1;
                  ind++;
               }
            }
            // Check that extended sum is not null
            if ( abs(ext_sum) > EPS*a_ii ){
               // Update all the weights
               for (int jj = 0; jj < ind; jj++){
                  int jpos = pos_kj[jj];
                  coef_P[ind_P+jpos] += a_ik*a_kj[jj] / ext_sum;
               }
            } else {
               // Update denominator
               denom += a_ik;
            }
         }

         // Adjust column index in ja_P and sort
         for (int i = ind_P; i < ind_P+n_int; i++){
            // Reset of WI
            WI[ja_P[i]] = 0;
            ja_P[i] = fcnodes[ja_P[i]];
            coef_P[i] = -coef_P[i] / denom;
         }
         ir_heapsort(&(ja_P[ind_P]),&(coef_P[ind_P]),n_int);
         ind_P += n_int;

      }

      // Update pointer to next row
      iat_P[inod-shift] = ind_P;

   } // Node_loop

   // Count number of non-zeroes
   nt_P = ind_P;

   // Deallocate scratches
   delete [] WI;
   delete [] ja_CC;
   delete [] ja_FS;
   delete [] pos_kj;
   delete [] coef_CC;
   delete [] coef_FS;
   delete [] a_kj;

   delete [] ind_bad_FS;
   delete [] ja_CC_2;  
   delete [] coef_CC_2;
   delete [] deg_CC_2;

   return 0;
}
