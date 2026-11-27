//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
#include <iostream>
#include <iomanip>
#define DBNODE  7
#define ENDNODE 7
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
/*****************************************************************************************
 *
 * Inner part of "EXTI_prolongation" that is used to distribute work among threads.
 * See "EXTI_prolongation" for parameters.
 *
*****************************************************************************************/

#include <stdlib.h>   // to use: NULL
#include <cmath>      // to use: abs
#include <algorithm>  // to use: fill_n,min
#include <limits>     // to use: numeric_limits<double>epsilon()
#include <new>
#include <vector>

#include "ir_heapsort.h"

const double ONE  = 1.0;
const double ZERO = 0.0;
const double EPS = std::numeric_limits<double>::epsilon();

int ONELP_ProlStripe_EXTI(const int firstrow, const int lastrow, const int nn_loc,
                          const int nn_A, const int nt_A,
                          const int *const iat_A, const int *const ja_A,
                          const double *const coef_A, const int *const coef_S,
                          const int *const fcnodes,
                          int &nn_P, int &nt_P, int *iat_P,
                          std::vector<int> &ja_P, std::vector<double> &coef_P){

   // Early exit if stripe is empty
   if (firstrow >= lastrow || nn_loc <= 0) {
      nn_P = 0;
      nt_P = 0;
      iat_P[0] = 0;
      return 0;
   }

   // Allocate scratch vectors based on global nn_A to avoid division by zero or thread-count blowup
   double avg_nnz = (nn_A > 0) ? (20.0 * static_cast<double>(nt_A) / static_cast<double>(nn_A)) : 50.0;
   int size_scr = std::max(64, static_cast<int>(avg_nnz));
   int *WI = new (std::nothrow) int [nn_A]();
   if (WI == nullptr) {
      return 1;
   }
   std::vector<int> ja_FS;
   std::vector<double> coef_FS;
   std::vector<int> list_weak;
   ja_FS.reserve(size_scr);
   coef_FS.reserve(size_scr);
   list_weak.reserve(size_scr);

   // Init pointer to the prolongation matrix
   nn_P = 0;
   int ind_P = 0;
   iat_P[0] = ind_P;

   // Loop over all the nodes
   int shift = firstrow - 1;
   for (int inod = firstrow; inod < lastrow ; inod ++){ // Node_loop

      int inod_coarse = fcnodes[inod];
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //cout << inod << " " << inod_coarse << " " << ind_P << endl;
      //if (inod == DBNODE) cout << "NODO: " << ((inod_coarse>=0) ? "CO":"FI") << endl;
      //if (inod == DBNODE){
      //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){
      //      cout << "jcol_A " << ja_A[j] << " coef_A " << coef_A[j] << " coef_S " << 
      //              coef_S[j] << " fcnode " << fcnodes[ja_A[j]] << endl;
      //   }
      //}
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      // Check whether inod is FINE OR COARSE
      if (inod_coarse >= 0){
         nn_P++;

         // It is a coarse node
         ja_P.push_back(inod_coarse);
         coef_P.push_back(ONE);
         ind_P++;

      } else if (inod_coarse == -1) {
         nn_P++;

         // NOTE: WI interpolatory set marker >  0 --> interpolatory node
         //                                   == 0 --> not assigned yet
         //                                   <  0 --> weak coarse node

         // First loop to classify distance 1 coarse neighbours of inod
         int n_int = 0;
         list_weak.clear();
         double denom = ZERO;
         double a_ii = ZERO;
         for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){ // Loop on dist_1 neigh
            int jcol = ja_A[j];
            if (jcol == inod){
               // Add the diagonal entry to denominator
               a_ii = coef_A[j];
               denom += a_ii;
            } else {
               if (fcnodes[jcol] >= 0){
                  // This neighbour is COARSE
                  if (coef_S[j] > 0){
                     // Strong connection

                     // Add this node to the interp. set and save its value
                     ja_P.push_back(jcol);
                     coef_P.push_back(coef_A[j]);
                     WI[jcol] = n_int+1;
                     n_int++;
                  } else {
                     // Weak connection

                     // The node is not in the interp. set (not yet)
                     denom += coef_A[j];
                     // Keep track of weak COARSE neighbours in case they are promoted
                     // to interpolatory points
                     WI[jcol] = -(j + 1);
                     list_weak.push_back(jcol);
                  }
               } else {
                  // This neighbour is FINE
                  if (coef_S[j] < 0){
                     // Weak connection, just adjust denominator
                     denom += coef_A[j];
                  } 
               }
            }
         } // End loop on dist_1 neig
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22
         //if (inod == DBNODE){
         //   cout << endl << endl;
         //   cout << "INTERP SET:" << endl;
         //   for (int i = 0; i < n_int; i++) cout << ja_P[ind_P+i] << " ";
         //   cout << endl;
         //   cout << "WI INTERP SET:" << endl;
         //   for (int i = 0; i < n_int; i++) cout << WI[ja_P[ind_P+i]] << " ";
         //   cout << endl;
         //   cout << "DIST 1 - NEIGH:" << endl;
         //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++)
         //      if (WI[ja_A[j]] > 0) cout << ja_A[j] << " ";
         //   cout << endl;
         //   cout << "WI DIST 1 - NEIGH:" << endl;
         //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++)
         //      if (WI[ja_A[j]] > 0) cout << WI[ja_A[j]] << " ";
         //   cout << endl;
         //   cout << "----------------------------------" << endl << endl;
            //double NUM = 0.0;
            //double DEN = 0.0;
            //for (int i = 0; i < n_int; i++) NUM+= coef_P[ind_P+i];
            //for (int j = iat_C[inod]; j < iat_C[inod+1]; j++){
                //if (fcnodes[ja_C[j]] < 0) DEN += coef_C[j];
            //}
            //cout << "SOMMA NUM: " << NUM << endl;
            //cout << "SOMMA DEN: " << DEN + DDD << endl;
            //cout << "----------------------------------" << endl << endl;
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22

         // Explore strong FINE neighbors directly using (iat_A, ja_A, coef_A) guarded by coef_S > 0
         for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){ // Loop on dist_1 neigh
            if (coef_S[j] <= 0) continue;
            int jcol = ja_A[j];
            if (jcol != inod){
               // Consider only FINE nodes
               if (fcnodes[jcol] < 0){
                  double a_ik = coef_A[j];
                  double a_ki_bar = ZERO;
                  double ext_sum = 0.0;
                  ja_FS.clear();
                  coef_FS.clear();
                  // Explore its strong neighbours to expand the interpolatory set
                  for (int k = iat_A[jcol]; k < iat_A[jcol+1]; k++){ // Strong neigh loop
                     if (coef_S[k] <= 0) continue;
                     int kcol = ja_A[k];
                     if (kcol == inod){
                        // This is the current FINE node --> save its coefficient
                        a_ki_bar = std::min(ZERO,coef_A[k]);
                        // Update the extended sum
                        ext_sum += a_ki_bar;
                     } else {
                        // Consider only COARSE neighbors
                        if (fcnodes[kcol] >= 0){
                           // Update the extended sum and store an auxiliary coef
                           double a_kj_bar = std::min(ZERO,coef_A[k]);
                           ja_FS.push_back(kcol);
                           coef_FS.push_back(a_kj_bar);
                           ext_sum += a_kj_bar;
                           // Check that the node has not been added yet to the interp. set
                           if ( WI[kcol] <= 0 ){
                              double value = ZERO;
                              if (WI[kcol] < 0){
                                 // The node was considered weak COARSE before
                                 value = coef_A[-WI[kcol] - 1];
                                 // Adjust denominator when the node is included
                                 denom -= value;
                              }
                              // kcol was not in the interp. set, add it and init its
                              // coefficient
                              ja_P.push_back(kcol);
                              coef_P.push_back(value);
                              WI[kcol] = n_int+1;
                              n_int++;
                           }
                        }
                     }
                  } // End Strong neigh loop

                  // Check that extended sum is not null
                  if ( std::abs(ext_sum) > EPS*std::abs(a_ii) ){
                     // ext_sum !=0

                     // Update weights
                     for (std::size_t k = 0; k < ja_FS.size(); k++){
                        int pos = WI[ja_FS[k]]-1;
                        coef_P[ind_P+pos] += (a_ik*coef_FS[k]) / ext_sum;
                     }
                     // Update denominator
                     denom += (a_ik * a_ki_bar) / ext_sum;

                  } else {
                     // ext_sum == 0
                  
                     // Treat the case of zero extended sum
                     double fac = ONE / (static_cast<double>(ja_FS.size() + 1));
                     // Update denominator
                     denom += a_ik*fac;
                     // Update all the weights
                     for (std::size_t k = 0; k < ja_FS.size(); k++){
                        int pos = WI[ja_FS[k]]-1;
                        coef_P[ind_P+pos] += a_ik*fac;
                     }

                  }
               }
            }
         } // End loop on dist_1 neigh
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         //if (inod == DBNODE){
         //   cout << inod << " DENOM FINALE " << denom << endl;
         //   cout << "INTERP SET:" << endl;
         //   for (int i = ind_P; i < ind_P+n_int; i++) cout << ja_P[i] << " ";
         //      cout << endl;
         //   cout << "COEF DENOM " << endl;
         //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){
         //      if (coef_S[j] < 0){
         //         int jnod = ja_A[j];
         //         cout << jnod << " " << coef_A[j] << " " << WI[jnod] << " " << (WI[jnod]>0 ? "INT":"NOI") << endl;
         //      }
         //   }
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         // Check for zero or near-zero denominator
         if (std::abs(denom) < EPS) {
            denom = (std::abs(a_ii) > EPS) ? a_ii : ONE;
         }

         // Adjust column index in ja_P, compute weights and sort
         for (int i = ind_P; i < ind_P+n_int; i++){
            // Reset WI
            WI[ja_P[i]] = 0;
            ja_P[i] = fcnodes[ja_P[i]];
            coef_P[i] = -coef_P[i] / denom;
         }
         if (n_int > 0) {
            ir_heapsort(&(ja_P[ind_P]),&(coef_P[ind_P]),n_int);
         }
         ind_P += n_int;

         // Reset WI for weak coarse nodes
         for (std::size_t i = 0; i < list_weak.size(); i++) WI[list_weak[i]] = 0;

      }

      // Update pointer to next row
      iat_P[nn_P] = ind_P;

   } // Node_loop

   // Count number of non-zeroes
   nt_P = ind_P;

   // Deallocate scratches
   delete [] WI;

   return 0;
}
