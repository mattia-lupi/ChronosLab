//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
#include <iostream>
#include <iomanip>
#define DBNODE  000
#define ENDNODE 200
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
using namespace std;

#include "ir_heapsort.h"

const double ONE  = 1.0;
const double ZERO = 0.0;
const double EPS = numeric_limits<double>::epsilon();

int TWOLP_ProlStripe_EXTI(const int firstrow, const int lastrow, const int nn_loc,
                          const int nn_A, const int nt_A, const int ntmax_P,
                          const int *const iat_A, const int *const ja_A,
                          const double *const coef_A, const int *const coef_S,
                          const int *const iat_C, const int *const ja_C,
                          const double *const coef_C, const int *const fcnodes,
                          int &nn_P, int &nt_P, int *iat_P, int *ja_P,
                          double *coef_P){

   // Allocate some scratch vectors
   // PER ESSERE SICURI BISOGNA AVERE IL MAX NUMERO DI NONZERI PER RIGA
   double avg_nnz = 20.0*static_cast<double>(nt_A) / static_cast<double>(lastrow-firstrow);
   int size_scr = static_cast<int>(avg_nnz);
   int *WI         = new int [nn_A](); if (WI == NULL) return 1;
   fill_n(WI,nn_A,0);
   int *ja_FS      = new int [size_scr](); if (ja_FS == NULL) return 1;
   double *coef_FS = new double [size_scr](); if (coef_FS == NULL) return 1;
   int *list_weak  = new int [size_scr](); if (list_weak == NULL) return 1;

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
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      // Check whether inod is FINE OR COARSE
      if (inod_coarse >= 0){
         nn_P++;

         // It is a coarse node
         ja_P[ind_P] = inod_coarse;
         coef_P[ind_P] = ONE;
         ind_P++;

      } else if (inod_coarse == -1) {
         nn_P++;

         // NOTE: WI interpolatory set marker >  0 --> interpolatory node
         //                                   == 0 --> not assigned yet
         //                                   <  0 --> weak coarse node

         // First loop to classify distance 1 coarse neighbours of inod
         int n_int = 0;
         int n_weak = 0;
         double denom = ZERO;
         for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){ // Loop on dist_1 neigh
            int jcol = ja_A[j];
            if (jcol == inod){
               // Add the diagonal entry to denominator
               double a_ii = coef_A[j];
               denom += a_ii;
            } else {
               if (fcnodes[jcol] >= 0){
                  // This neighbour is COARSE
                  if (coef_S[j] > 0){
                     // Strong connection

                     if (WI[jcol] <= 0){
                        // Add this node to the interp. set and save its value
                        ja_P[ind_P+n_int] = jcol;
                        coef_P[ind_P+n_int] = coef_A[j];
                        WI[jcol] = n_int+1;
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        //if (inod == DBNODE){
                        //   cout << "XXXXXXXXXAGGIUNTO CC: " << jcol << " " << WI[jcol] << " " << coef_A[j] << endl;
                        //}
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        n_int++;
                     }
                  } else {
                     // Weak connection

                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     //if (inod == DBNODE){
                     //   cout << "*********WEAK COARSE 1  " << jcol << " " << WI[jcol] << endl;
                     //}
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     // The node is not in the interp. set (not yet)
                     denom += coef_A[j];
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     //if (inod == DBNODE){
                     //    cout << "+++++++AGGIUNGO XXX " << jcol << " " << coef_A[j] << endl;
                     //}
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     // Keep track of weak COARSE neighbours in case they are promoted
                     // to interpolatory points
                     WI[jcol] = -1;
                     list_weak[n_weak] = jcol;
                     n_weak++;
                  }
               } else {
                  // This neighbour is FINE
                  if (coef_S[j] < 0){
                     // Weak connection, just adjust denominator
                     denom += coef_A[j];
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     //if (inod == DBNODE){
                     //    cout << "+++++++AGGIUNGO XXX " << jcol << " " << coef_A[j] << endl;
                     //}
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
         //   //double NUM = 0.0;
         //   //double DEN = XXX;
         //   //for (int i = 0; i < n_int; i++) NUM+= coef_P[ind_P+i];
         //   //for (int j = iat_C[inod]; j < iat_C[inod+1]; j++){
         //   //    if (fcnodes[ja_C[j]] < 0) DEN += coef_C[j];
         //   //}
         //   //cout << "SOMMA NUM: " << NUM << endl;
         //   //cout << "SOMMA DEN: " << DEN + DDD << endl;
         //   //cout << "----------------------------------" << endl << endl;
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22

         // Explore strong FINE neighbors to extend the interpolatory set (using
         // the compressed matrix C)
         for (int j = iat_C[inod]; j < iat_C[inod+1]; j++){ // Loop on dist_1 neigh
            int jcol = ja_C[j];
            if (jcol != inod){ //This check is necessary only if C has diagonals
               // Consider only FINE nodes
               if (fcnodes[jcol] < 0){
                  double a_ik = coef_C[j];
                  double a_ik_bar = min(ZERO,a_ik);
                  double ext_sum = a_ik_bar;
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  //if (inod == DBNODE){
                  //   cout << "FNEIGH " << jcol << " INIT ESUM " << a_ik_bar << endl;
                  //}
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  int n_FS = 0;
                  // Explore its strong neighbours to expand the interpolatory set
                  for (int k = iat_C[jcol]; k < iat_C[jcol+1]; k++){
                     int kcol = ja_C[k];
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     //if (inod == DBNODE){
                     //   cout << ja_C[k] << "  " << ((fcnodes[kcol] >= 0) ? "CO":"FI") << endl;
                     //}
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     // Consider only COARSE neighbors
                     if (fcnodes[kcol] >= 0){
                        // Update the extended sum and store an auxiliary coef
                        double a_kj_bar = min(ZERO,coef_C[k]);
                        ja_FS[n_FS] = kcol;
                        coef_FS[n_FS] = a_kj_bar;
                        n_FS++;
                        ext_sum += a_kj_bar;
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        //if (inod == DBNODE){
                        //   cout << "ESUM ADD1 " << kcol << " VAL " << a_kj_bar << endl;
                        //}
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                        // Check that the node has not been added yet to the interp. set
                        if ( WI[kcol] <= 0 ){
                           //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                           //if (inod == DBNODE){
                           //   cout << "+++++++AGGIUNGO WC INT " << kcol << " " << WI[kcol] << endl;
                           //   if (WI[kcol] < 0) cout << "+++++++TOLGO DEN " << " " << kcol << coef_A[-WI[kcol]] << endl;
                           //}
                           //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                           // kcol was not in the interp. set, add it and init its
                           // coefficient
                           ja_P[ind_P+n_int] = kcol;
                           coef_P[ind_P+n_int] = ZERO;
                           WI[kcol] = n_int+1;
                           n_int++;
                        }
                     }
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  //if (inod == DBNODE){
                  //   cout << "EXT_SUM FINAL" << ext_sum << endl;
                  //   cout << "WGT CORR:  " << endl;
                  //   for (int k = 0; k < n_FS; k++){
                  //      int pos = WI[ja_FS[k]];
                  //      pos = (pos>0) ? (pos-1):(-pos-2);
                  //      double sav = coef_FS[k];
                  //      cout << pos << " " << ja_FS[k] << " " << (a_ik*coef_FS[k]) << endl;
                  //   }
                  //}
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  // Update weights
                  for (int k = 0; k < n_FS; k++){
                     int pos = WI[ja_FS[k]]-1;
                     coef_P[ind_P+pos] += (a_ik*coef_FS[k]) / ext_sum;
                  }
                  // Update denominator
                  denom += (a_ik * a_ik_bar) / ext_sum;
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

         // Adjust column index in ja_P, compute weights and sort
         for (int i = ind_P; i < ind_P+n_int; i++){
            // Reset WI
            WI[ja_P[i]] = 0;
            ja_P[i] = fcnodes[ja_P[i]];
            coef_P[i] = -coef_P[i] / denom;
         }
         ir_heapsort(&(ja_P[ind_P]),&(coef_P[ind_P]),n_int);
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22
         //double sum = 0.0;
         //for (int i = ind_P; i < ind_P+n_int; i++) sum += coef_P[i];
         //for (int i = ind_P; i < ind_P+n_int; i++) coef_P[i] /= sum;
         //if (inod >= DBNODE && inod <=ENDNODE){
            //cout << endl;
            //cout << "WGTS:" << endl;
            //for (int i = ind_P; i < ind_P+n_int; i++) cout << ja_P[i] << " ";
            //cout << endl;
            //for (int i = ind_P; i < ind_P+n_int; i++) cout << coef_P[i] << " ";
            //cout << endl;
            //double sum = 0.0;
            //for (int i = ind_P; i < ind_P+n_int; i++){
            //   sum += coef_P[i];
            //}
            //cout << "SOMMA " << sum << endl;
            //cout << "XXXX " << XXX << endl;
         //}
         //if (inod >= DBNODE && inod <= ENDNODE){
         //   double sum_P = 0.0;
         //   for (int i = ind_P; i < ind_P+n_int; i++){
         //      sum_P += coef_P[i];
         //   }
         //   double sum_A = 0.0;
         //   for (int i = iat_A[inod]; i < iat_A[inod+1]; i++){
         //      sum_A += coef_A[i];
         //   }
         //   cout << "SOMMA A | SOMMA P: " << setw(12) << sum_A/DDD << "  " << sum_P << endl;
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@22
         ind_P += n_int;

         // Reset WI for weak coarse nodes
         for (int i = 0; i < n_weak; i++) WI[list_weak[i]] = 0;

      }

      // Update pointer to next row
      iat_P[nn_P] = ind_P;

   } // Node_loop

   // Count number of non-zeroes
   nt_P = ind_P;

   // Deallocate scratches
   delete [] WI;
   delete [] ja_FS;
   delete [] coef_FS;
   delete [] list_weak;

   return 0;
}
