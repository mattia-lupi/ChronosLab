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

#include "ir_heapsort.h"

const double ONE  = 1.0;
const double ZERO = 0.0;
const double EPS = std::numeric_limits<double>::epsilon();

int ProlStripe_EXTI(const int firstrow, const int lastrow, const int nn_loc,
                    const int nn_A, const int nt_A, const int ntmax_P,
                    const int *const iat_A, const int *const ja_A,
                    const double *const coef_A, const int *const coef_S,
                    const int *const fcnodes, int &nn_P, int &nt_P, int *iat_P, int *ja_P,
                    double *coef_P){

   // Allocate some scratch vectors
   double avg_nnz = 1.2*static_cast<double>(nt_A) / static_cast<double>(lastrow-firstrow);
   int size_scr = static_cast<int>(avg_nnz)*static_cast<int>(avg_nnz);
   int *WI         = new int [nn_A](); if (WI == NULL) return 1;
   std::fill_n(WI,nn_A,0);
   int *ja_CC      = new int [size_scr](); if (ja_CC == NULL) return 1;
   int *ja_FS      = new int [size_scr](); if (ja_FS == NULL) return 1;
   int *pos_kj     = new int [size_scr](); if (pos_kj == NULL) return 1;
   double *coef_CC = new double [size_scr](); if (coef_CC == NULL) return 1;
   double *coef_FS = new double [size_scr](); if (coef_FS == NULL) return 1;
   double *a_kj    = new double [size_scr](); if (a_kj == NULL) return 1;

   // Init pointer to the prolongation matrix
   nn_P=0;
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
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         //if (inod == DBNODE){
         //    for (int i = 0; i < nn_A; i++) if (WI[i] != 0) cout <<"MERDONE!!!!"<< endl;
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
                     if (WI[jcol] > 0){
                        // The node has already been added as second level, just correct
                        // the weight
                        int pos = WI[jcol]-1;
                        coef_CC[pos] = coef_A[j];
                     } else if (WI[jcol] == 0){
                        // This node was not in the interpolatory set
                        ja_CC[n_CC] = jcol;
                        coef_CC[n_CC] = coef_A[j];
                        WI[jcol] = n_CC+1;
                        n_CC++;
                     } else {
                        // This node was considered weak before, now it is promoted in
                        // the interpolatory set
                        WI[jcol] = -WI[jcol];
                     }
                  } else {
                     // This neighbour is FINE --> store coefficient and index in an
                     // auxiliary list
                     ja_FS[n_FS] = jcol;
                     coef_FS[n_FS] = coef_A[j];
                     n_FS++;
                     // Explore its strong neighbours to expand the interpolatory set
                     for (int k = iat_A[jcol]; k < iat_A[jcol+1]; k++){
                        if (coef_S[k] > 0){
                           int kcol = ja_A[k];
                           if (fcnodes[kcol] >= 0){
                              // This node has to be added to the interpolatory set
                              if (WI[kcol] == 0){
                                 // First time the node is visited
                                 ja_CC[n_CC] = kcol;
                                 coef_CC[n_CC] = 0.0;
                                 // Mark the position + 2 of node of second level in order
                                 // to recognize them if encountered in the first level
                                 WI[kcol] = n_CC+1;
                                 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                 //if (inod == DBNODE){
                                 //   cout << "AGGIUNTO dist 2 " << kcol <<  endl;
                                 //}
                                 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                 n_CC++;
                              } else if (WI[kcol] < 0){
                                 // Correct denominator if the node was considered weak before
                                 int pos = -(WI[kcol]+1);
                                 // //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                 //if (inod == DBNODE) cout << "AGGIUSTO DEN " << endl;
                                 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                 denom -= coef_CC[pos];
                                 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                 //if (inod == DBNODE){
                                 //   cout << "*******TOLGO DEN " << " " << kcol << coef_CC[pos] << endl;
                                 //}
                                 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                                 // Mark it as interpolatory
                                 WI[kcol] = -WI[kcol];
                              }
                           }
                        }
                     }
                  }
               } else {
                  // Weak connection increase the denominator with the matrix entry
                  // if it is not yet in the interpolatory set
                  if (WI[jcol] == 0){
                     denom += coef_A[j];
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     //if (inod == DBNODE){
                     //    cout << "+++++++AGGIUNGO XXX " << jcol << " " << coef_A[j] << endl;
                     //}
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     if (fcnodes[jcol] >= 0){
                        // Keep track of weak COARSE neighbours in case they are promoted
                        // to interpolatory points
                        ja_CC[n_CC] = jcol;
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        //if (inod == DBNODE){
                        //   cout << "COARSE WEAK " << jcol <<  endl;
                        //}
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        coef_CC[n_CC] = coef_A[j];
                        WI[jcol] = -(n_CC+1);
                        n_CC++;
                     }
                  } else {
                     // The node has already been added as second level, just
                     // correct the weight
                     int pos = WI[jcol] - 1;
                     coef_CC[pos] = coef_A[j];
                  }
               }
            }
         } // End of loop on dist_1 neigh
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         //if (inod == DBNODE){
         //   cout << endl << endl;
         //   cout << "DIST 1 - NEIGH:" << endl;
         //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++)
         //      if (WI[ja_A[j]] > 0) cout << ja_A[j] << " ";
         //    cout << endl;
         //   cout << "INTERP SET:" << endl;
         //   for (int j = 0; j < n_CC; j++)
         //       if (WI[ja_CC[j]] > 0) cout << ja_CC[j] << " ";
         //   cout << endl;
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         // Init weights of interpolatory points and sparse reset of WI relative to
         // weak COARSE neighbours
         int n_int = 0;
         for (int i = 0; i < n_CC; i++){
            int jcol = ja_CC[i];
            if (WI[jcol] > 0){
               // This is an interpolatory point
               ja_P[ind_P+n_int] = jcol;
               coef_P[ind_P+n_int] = coef_CC[i];
               WI[jcol] = n_int+1;
               n_int++;
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               //if (inod == DBNODE){
               //   cout << "XXXXXXXXXAGGIUNTO CC: " << jcol << " " << WI[jcol] << " " << coef_CC[i] << endl;
               //}
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            } else {
               // This is a weak COARSE neighbour --> reset it
               WI[jcol] = 0;
            }
         }
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         //if (inod == DBNODE){
         //   cout << endl << endl;
         //   cout << "----------------------------------" << endl << endl;
         //   double NUM = 0.0;
         //   double DEN = XXX;
         //   for (int i = 0; i < n_int; i++) NUM+= coef_P[ind_P+i];
         //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){
         //       if (coef_S[j] > 0 && fcnodes[ja_A[j]] < 0) DEN += coef_A[j];
         //   }
         //   cout << "SOMMA NUM: " << NUM << endl;
         //   cout << "SOMMA DEN: " << DEN + DDD << endl;
         //   cout << "----------------------------------" << endl << endl;
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
         // Correct the interpolation weight to take into account strong FINE neighs
         // Loop over strong FINE neighbors
         for (int k = 0; k < n_FS; k++){
            int knod = ja_FS[k];
            double a_ik = coef_FS[k];
            // Init the extended sum with the entry corresponding to inod
            double ext_sum = std::min(ZERO,a_ik);
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            //if (inod == DBNODE){
            //   cout << "FNEIGH " << knod << " INIT ESUM " << ext_sum << endl;
            //}
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            // Explore closest neighbors of knod considering only nodes in the
            // interpolatory set
            int ind = 0;
            for (int l = iat_A[knod]; l < iat_A[knod+1]; l++){
               int lcol = ja_A[l];
               if ( WI[lcol] > 0){
                  // This is a node in the interpolatory set
                  double a_kl = std::min(ZERO,coef_A[l]);
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  //if (inod == DBNODE){
                  //   cout << "ESUM ADD " << lcol << " VAL " << a_kl << endl;
                  //}
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  ext_sum += a_kl;
                  a_kj[ind] = a_kl;
                  pos_kj[ind] = WI[lcol] - 1;
                  ind++;
               }
            }
            // Check that extended sum is not null
            if ( abs(ext_sum) > EPS*a_ii ){
               // Update denominator
               denom += a_ik*std::min(ZERO,a_ik) / ext_sum;
               // Update all the weights
               for (int jj = 0; jj < ind; jj++){
                  int jpos = pos_kj[jj];
                  coef_P[ind_P+jpos] += a_ik*a_kj[jj] / ext_sum;
               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               //if (inod == DBNODE){
                  //cout << "EXT_SUM FINAL" << ext_sum << endl;
                  //cout << "WGT CORR:  " << endl;
                  //for (int jj = 0; jj < ind; jj++){
                  //   int pos = pos_kj[jj];
                  //   cout << pos << " " << ja_P[ind_P+pos] << " " << a_ik*a_kj[jj] << endl;
                  //}
               //}
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            } else {
               // Treat the case of zero extended sum
               double fac = ONE / (static_cast<double>(n_int+1));
               // Update denominator
               denom += a_ik*fac;
               // Update all the weights
               for (int jj = 0; jj < ind; jj++){
                  int jpos = pos_kj[jj];
                  coef_P[ind_P+jpos] += a_ik*fac;
               }
            }
         }
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         //if (inod == DBNODE){
         //   cout << inod << " DENOM " << denom << endl;
         //   cout << "COEF DENOM " << endl;
         //   for (int j = iat_A[inod]; j < iat_A[inod+1]; j++){
         //      if (coef_S[j] < 0){
         //         int jnod = ja_A[j];
         //         cout << jnod << " " << coef_A[j] << " " << (WI[jnod]>0 ? "INT":"NOI") << endl;
         //      }
         //   }
         //}
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         // Adjust column index in ja_P and sort
         for (int i = ind_P; i < ind_P+n_int; i++){
            // Reset of WI
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
            //for (int i = ind_P; i < ind_P+n_int; i++) sum += coef_P[i];
            //cout << "SOMMA  " << sum << endl;
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
         //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
         ind_P += n_int;

      }

      // Update pointer to next row
      //iat_P[inod-shift] = ind_P;
      iat_P[nn_P] = ind_P;

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

   return 0;
}
