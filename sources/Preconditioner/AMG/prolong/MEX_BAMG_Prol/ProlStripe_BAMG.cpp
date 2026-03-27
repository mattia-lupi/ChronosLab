#include <vector>
#include <math.h>    // tu use: abs in clapack
#include <algorithm> // tu use: max
#include <lapacke.h>

#include "DebEnv.h"
#include "precision.h"    // to use: iReg, rExt
#include "BAMG_params.h"  // to use: BAMG_params
#include "linsol_error.h" // to thorw errors

#include "check_neigh_cond.h"
#include "add_new_neighs.h"
#include "inl_blas1.h"
#include "maxVol.h"
#include "heapsort.h"
#include "Reduce_IntSet.h"

//----------------------------------------------------------------------------------------

const rExt ONE = 1.0;

// Inner part of BAMG prolongation that is used to distributed work among threads.
void ProlStripe_BAMG(const BAMG_params& params, iReg firstrow_0, iReg firstrow, iReg lastrow,
                     iReg nn_S, iReg ntvecs, const iExt *const iat_S, const iReg *const ja_S,
                     const iReg *const fcnodes,
                     const rExt *const *const TV, iExt &nt_P, iExt *iat_P, iReg *ja_P,
                     rExt *coef_P, iReg *c_mark, iReg *dist_count){

   // Get verbosity flag
   bool VERB_FLAG = params.verbosity >= VLEV_MEDIUM;

   // Extract parameters from params
   iReg itmax_vol = params.itmax_vol;
   iReg dist_min  = params.dist_min;
   iReg dist_max  = params.dist_max;
   iReg mmax      = params.mmax;
   rExt maxcond   = params.maxcond;
   rExt maxrownrm = params.maxrownrm;
   rExt tol_vol   = params.tol_vol;
   rExt eps       = params.eps;
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //int myidN = omp_get_thread_num();
   //char label[6];  
   //sprintf(label, "%05d", myidN);
   //string slabel = label;
   //string name = "OUT_" + slabel;
   //FILE *of = fopen(name.c_str(),"w");
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Query workspace for dgels
   rExt *dummy_double = nullptr;
   iReg optimal_lwork;
   rExt db_lwork;
   lapack_int info = LAPACKE_dgels_work(LAPACK_COL_MAJOR,'N',ntvecs,ntvecs,1,
                                        dummy_double,ntvecs,dummy_double,ntvecs,
                                        &db_lwork,-1);
   if (info != 0)
      throw linsol_error ("ProlStripe_BAMG","Quering LAPACK DGELS workspace");
   optimal_lwork = static_cast<int>(db_lwork);

   // Allocate some scratch vectors
   std::vector<iReg> vec_int_list;
   std::vector<iReg> vec_neigh;
   std::vector<iReg> vec_WI;
   std::vector<rExt> vec_WR;
   try {
      vec_int_list.resize(nn_S);
      vec_neigh.resize(nn_S);
      vec_WI.resize(nn_S);
      vec_WR.resize(optimal_lwork);
   } catch (linsol_error) {
      throw linsol_error ("ProlStripe_BAMG","allocating scratches");
   }
   iReg *int_list = vec_int_list.data();
   iReg *neigh = vec_neigh.data();
   iReg *WI = vec_WI.data();
   rExt *WR = vec_WR.data();

   // Allocate temporary space for "compressed" test space (restricted to neighbours only)
   std::vector<rExt*> vec_TVcomp;
   std::vector<rExt>  vec_TVbuf;
   try {
      vec_TVcomp.resize(nn_S+1);
      vec_TVbuf.resize(ntvecs*(nn_S+1));
   } catch (linsol_error) {
      throw linsol_error ("ProlStripe_BAMG","allocating compressed TV");
   }
   rExt **TVcomp = vec_TVcomp.data();
   rExt  *TVbuf  = vec_TVbuf.data();

   iReg kk = 0;
   for (iReg i = 0; i < nn_S+1; i++){
      TVcomp[i] = &(TVbuf[kk]);
      kk += ntvecs;
   }

   // Init WI
   for (iReg i = 0; i < nn_S; i++) WI[i] = 0;

   // Init pointer to the prolongation matrix
   iExt ind_P = 0;
   iat_P[0] = ind_P;

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (DEBUG && BAMG_DEBUG){
      type_OMP_iReg myid = omp_get_thread_num();
      fprintf(DebEnv.t_logfile[myid],"---------------- BEGIN OF PROLONGATION -----------------\n");
      fprintf(DebEnv.t_logfile[myid],"OFFSET of current partition: %d\n",firstrow_0);
      fprintf(DebEnv.t_logfile[myid]," %9s %6s %6s %6s %6s %9s %9s\n","INOD","C/F",
              "dist","n_int","rank","row_norm","res_nrm");
      fflush(DebEnv.t_logfile[myid]);
   }
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Loop over all the nodes
   iReg shift = firstrow - 1;
   for (iReg inod = firstrow; inod < lastrow ; inod ++){ // Node_loop

      iReg inod_coarse = fcnodes[inod];

      // Check whether inod is FINE OR COARSE
      if (inod_coarse >= 0){

         // It is a coarse node
         ja_P[ind_P] = inod_coarse;
         coef_P[ind_P] = ONE;
         ind_P++;

         //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if (DEBUG && BAMG_DEBUG){
            type_OMP_iReg myid = omp_get_thread_num();
            fprintf(DebEnv.t_logfile[myid]," %9d COARSE %6s %6s %6s %9s %9s %8s\n",
                    inod-firstrow_0,"0","0","0","0","0","-");
            fflush(DebEnv.t_logfile[myid]);
         }
         //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      } else {

         // It is a fine node
         iReg row_rank;

         // Add distance zero neighbours (inod itself)
         iReg distance = 0;
         iReg n_neigh = 1;
         neigh[0] = inod;
         WI[inod] = 1;

         // Add the list of neighbour up to distance dist_min
         iReg istart_neigh = 0;
         iReg iend_neigh = 1;
         while (distance < dist_min){
            distance++;
            try {
               add_new_neighs(iat_S,ja_S,istart_neigh,iend_neigh,nn_S,n_neigh,neigh,WI);
            } catch (linsol_error) {
               throw linsol_error ("ProlStripe_BAMG",
                                   "add the list of neighbour up to distance dist_min");
            }
            if (n_neigh == iend_neigh){
               // No new neighs have been added ==> break
               break;
            }
            istart_neigh = iend_neigh;
            iend_neigh = n_neigh;
         }

         // Compute norm of current row
         rExt bnorm = inl_dnrm2(ntvecs,TV[inod],1);
         rExt res_ass = eps*bnorm;
         rExt res_row;
         rExt row_nrm;

         bool UPD_WEIGHTS = true;
         while ( UPD_WEIGHTS ){ // UPD_WEIGHTS Loop

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (DEBUG && BAMG_DEBUG){
               type_OMP_iReg myid = omp_get_thread_num();
               fprintf(DebEnv.t_logfile[myid]," %9d   FINE",inod-firstrow_0);
               fflush(DebEnv.t_logfile[myid]);
            }
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            // Select possible interpolatory points from the list of neighbours
            iReg n_int = 0;
            for (iReg i = 1; i < n_neigh; i++){
               iReg i_neigh = neigh[i];
               if (fcnodes[i_neigh] >= 0){
                  int_list[n_int] = i_neigh;
                  for (iReg j = 0; j < ntvecs; j++) TVcomp[n_int][j] = TV[i_neigh][j];
                  n_int++;
               }
            }
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            //fprintf(of,"%5d %5d\n",n_neigh,n_int);
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            if (n_int > 0){

               // Select the best basis using maxVol
               maxVol(mmax,maxcond,itmax_vol,tol_vol,n_int,ntvecs,TVcomp,row_rank,int_list);

               // Sort the list
               heapsort(int_list,row_rank);

               // Load the basis in TVcomp
               for (iReg i = 0; i < row_rank; i++){
                  iReg i_neigh = int_list[i];
                  for (iReg j = 0; j < ntvecs; j++) TVcomp[i][j] = TV[i_neigh][j];
               }

               // Load the vector of inod into the coef_P
               for (iReg j = 0; j < ntvecs; j++) coef_P[ind_P+j] = TV[inod][j];

               // Compute weights
               info = LAPACKE_dgels_work(LAPACK_COL_MAJOR,'N',ntvecs,row_rank,1,
                                         &(TVcomp[0][0]),ntvecs,&(coef_P[ind_P]),ntvecs,
                                         WR,optimal_lwork);
               if(info != 0){throw linsol_error ("ProlStripe_BAMG","error in LAPACKE_dgels");}

               // Compute residual
               res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);

               // Compute row norm
               row_nrm = inl_dnrm2(row_rank,&(coef_P[ind_P]),1);

            } else {

               // There is no interpolatory set
               row_rank = 0;
               res_row = bnorm;
               row_nrm = 0.0;

            }

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (DEBUG && BAMG_DEBUG){
               type_OMP_iReg myid = omp_get_thread_num();
               if (n_int > 0){
                  fprintf(DebEnv.t_logfile[myid]," %6d %6d %6d %9.2e %9.2e",distance,n_int,
                          row_rank,row_nrm,res_row/bnorm);
               } else {
                  fprintf(DebEnv.t_logfile[myid]," %6d %6d %6d %9.2e %9.2e",distance,n_int,
                          row_rank,row_nrm,1.0);
               }
               iReg inod0 = inod-firstrow_0;
               if (inod0 == 6){
                  fprintf(DebEnv.r_logfile,"row_nrm %e | res_row %e\n",row_nrm,res_row);
                  fprintf(DebEnv.r_logfile,"NEIGHS: %d\n",n_int);
                  for (int i = 0; i < n_int; i++){
                     fprintf(DebEnv.r_logfile," %d %d\n",int_list[i],fcnodes[int_list[i]]);
                  }
                  fprintf(DebEnv.r_logfile,"---------------\n");
                  for (int i = 0; i < row_rank; i++){
                     fprintf(DebEnv.r_logfile," %d",fcnodes[int_list[i]]);
                  }
                  fprintf(DebEnv.r_logfile,"\n");
                  for (int i = 0; i < row_rank; i++){
                     fprintf(DebEnv.r_logfile," %e",coef_P[ind_P+i]);
                  }
                  fprintf(DebEnv.r_logfile,"\n");
               }
            }
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            // Check residual
            if ( row_nrm < maxrownrm && ((res_row < res_ass) || (row_rank == ntvecs)) ){

               // The residual is low or the number of interpolatory points equals
               // the basis dimension ==> this node has been successfully interpolated
               UPD_WEIGHTS = false;
               //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if (DEBUG && BAMG_DEBUG){
                  type_OMP_iReg myid = omp_get_thread_num();
                  fprintf(DebEnv.t_logfile[myid]," %10s\n","OK");
                  fflush(DebEnv.t_logfile[myid]);
               }
               if (VERB_FLAG){
                  // Increase count of nodes interpolated at a given distance
                  #pragma omp atomic update
                  dist_count[distance-1]++;
               }
               //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            } else {

               // The residual is still high ==> try to increase distance if possible
               if (distance < dist_max){
                  distance++;
                  try {
                     add_new_neighs(iat_S,ja_S,istart_neigh,iend_neigh,nn_S,n_neigh,neigh,WI);
                  } catch (linsol_error) {
                     throw linsol_error ("ProlStripe_BAMG","try to increase distance if possible");
                  }
                  if (n_neigh == iend_neigh){
                     // No new neighs have been added ==> break
                     UPD_WEIGHTS = false;
                     if (n_int == 0){
                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        if (DEBUG && BAMG_DEBUG){
                           type_OMP_iReg myid = omp_get_thread_num();
                           fprintf(DebEnv.t_logfile[myid]," %10s\n","NO NEIGH");
                           fflush(DebEnv.t_logfile[myid]);
                        }
                        if (VERB_FLAG){
                           // Increase the count of fine nodes without neighbours
                           #pragma omp atomic update
                           dist_count[dist_max+2]++;
                        }
                        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     } else {
                        // The interpolatory set is not empty, try to interpolate with a
                        // smaller row norm, if necessary
                        if (row_nrm > RELAX_FAC*maxrownrm){
                           Reduce_IntSet(RELAX_FAC*maxrownrm,itmax_vol,tol_vol,maxcond,
                                         optimal_lwork,inod,n_neigh,ntvecs,row_rank,
                                         fcnodes,int_list,neigh,TV,TVcomp,WR,
                                         &(coef_P[ind_P]),row_nrm);
                        }
                        res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);
                        if (res_row < res_ass){
                           if (VERB_FLAG){
                              // Increase the count of fine nodes with large norm
                              #pragma omp atomic update
                              dist_count[dist_max+1]++;
                           }
                        } else {
                           // Mark the node as one to be promoted
                           c_mark[inod] = 1;
                           if (VERB_FLAG){
                              // Increase the count of fine nodes with large error
                              #pragma omp atomic update
                              dist_count[dist_max]++;
                           }
                        }
                     }

                  } else {
                     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     if (DEBUG && BAMG_DEBUG){
                        type_OMP_iReg myid = omp_get_thread_num();
                        fprintf(DebEnv.t_logfile[myid]," %10s\n","   CYCLE");
                        fflush(DebEnv.t_logfile[myid]);
                     }
                     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  }
                  istart_neigh = iend_neigh;
                  iend_neigh = n_neigh;

               } else {

                  // Maximal distance have been already reached ==> break
                  UPD_WEIGHTS = false;

                  // Try to interpolate with a smaller row norm if the norm is large
                  if (row_nrm > RELAX_FAC*maxrownrm)
                     Reduce_IntSet(RELAX_FAC*maxrownrm,itmax_vol,tol_vol,maxcond,
                                   optimal_lwork,inod,n_neigh,ntvecs,row_rank,fcnodes,
                                   int_list,neigh,TV,TVcomp,WR,&(coef_P[ind_P]),row_nrm);
                  // Compute residual
                  res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);
                  // Mark the node as one to be promoted to coarse
                  if (res_row > res_ass) c_mark[inod] = 1;

                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
                  if (DEBUG && BAMG_DEBUG){
                     iReg inod0 = inod-firstrow_0;
                     if (inod0 == 6){
                        type_OMP_iReg myid = omp_get_thread_num();
                        rExt res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[ind_P+row_rank]),1);
                        fprintf(DebEnv.t_logfile[myid],"\nres_row %e row_nrm %e row_rank %d\n",
                                res_row,row_nrm,row_rank);
                     }
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2

                  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if (DEBUG && BAMG_DEBUG){
                     type_OMP_iReg myid = omp_get_thread_num();
                     if (res_row > res_ass){
                        fprintf(DebEnv.t_logfile[myid]," %10s\n","REDUC. SET");
                        fprintf(DebEnv.t_logfile[myid]," %9d   FINE",inod-firstrow_0);
                        fprintf(DebEnv.t_logfile[myid]," %6d %6d %6d %9.2e %9.2e",distance,
                                n_int,row_rank,row_nrm,res_row/bnorm);
                        fprintf(DebEnv.t_logfile[myid]," %10s\n","LARGE ERR");
                        fprintf(DebEnv.t_logfile[myid],"NEIGHLIST: ");
                        for (iReg i = 0; i < n_neigh; i++){
                          if (fcnodes[neigh[i]] >= 0)
                             fprintf(DebEnv.t_logfile[myid]," %8d",neigh[i]);
                        }
                        fprintf(DebEnv.t_logfile[myid],"\n");
                     } else {
                        fprintf(DebEnv.t_logfile[myid]," %10s\n","REDUC. SET");
                        fprintf(DebEnv.t_logfile[myid]," %9d   FINE",inod-firstrow_0);
                        fprintf(DebEnv.t_logfile[myid]," %6d %6d %6d %9.2e %9.2e",distance,
                                n_int,row_rank,row_nrm,res_row/bnorm);
                        fprintf(DebEnv.t_logfile[myid]," %10s\n","LARGE NRM");
                        fprintf(DebEnv.t_logfile[myid],"NEIGHLIST: ");
                        for (iReg i = 0; i < n_neigh; i++){
                          if (fcnodes[neigh[i]] >= 0)
                             fprintf(DebEnv.t_logfile[myid]," %8d",neigh[i]);
                        }
                        fprintf(DebEnv.t_logfile[myid],"\n");

                     }
                     fflush(DebEnv.t_logfile[myid]);
                  }
                  if (VERB_FLAG){
                     if (res_row > res_ass){
                        // Increase count of nodes with high error
                        #pragma omp atomic update
                        dist_count[dist_max]++;
                     } else {
                        // Increase count of nodes with large norm
                        #pragma omp atomic update
                        dist_count[dist_max+1]++;
                     }
                  }
                  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               }

            }

         } // end UPD_WEIGHTS Loop
         if (VERB_FLAG){
            if (row_rank < ntvecs){
               #pragma omp atomic update
               dist_count[dist_max+3]++;
            }
         }

         // Weights are already in coef_S, store column indices in P
         for (iReg i = 0; i < row_rank; i++){
            ja_P[ind_P] = fcnodes[int_list[i]];
            ind_P++;
         }

         // Sparse reset of WI
         for (iReg i = 0; i < n_neigh; i++) WI[neigh[i]] = 0;

      }

      // Update pointer to next row
      iat_P[inod-shift] = ind_P;

   } // Node_loop

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (DEBUG && BAMG_DEBUG){
      type_OMP_iReg myid = omp_get_thread_num();
      fprintf(DebEnv.t_logfile[myid],"---------------- END OF PROLONGATION -------------------\n");
      fflush(DebEnv.t_logfile[myid]);
   }
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // Count number of non-zeroes
   nt_P = ind_P;
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //fclose(of);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@

}
