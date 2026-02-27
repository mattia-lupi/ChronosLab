#include "lapacke.h"

#include "precision.h"
#include "DebEnv.h"
#include "maxVol.h"
#include "heapsort.h"

//----------------------------------------------------------------------------------------

// Reduce the interpolatory set to allow for a smaller norm of the prolongation row
void Reduce_IntSet(const rExt maxrownrm, const iReg itmax_vol, const rExt tol_vol,
                   const rExt maxcond, const iReg optimal_lwork, const iReg inod,
                   const iReg n_neigh, const iReg ntvecs, iReg &row_rank,
                   const iReg *const fcnodes, iReg *int_list, const iReg *neigh,
                   const rExt *const *const TV, rExt **TVcomp, rExt *WR, rExt *coef_P,
                   rExt &row_nrm){

   while (row_nrm > maxrownrm) {
      if (row_rank == 1){
         row_rank--;
         break;
      }

      // Reload interpolatory points into TVcomp
      iReg n_int = 0;
      for (iReg i = 1; i < n_neigh; i++){
         iReg i_neigh = neigh[i];
         if (fcnodes[i_neigh] >= 0){
            int_list[n_int] = i_neigh;
            for (iReg j = 0; j < ntvecs; j++)
               TVcomp[n_int][j] = TV[i_neigh][j];
            n_int++;
         }
      }

      // Reduce the number of interpolatory points by 1
      int cmax = row_rank-1;

      // Select the best basis using maxVol
      maxVol(cmax,maxcond,itmax_vol,tol_vol,n_int,ntvecs,TVcomp,
             row_rank,int_list);
      // Sort the list
      heapsort(int_list,row_rank);

      // Load the basis in TVcomp
      for (int i = 0; i < row_rank; i++){
         int i_neigh = int_list[i];
         for (int j = 0; j < ntvecs; j++)
            TVcomp[i][j] = TV[i_neigh][j];
      }

      // Load the vector of inod into the coef_P
      for (int j = 0; j < ntvecs; j++)
         coef_P[j] = TV[inod][j];

      // Compute weights
      lapack_int info = LAPACKE_dgels_work(LAPACK_COL_MAJOR,'N',ntvecs,row_rank,1,&(TVcomp[0][0]),
                                ntvecs,coef_P,ntvecs,WR,optimal_lwork);
      if(info != 0){throw linsol_error ("ProlStripe_BAMG","error in LAPACKE_dgels");}

      // Compute row norm
      row_nrm = inl_dnrm2(row_rank,coef_P,1);

      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      /*
      if (DEBUG && BAMG_DEBUG){
         type_OMP_iReg myid = omp_get_thread_num();
         res_row = inl_dnrm2(ntvecs-row_rank,&(coef_P[row_rank]),1);
         fprintf(DebEnv.t_logfile[myid],"REDUCING NORM: ");
         fprintf(DebEnv.t_logfile[myid],"nn %5d row_nrm %15.6e  resid %15.6e\n",cmax,
                 row_nrm,res_row);
         fflush(DebEnv.t_logfile[myid]);
      }
      */
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   }
}
