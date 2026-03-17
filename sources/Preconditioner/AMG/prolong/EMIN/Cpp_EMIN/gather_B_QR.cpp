#include <stdlib.h>
#include <omp.h>
#include <limits>
#include "emin_blas.h"     // to use: DGEMV
//////////////////////////////////
#include "inl_blas1.h"
#include <iostream>
#include <stdio.h>
//////////////////////////////////

#include "DebEnv.h"
#include "parm_EMIN.h"

#define PRINT_LOC_INFO 0
#define PRINT_LOC_INFO_ALL 0

/*****************************************************************************************
 *
 * This function gathers from the prolongation pattern (given by rows) and the test vector
 * array TV, the block diagonal matrix B that is used to enforce the constraint. B is
 * immediately factorized with QR and only Q is returned.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> allocation error for private scratches
 * 3 ---> allocation error for the final output
 * 4 ---> lapack error
 *
*****************************************************************************************/
int gather_B_QR(const int np, const double condmax, const int nn, const int nn_C,
                const int ntv, const int *fcnode, const int *iat_patt, const int *ja_patt,
                const double *const *TV, double *&mat_Q, double *coef_P0)
{
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// g E tau POSSONO ESSERE DEI VETTORI LOCALI DI DIMENSIONE NTV CHE POI VENGONO CANCELLATI
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //FILE *bbf = fopen("LOG_gather","w");
   //FILE *of = fopen("P0_prima","w");
   //for (int i = 0; i < iat_patt[nn]; i++) fprintf(of,"%20.11e\n",coef_P0[i]);
   //fclose(of);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Init error code
   int ierr = 0;

   // Allocate Q
   int nrows_Q = iat_patt[nn];
   int nterm_Q = ntv*nrows_Q;
   mat_Q = (double*) malloc( nterm_Q*sizeof(double) );
   if (mat_Q == nullptr) return ierr = 3;

   /* form of op(A) & op(B) to use in matrix vector multiplication */
   char const *chn = "N", *cht = "T";
   /* scalar values to use in dgemv */
   double const one = 1.0, mone = -1.0, zero = 0.0;
   lapack_int const oneint = 1;

   // Allocate shared scratches
   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   if (c2glo == nullptr)  return ierr = 1;

   #pragma omp parallel num_threads(np)
   {
      // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstcol, ncolth, lastcol;
      if (mythid <= resto) {
         ncolth = bsize+1;
         firstcol = mythid*ncolth;
         if (mythid == resto) ncolth--;
      } else {
         ncolth = bsize;
         firstcol = mythid*bsize + resto;
      }
      lastcol = firstcol + ncolth;

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double *BB_blk;
      double *QQ_blk;
      double *v1;
      double *v2;
      double *v3;
      FILE *log;
      FILE *reslog;
      if (DEBUG){
         int NRMAX = 0;
         for (int i = 0; i < nn; i++) NRMAX = std::max(NRMAX,iat_patt[i+1]-iat_patt[i]);
         BB_blk = (double*) malloc(NRMAX*6*sizeof(double));
         QQ_blk = (double*) malloc(NRMAX*6*sizeof(double));
         v1 = (double*) malloc(NRMAX*sizeof(double));
         v2 = (double*) malloc(NRMAX*sizeof(double));
         v3 = (double*) malloc(NRMAX*sizeof(double));
         char filename[100];
         sprintf(filename, "LOGGO_%02d",mythid);
         log = fopen(filename,"w");
         sprintf(filename, "RANDRES_%02d",mythid);
         reslog = fopen(filename,"w");
      }
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      // Estimate number of rows and first position for this chunk of columns
      int pos_g = iat_patt[firstcol];
      int pos_Q = pos_g*ntv;

      // Set proper position in Q and g
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if (DEBUG){
         fprintf(DebEnv.t_logfile[mythid],"XX-->pos_g %d size nrows_Q %d true %d\n",pos_g,nrows_Q,nn*ntv);
      }
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      double *g_scr = (double*) malloc( ntv*sizeof(double) );
      double *BB_scr = &(mat_Q[pos_Q]);

      // Find the largest number of rows in a block
      int nrmax_blk = 0;
      int istart, iend;
      iend = iat_patt[firstcol];
      for (int i = firstcol+1; i <= lastcol; i++){
         istart = iend;
         iend = iat_patt[i];
         nrmax_blk = std::max(nrmax_blk,iend-istart);
      }
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if (DEBUG){
         fprintf(DebEnv.t_logfile[mythid],"NEL CALCOLO LWORK\n");
         fprintf(DebEnv.t_logfile[mythid],"lnrmax_blk %d\n", nrmax_blk);
         fflush(DebEnv.t_logfile[mythid]);
      }
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      // Query work space for DGEQRF and DORGQR
      lapack_int ierr_lapack;
      lapack_int l_nn = static_cast<lapack_int>(std::max(ntv,nrmax_blk));
      lapack_int l_mm = static_cast<lapack_int>(ntv);
      lapack_int l_kk = static_cast<lapack_int>(ntv);
      double query_work;
      double *SIGMA = nullptr;      // (DI DIMENSIONE l_nn == nr_BB_loc)
      double *dummy_U = nullptr;    // (FASULLO NON VIENE USATO)
      double *VT = nullptr;         // (DI DIMENSIONE l_mm*l_mm = ntv*ntv)  
      double *tau = nullptr;
      double *work = nullptr;

      lapack_int lwork = -1;
      char const *cho = "o", *chs = "s";
      dgesvd(cho,chs,&l_nn,&l_mm,BB_scr,&l_nn,SIGMA,dummy_U,&l_nn,VT,&l_mm,
             &query_work,&lwork,&ierr_lapack);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      if (ierr > 0) goto exit_pragma;
      lwork = static_cast<lapack_int>(std::max(query_work,static_cast<double>(ntv)));

      // Allocate private workspace
      SIGMA = (double*) malloc( std::max(ntv,nrmax_blk)*sizeof(double) );
      VT = (double*) malloc( ntv*ntv*sizeof(double) );
      tau = (double*) malloc( (ncolth*ntv)*sizeof(double) );
      work = (double*) malloc( lwork*sizeof(double) );
      if (SIGMA == nullptr || VT == nullptr || tau == nullptr || work == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      if (ierr > 0) goto exit_pragma;

      // Create mapping from coarse node numbering to global (original) numbering
      #pragma omp for
      for (int i = 0; i < nn; i++){
         int k = fcnode[i];
         if (k >= 0) c2glo[k] = i;
      }

      // Loop over the current chunk of columns
      int ind_BB, ind_tau, nnz_BB;
      ind_BB = 0;
      ind_tau = 0;
      int istart_patt, iend_patt;
      for (int icol = firstcol; icol < lastcol; icol++){
         // Check that this is a FINE node
         if (fcnode[icol] < 0){

            istart_patt = iat_patt[icol];
            iend_patt = iat_patt[icol+1];
            int nr_BB_loc = iend_patt-istart_patt;

            // Check that the row is not empty
            if (nr_BB_loc > 0){
               // Copy TV into BB_scr
               int k = ind_BB;
               for (int i = istart_patt; i < iend_patt; i++){
                  int i_F = c2glo[ja_patt[i]];
                  int kk = k;
                  for (int j = 0; j < ntv; j++){
                     BB_scr[kk] = TV[i_F][j];
                     kk += nr_BB_loc;
                  }
                  k++;
               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               if (DEBUG){
                  for (int k = 0; k < nr_BB_loc*ntv; k++) BB_blk[k] = BB_scr[ind_BB+k];
               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               // Copy TV into g
               for (int j = 0; j < ntv; j++) g_scr[j] = TV[icol][j];

               // Compute g = g - BB^T*coef_P0
               lapack_int b_m = static_cast<lapack_int>( nr_BB_loc );
               lapack_int b_n = static_cast<lapack_int>( ntv );
               dgemv(cht, &b_m, &b_n, &mone, &(BB_scr[ind_BB]), &b_m,
                     &(coef_P0[istart_patt]), &oneint, &one, g_scr, &oneint);
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               if (DEBUG){
                  fprintf(DebEnv.t_logfile[mythid],"-----------------------------------\n");
                  fprintf(DebEnv.t_logfile[mythid],"GLOBAL ROW OF A (ICOL): %d %d %d\n",icol,nn,ncolth);
                  fprintf(DebEnv.t_logfile[mythid],"Size of B_loc: %d %d\n",nr_BB_loc,ntv);
                  fprintf(DebEnv.t_logfile[mythid],"\nB_loc:\n");
                  for (int i = 0; i < nr_BB_loc; i++){
                     for (int j = 0; j < ntv; j++)
                        fprintf(DebEnv.t_logfile[mythid]," %17.10e",BB_scr[ind_BB+j*nr_BB_loc+i]);
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                  }
                  fprintf(DebEnv.t_logfile[mythid],"\nrhs init:\n");
                  for (int j = 0; j < ntv; j++)
                     fprintf(DebEnv.t_logfile[mythid]," %15.6e",TV[icol][j]);
                  fprintf(DebEnv.t_logfile[mythid],"\nP0 init:\n");
                  for (int j = 0; j < nr_BB_loc; j++)
                     fprintf(DebEnv.t_logfile[mythid]," %15.6e",coef_P0[istart_patt+j]);
                  fprintf(DebEnv.t_logfile[mythid],"\nrhs:\n");
                  for (int j = 0; j < ntv; j++)
                     fprintf(DebEnv.t_logfile[mythid]," %15.6e",g_scr[j]);
                  fprintf(DebEnv.t_logfile[mythid],"\n");
                  fflush(DebEnv.t_logfile[mythid]);
               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

               //----------------------------------+
               // Try the solution of full-rank BB |
               //----------------------------------+

               bool FAIL_QR = false;
               lapack_int i_lpk_1, i_lpk_2, i_lpk_3;
               if (nr_BB_loc >= ntv){

                  // Perform QR on BB
                  l_nn = static_cast<lapack_int>(nr_BB_loc);
                  lapack_int l_ll = l_nn;
                  dgeqrf(&l_nn,&l_mm,&(BB_scr[ind_BB]),&l_ll,&(tau[ind_tau]),work,&lwork,&i_lpk_1);

                  // Check conditioning of the resulting R
                  double max_DR = 0.0;
                  double min_DR = std::numeric_limits<double>::max();
                  for (int kk = 0; kk < l_mm; kk++){
                     max_DR = std::max(max_DR,fabs(BB_scr[ind_BB+kk*l_ll+kk]));
                     min_DR = std::min(min_DR,fabs(BB_scr[ind_BB+kk*l_ll+kk]));
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
                  if (DEBUG){
                     fprintf(DebEnv.t_logfile[mythid],"DIAG R: ");
                     for (int kk = 0; kk < l_mm; kk++)
                        fprintf(DebEnv.t_logfile[mythid]," %15.6e",fabs(BB_scr[ind_BB+kk*l_ll+kk]));
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                     fprintf(DebEnv.t_logfile[mythid],
                             "%6d max_DR %15.6e min_DR %15.6e COND_1 %15.6e\n",icol,
                             max_DR,min_DR,max_DR/min_DR);
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
                  if (max_DR / min_DR > condmax){

                     //------------------------------------+ 
                     // BB is rank deficient use SVD below |
                     //------------------------------------+ 

                     FAIL_QR = true;
                     // Reload TV in BB_scr
                     int k = ind_BB;
                     for (int i = istart_patt; i < iend_patt; i++){
                        int i_F = c2glo[ja_patt[i]];
                        int kk = k;
                        for (int j = 0; j < ntv; j++){
                           BB_scr[kk] = TV[i_F][j];
                           kk += nr_BB_loc;
                        }
                        k++;
                     }
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     if ( PRINT_LOC_INFO ) std::cout << icol <<
                        " conditioning larger than threshold: "
                        << max_DR / min_DR << " > " << condmax << std::endl;
                     if (DEBUG) fprintf(DebEnv.t_logfile[mythid],
                            "%6d COND_LRG: max %15.6e min %15.6e cond %15.6e\n",
                            icol,max_DR,min_DR,max_DR / min_DR);
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  } else {

                     //-----------------+
                     // BB is full rank |
                     //-----------------+

                     // Solve transposed triangular system g = inv(RR')*g
                     char const *chu = "U";
                     dtrtrs(chu,cht,chn,&l_mm,&oneint,&(BB_scr[ind_BB]),&l_ll,g_scr,&l_mm,&i_lpk_2);

                     // Transform QQ from Householder rotation to standard form
                     dorgqr(&l_nn,&l_mm,&l_kk,&(BB_scr[ind_BB]),&l_ll,&(tau[ind_tau]),work,&lwork,&i_lpk_3);
                     if (i_lpk_1 || i_lpk_2 || i_lpk_3){
                        #pragma omp atomic write
                        ierr = 4;
                        goto exit_loop_icol;
                     }

                     // Update coef_P0 to ensure the TV constraint: coef_P0 += QQ*g
                     //XXXXX METTERE UN CONTROLLO CHE COEF_P0 NON DIVENTI TROPPO GRANDE
                     b_m = static_cast<lapack_int>( nr_BB_loc );
                     b_n = static_cast<lapack_int>( ntv );
                     dgemv(chn, &b_m, &b_n, &one, &(BB_scr[ind_BB]), &b_m,
                           g_scr, &oneint, &one, &(coef_P0[istart_patt]), &oneint);

                     // Compute the number of entries of Q
                     nnz_BB = nr_BB_loc*ntv;
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     if (DEBUG){
                        fprintf(DebEnv.t_logfile[mythid],"\nQ:\n");
                        for (int i = 0; i < nr_BB_loc; i++){
                           for (int j = 0; j < std::min(ntv,nr_BB_loc); j++)
                              fprintf(DebEnv.t_logfile[mythid]," %17.10e",BB_scr[ind_BB+j*nr_BB_loc+i]);
                           fprintf(DebEnv.t_logfile[mythid],"\n");
                        }
                        fprintf(DebEnv.t_logfile[mythid],"\nDP:\n");
                        for (int j = 0; j < nr_BB_loc; j++)
                          fprintf(DebEnv.t_logfile[mythid]," %15.6e",coef_P0[istart_patt+j]);
                        fprintf(DebEnv.t_logfile[mythid],"\n");
                        fflush(DebEnv.t_logfile[mythid]);
                     }
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  }

               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               int TRUERANK = ntv;
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

               //--------------------------------------------+
               // Perform least square solution if necessary |
               //--------------------------------------------+

               if ( (nr_BB_loc < ntv) || FAIL_QR){

                  if ( PRINT_LOC_INFO )
                     std::cout << icol << " Solution for RANK deficient system " << std::endl;

                  // Compute SVD of BB
                  l_nn = static_cast<lapack_int>(nr_BB_loc);
                  lapack_int l_ll = l_nn;
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  if (DEBUG){
                     fprintf(DebEnv.t_logfile[mythid],"PRIMA DI DGESVD\n");
                     fprintf(DebEnv.t_logfile[mythid],"l_nn %ld l_mm %ld l_ll %ld lwork %ld\n",
                             (long int) l_nn,(long int) l_mm,(long int) l_ll,(long int) lwork);
                     fflush(DebEnv.t_logfile[mythid]);
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  dgesvd(cho,chs,&l_nn,&l_mm,&(BB_scr[ind_BB]),&l_ll,SIGMA,dummy_U,&l_ll,VT,
                         &l_mm,work,&lwork,&i_lpk_1);
                  if (i_lpk_1){
                     #pragma omp atomic write
                     ierr = 4;
                     goto exit_loop_icol;
                  }

                  // Compute the rank of BB using SIGMA
                  int rank_BB = 1;
                  if (nr_BB_loc > 1){
                     while (SIGMA[0] < condmax *abs(SIGMA[rank_BB])){
                        rank_BB++;
                        if (rank_BB == std::min(ntv,nr_BB_loc)){
                           break;
                        }
                     }
                  }

                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  TRUERANK = rank_BB;
                  if (DEBUG){
                     if (nr_BB_loc < ntv) fprintf(DebEnv.t_logfile[mythid],"FAT B\n");
                     //fprintf(DebEnv.t_logfile[mythid],"\nRANK B: %d\n",rank_BB);
                     fprintf(DebEnv.t_logfile[mythid],"\n%6d COND_RID %15.6e RANK B %3d\n",icol,
                             SIGMA[0]/SIGMA[rank_BB-1],rank_BB);
                     fprintf(DebEnv.t_logfile[mythid],"\nU = Q:\n");
                     for (int i = 0; i < nr_BB_loc; i++){
                        for (int j = 0; j < std::min(ntv,nr_BB_loc); j++)
                           fprintf(DebEnv.t_logfile[mythid]," %15.6e",BB_scr[ind_BB+j*nr_BB_loc+i]);
                        fprintf(DebEnv.t_logfile[mythid],"\n");
                     }
                     fprintf(DebEnv.t_logfile[mythid],"\nSIGMA:\n");
                     for (int j = 0; j < std::min(ntv,nr_BB_loc); j++)
                        fprintf(DebEnv.t_logfile[mythid]," %15.6e",SIGMA[j]);
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                     fprintf(DebEnv.t_logfile[mythid],"\nVT:\n");
                     for (int i = 0; i < std::min(nr_BB_loc,ntv); i++){
                        for (int j = 0; j < ntv; j++)
                           fprintf(DebEnv.t_logfile[mythid]," %15.6e",VT[j*ntv+i]);
                        fprintf(DebEnv.t_logfile[mythid],"\n");
                     }
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  // Compute w = V^T * g, store w temporarily in work
                  b_m = static_cast<lapack_int>( rank_BB );
                  b_n = static_cast<lapack_int>( ntv );
                  dgemv(chn, &b_m, &b_n, &one, VT, &b_n,
                        g_scr, &oneint, &zero, work, &oneint);

                  // Scale w with the inverse of singular values
                  for (int i = 0; i < rank_BB; i++) work[i] /= SIGMA[i];

                  // Compute coef_P0 += V * w
                  //XXXXX METTERE UN CONTROLLO CHE COEF_P0 NON DIVENTI TROPPO GRANDE
                  b_m = static_cast<lapack_int>( nr_BB_loc );
                  b_n = static_cast<lapack_int>( rank_BB );
                  dgemv(chn, &b_m, &b_n, &one, &(BB_scr[ind_BB]), &b_m,
                        work, &oneint, &one, &(coef_P0[istart_patt]), &oneint);
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  if (DEBUG){
                     fprintf(DebEnv.t_logfile[mythid],"\nDP:\n");
                     for (int j = 0; j < nr_BB_loc; j++)
                       fprintf(DebEnv.t_logfile[mythid]," %15.6e",coef_P0[istart_patt+j]);
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                     fflush(DebEnv.t_logfile[mythid]);
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  if (true) { //@@@@ DA VALUTARE SE METTERLO TRA I PARAMETRI
                     // Pad BB_scr with zeroes
                     int npad = (ntv-rank_BB)*nr_BB_loc;
                     for (int k = 0; k < npad; k++) BB_scr[ind_BB+nr_BB_loc*rank_BB+k] = 0.0;
                  }

                  // Compute the number of entries of Q
                  nnz_BB = nr_BB_loc*ntv;

               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               if (DEBUG){
                  //fprintf(log,"nr_BB_loc %d\n",nr_BB_loc);
                  for (int k = 0; k < nr_BB_loc*ntv; k++) QQ_blk[k] = BB_scr[ind_BB+k];
                  for (int k = 0; k < nr_BB_loc; k++) v1[k] = (double) rand() / (double) RAND_MAX;
                  // Compute v2 = QT*v1
                  for (int i = 0; i < ntv; i++){
                     v2[i] = 0.0;
                     for (int j = 0; j < nr_BB_loc; j++){
                        v2[i] += QQ_blk[i*nr_BB_loc+j]*v1[j];
                     }
                  }
                  // Compute v3 = v1 - Q*v2
                  for (int i = 0; i < nr_BB_loc; i++){
                     v3[i] = v1[i];
                     for (int j = 0; j < ntv; j++){
                        v3[i] -= QQ_blk[j*nr_BB_loc+i]*v2[j];
                     }
                  }
                  // Compute v2 = BT*v3
                  for (int i = 0; i < ntv; i++){
                     v2[i] = 0.0;
                     for (int j = 0; j < nr_BB_loc; j++){
                        v2[i] += BB_blk[i*nr_BB_loc+j]*v3[j];
                     }
                  }
                  // Compute norm of v2
                  double norm_v2 = 0.0;
                  for (int i = 0; i < ntv; i++) norm_v2 += v2[i]*v2[i];
                  norm_v2 = sqrt(norm_v2);
                  rExt rownorm = 0.0;
                  for (int j = 0; j < nr_BB_loc; j++)
                     rownorm += coef_P0[istart_patt+j]*coef_P0[istart_patt+j];
                  rownorm = sqrt(rownorm);
                  fprintf(log,"%8d NORM RANDRES %d %20.10e %20.10e\n",icol,TRUERANK,norm_v2,rownorm);
                  fprintf(reslog,"%8d %d %20.10e %20.10e\n",icol,TRUERANK,norm_v2,rownorm);
                  if (norm_v2 > 1.e-8){
                     fprintf(log,"norm_res %20.10e\n",norm_v2);
                     fprintf(log,"v rand:\n");
                     for (int i = 0; i < nr_BB_loc; i++) fprintf(log,"%20.10e\n",v1[i]);
                     fprintf(log,"BB:\n");
                     for (int i = 0; i < nr_BB_loc; i++){
                        for (int j = 0; j < ntv; j++){
                           fprintf(log," %20.10e",BB_blk[j*nr_BB_loc+i]);
                        }
                        fprintf(log,"\n");
                     }
                     fprintf(log,"QQ:\n");
                     for (int i = 0; i < nr_BB_loc; i++){
                        for (int j = 0; j < ntv; j++){
                           fprintf(log," %20.10e",QQ_blk[j*nr_BB_loc+i]);
                        }
                        fprintf(log,"\n");
                     }
                  }
                  rExt norm = 0.0;
                  for (int j = 0; j < nr_BB_loc; j++)
                     norm += coef_P0[istart_patt+j]*coef_P0[istart_patt+j];
                  fprintf(DebEnv.t_logfile[mythid],"ICOL %d ROWNORM %e\n",icol,sqrt(norm));
               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

               // Update pointers in BB, g and tau
               ind_BB += nnz_BB;
               ind_tau += ntv;

            }

         }

      } // End loop over columns
      exit_loop_icol: ;

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if (DEBUG){
         free(BB_blk);
         free(QQ_blk);
         free(v1);
         free(v2);
         free(v3);
         fclose(log);
         fclose(reslog);
      }
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      // Free private workspace
      free(work);
      free(tau);
      free(VT);
      free(SIGMA);
      free(g_scr);

      // Exit point
      exit_pragma: ;

   } // End of parallel region
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   {
#if 0
      FILE *oof = fopen("nr_Q","w");
      FILE *of = fopen("QQ","w");
      int ind = 0;
      int irow = 0;
      int jcol = 0;
      for (int i = 0; i < nn; i++){
         int nr_B = iat_patt[i+1] - iat_patt[i];
         fprintf(oof,"%4d\n",nr_B);
         for (int ii = 0; ii < nr_B; ii++){
            for (int jj = 0; jj < ntv; jj++){
               fprintf(of,"%8d %8d %20.10e\n",irow+ii+1,jcol+jj+1,mat_Q[ind+jj*nr_B+ii]);
            }
         }
         irow += nr_B;
         jcol += ntv;
         ind += nr_B*ntv;
      }
      fclose(of);
      fclose(oof);
#endif
   }
   //of = fopen("P0_dopo","w");
   //for (int i = 0; i < iat_patt[nn]; i++) fprintf(of,"%20.11e\n",coef_P0[i]);
   //fclose(of);
   //cout << "FATTO" << endl;
   //exit(0);
   //fclose(bbf);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Free shared scratches
   free(c2glo);

   return ierr;

}
