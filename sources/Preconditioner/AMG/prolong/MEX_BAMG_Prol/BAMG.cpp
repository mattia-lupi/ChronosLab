#include <omp.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <math.h>
#include <lapacke.h>

#include "precision.h"
#include "DebEnv.h"
#include "BAMG_params.h"
#include "linsol_error.h"
#include "check_neigh_cond.h"
#include "add_new_neighs.h"
#include "inl_blas1.h"
#include "maxVol.h"
#include "heapsort.h"
#include "Reduce_IntSet.h"

int BAMG(const BAMG_params &params, const iReg nthreads, iReg nn_L, iReg nn_C,
         iReg nn_S, const iExt *const iat_S, const iReg *const ja_S, const iReg ntv,
         const iReg *const fcnodes, const rExt *const *const TV,
         iExt &nt_I, std::vector<iExt> &vec_iat_I, std::vector<iReg> &vec_ja_I,
         std::vector<rExt> &vec_coef_I, std::vector<iReg> &vec_c_mark) {

   const rExt ONE = 1.0;
   DebEnv.nthreads = nthreads;

   try {
      vec_c_mark.assign(nn_C, 0);
   } catch (const std::exception &) {
      linsol_error("BAMG", "allocating c_mark");
      return 1;
   }
   iReg *c_mark = vec_c_mark.data();

   const iReg len_count = params.dist_max + 4;
   std::vector<iReg> dist_count;
   if (params.verbosity >= VLEV_MEDIUM) {
      dist_count.assign(len_count, 0);
   }

   const iReg CHUNK = 32;
   const iReg nchunks = (nn_C > 0) ? ((nn_C + CHUNK - 1) / CHUNK) : 0;
   std::vector<std::vector<iExt>> chunk_iat(nchunks);
   std::vector<std::vector<iReg>> chunk_ja(nchunks);
   std::vector<std::vector<rExt>> chunk_coef(nchunks);
   std::vector<iExt> chunk_nt(nchunks, 0);

   int ierr = 0;

   // Dynamic 32-row chunks balance variable BFS/maxVol work across irregular fine rows
   #pragma omp parallel num_threads(nthreads)
   {
      int ierr_L = 0;
      std::vector<iReg> dist_count_local;
      if (params.verbosity >= VLEV_MEDIUM) {
         dist_count_local.assign(len_count, 0);
      }

      // Hoist LAPACK workspace query and thread-local scratch once per thread
      iReg optimal_lwork = 1;
      if (ntv > 0) {
         rExt db_lwork = 0.0;
         lapack_int qinfo = LAPACKE_dgels_work(LAPACK_COL_MAJOR, 'N', ntv, ntv, 1,
                                               nullptr, ntv, nullptr, ntv, &db_lwork, -1);
         if (qinfo != 0) {
            ierr_L = 2;
         } else {
            optimal_lwork = std::max(static_cast<iReg>(1), static_cast<iReg>(db_lwork));
         }
      }

      std::vector<iReg> vec_int_list(nn_S), vec_neigh(nn_S), vec_WI(nn_S, 0);
      std::vector<rExt> vec_WR(optimal_lwork), vec_wrot(ntv), vec_vcol(nn_S), vec_vrow(ntv);
      std::vector<rExt*> vec_TVcomp(nn_S + 1);
      std::vector<rExt>  vec_TVbuf(static_cast<std::size_t>(ntv) * (nn_S + 1));
      for (iReg i = 0; i < nn_S + 1; ++i) {
         vec_TVcomp[i] = &vec_TVbuf[static_cast<std::size_t>(i) * ntv];
      }

      iReg *int_list = vec_int_list.data();
      iReg *neigh    = vec_neigh.data();
      iReg *WI       = vec_WI.data();
      rExt *WR       = vec_WR.data();
      rExt *wrot_buf = vec_wrot.data();
      rExt *vcol_buf = vec_vcol.data();
      rExt *vrow_buf = vec_vrow.data();
      rExt **TVcomp  = vec_TVcomp.data();

      const bool VERB_FLAG = params.verbosity >= VLEV_MEDIUM;
      const iReg itmax_vol = params.itmax_vol;
      const iReg dist_min  = params.dist_min;
      const iReg dist_max  = params.dist_max;
      const iReg mmax      = params.mmax;
      const rExt maxcond   = params.maxcond;
      const rExt maxrownrm = params.maxrownrm;
      const rExt tol_vol   = params.tol_vol;
      const rExt eps       = params.eps;

      #pragma omp for schedule(dynamic, 1)
      for (iReg c = 0; c < nchunks; ++c) {
         if (ierr_L != 0) continue;
         const iReg c_start = nn_L + c * CHUNK;
         const iReg c_end   = std::min(nn_L + nn_C, c_start + CHUNK);
         const iReg c_len   = c_end - c_start;

         auto &loc_iat  = chunk_iat[c];
         auto &loc_ja   = chunk_ja[c];
         auto &loc_coef = chunk_coef[c];
         loc_iat.resize(c_len + 1);
         const iExt cap = std::max(params.mmax, ntv) * c_len;
         loc_ja.resize(cap);
         loc_coef.resize(cap);

         iExt ind_P = 0;
         loc_iat[0] = 0;

         try {
            for (iReg inod = c_start; inod < c_end; ++inod) {
               iReg inod_coarse = fcnodes[inod];
               if (inod_coarse >= 0) {
                  loc_ja[ind_P] = inod_coarse;
                  loc_coef[ind_P] = ONE;
                  ind_P++;
               } else {
                  iReg row_rank = 0;
                  iReg distance = 0;
                  iReg n_neigh = 1;
                  neigh[0] = inod;
                  WI[inod] = 1;

                  iReg istart_neigh = 0, iend_neigh = 1;
                  while (distance < dist_min) {
                     distance++;
                     add_new_neighs(iat_S, ja_S, istart_neigh, iend_neigh, nn_S, n_neigh, neigh, WI);
                     if (n_neigh == iend_neigh) break;
                     istart_neigh = iend_neigh;
                     iend_neigh = n_neigh;
                  }

                  rExt bnorm = inl_dnrm2(ntv, TV[inod], 1);
                  rExt res_ass = eps * bnorm;
                  rExt res_row = 0.0, row_nrm = 0.0;
                  bool UPD_WEIGHTS = true;

                  while (UPD_WEIGHTS) {
                     iReg n_int = 0;
                     for (iReg i = 1; i < n_neigh; ++i) {
                        iReg i_neigh = neigh[i];
                        if (fcnodes[i_neigh] >= 0) {
                           int_list[n_int] = i_neigh;
                           for (iReg j = 0; j < ntv; ++j) TVcomp[n_int][j] = TV[i_neigh][j];
                           n_int++;
                        }
                     }
                     if (n_int > 0) {
                        maxVol(mmax, maxcond, itmax_vol, tol_vol, n_int, ntv, TVcomp, row_rank,
                               int_list, wrot_buf, vcol_buf, vrow_buf);
                        heapsort(int_list, row_rank);
                        for (iReg i = 0; i < row_rank; ++i) {
                           iReg i_neigh = int_list[i];
                           for (iReg j = 0; j < ntv; ++j) TVcomp[i][j] = TV[i_neigh][j];
                        }
                        for (iReg j = 0; j < ntv; ++j) loc_coef[ind_P + j] = TV[inod][j];
                        lapack_int info = LAPACKE_dgels_work(LAPACK_COL_MAJOR, 'N', ntv, row_rank, 1,
                                                             &(TVcomp[0][0]), ntv, &(loc_coef[ind_P]), ntv,
                                                             WR, optimal_lwork);
                        if (info != 0) throw linsol_error("BAMG", "error in LAPACKE_dgels");
                        res_row = inl_dnrm2(ntv - row_rank, &(loc_coef[ind_P + row_rank]), 1);
                        row_nrm = inl_dnrm2(row_rank, &(loc_coef[ind_P]), 1);
                     } else {
                        row_rank = 0;
                        res_row = bnorm;
                        row_nrm = 0.0;
                     }

                     if (row_nrm < maxrownrm && ((res_row < res_ass) || (row_rank == ntv))) {
                        UPD_WEIGHTS = false;
                        if (VERB_FLAG && distance > 0) dist_count_local[distance - 1]++;
                     } else if (distance < dist_max) {
                        distance++;
                        add_new_neighs(iat_S, ja_S, istart_neigh, iend_neigh, nn_S, n_neigh, neigh, WI);
                        if (n_neigh == iend_neigh) {
                           UPD_WEIGHTS = false;
                           if (n_int == 0) {
                              c_mark[inod - nn_L] = 1;
                              if (VERB_FLAG) dist_count_local[dist_max + 2]++;
                           } else {
                              if (row_nrm > RELAX_FAC * maxrownrm) {
                                 Reduce_IntSet(RELAX_FAC * maxrownrm, itmax_vol, tol_vol, maxcond,
                                               optimal_lwork, inod, n_neigh, ntv, row_rank,
                                               fcnodes, int_list, neigh, TV, TVcomp, WR,
                                               &(loc_coef[ind_P]), row_nrm, wrot_buf, vcol_buf, vrow_buf);
                              }
                              res_row = (row_rank == 0) ? bnorm : inl_dnrm2(ntv - row_rank, &(loc_coef[ind_P + row_rank]), 1);
                              if (res_row < res_ass) {
                                 if (VERB_FLAG) dist_count_local[dist_max + 1]++;
                              } else {
                                 c_mark[inod - nn_L] = 1;
                                 if (VERB_FLAG) dist_count_local[dist_max]++;
                              }
                           }
                        }
                        istart_neigh = iend_neigh;
                        iend_neigh = n_neigh;
                     } else {
                        UPD_WEIGHTS = false;
                        if (row_nrm > RELAX_FAC * maxrownrm) {
                           Reduce_IntSet(RELAX_FAC * maxrownrm, itmax_vol, tol_vol, maxcond,
                                         optimal_lwork, inod, n_neigh, ntv, row_rank, fcnodes,
                                         int_list, neigh, TV, TVcomp, WR, &(loc_coef[ind_P]),
                                         row_nrm, wrot_buf, vcol_buf, vrow_buf);
                        }
                        res_row = (row_rank == 0) ? bnorm : inl_dnrm2(ntv - row_rank, &(loc_coef[ind_P + row_rank]), 1);
                        if (res_row > res_ass) c_mark[inod - nn_L] = 1;
                        if (VERB_FLAG) dist_count_local[dist_max + (res_row > res_ass ? 0 : 1)]++;
                     }
                  }
                  if (VERB_FLAG && row_rank < ntv) dist_count_local[dist_max + 3]++;
                  for (iReg i = 0; i < row_rank; ++i) {
                     loc_ja[ind_P] = fcnodes[int_list[i]];
                     ind_P++;
                  }
                  for (iReg i = 0; i < n_neigh; ++i) WI[neigh[i]] = 0;
               }
               loc_iat[(inod - c_start) + 1] = ind_P;
            }
            chunk_nt[c] = ind_P;
         } catch (const std::exception &) {
            linsol_error("BAMG", "in computing BAMG chunk");
            ierr_L = 2;
         }
      }

      if (VERB_FLAG) {
         #pragma omp critical (dist_count_reduce)
         {
            for (iReg i = 0; i < len_count; ++i) dist_count[i] += dist_count_local[i];
         }
      }

      #pragma omp atomic update
      ierr += ierr_L;
   }

   if (ierr != 0) return ierr;

   std::vector<iExt> chunk_off(nchunks + 1, 0);
   for (iReg c = 0; c < nchunks; ++c) chunk_off[c + 1] = chunk_off[c] + chunk_nt[c];
   nt_I = chunk_off[nchunks];

   try {
      vec_iat_I.resize(nn_S + 1);
      vec_ja_I.resize(nt_I);
      vec_coef_I.resize(nt_I);
   } catch (const std::exception &) {
      linsol_error("BAMG", "allocating prolongation matrix");
      return 1;
   }

   iExt *iat_I  = vec_iat_I.data();
   iReg *ja_I   = vec_ja_I.data();
   rExt *coef_I = vec_coef_I.data();

   #pragma omp parallel for schedule(static) num_threads(nthreads)
   for (iReg c = 0; c < nchunks; ++c) {
      const iReg c_start = nn_L + c * CHUNK;
      const iReg c_end   = std::min(nn_L + nn_C, c_start + CHUNK);
      const iReg c_len   = c_end - c_start;
      iExt dst = chunk_off[c];
      const auto &loc_iat  = chunk_iat[c];
      const auto &loc_ja   = chunk_ja[c];
      const auto &loc_coef = chunk_coef[c];
      for (iReg r = 0; r < c_len; ++r) {
         iat_I[c_start + r] = dst;
         const iExt s0 = loc_iat[r], s1 = loc_iat[r + 1], len = s1 - s0;
         for (iExt k = 0; k < len; ++k) {
            ja_I[dst + k] = loc_ja[s0 + k];
            coef_I[dst + k] = loc_coef[s0 + k];
         }
         dst += len;
      }
   }
   for (iReg i = 0; i < nn_L; ++i) iat_I[i] = 0;
   for (iReg i = nn_L + nn_C; i <= nn_S; ++i) iat_I[i] = nt_I;

   if (params.verbosity >= VLEV_MEDIUM) {
      iReg tot_count = 0;
      for (iReg i = 0; i <= params.dist_max + 2; ++i) tot_count += dist_count[i];
      rExt r_tot_count = (tot_count > 0) ? (100.0 / static_cast<rExt>(tot_count)) : 0.0;
      for (iReg i = 0; i < params.dist_max; ++i)
         std::cout << "# of nodes interp at dist " << std::setw(4) << i + 1 << ":  " <<
         std::setw(12) << dist_count[i] << " | " << std::setprecision(2) << std::setw(6) <<
         static_cast<rExt>(dist_count[i]) * r_tot_count << "%" << std::endl;
      std::cout << "# of nodes with high error:      " <<
      std::setw(12) << dist_count[params.dist_max] << " | " << std::setprecision(2) <<
      std::setw(6) << static_cast<rExt>(dist_count[params.dist_max]) * r_tot_count <<
      "%" << std::endl;
      std::cout << "# of nodes with large weights:   " << std::setw(12) <<
      dist_count[params.dist_max + 1] << " | " << std::setprecision(2) << std::setw(6) <<
      static_cast<rExt>(dist_count[params.dist_max + 1]) * r_tot_count << "%" << std::endl;
      std::cout << "# of nodes without neighbours:   " << std::setw(12) <<
      dist_count[params.dist_max + 2] << " | " << std::setprecision(2) << std::setw(6) <<
      static_cast<rExt>(dist_count[params.dist_max + 2]) * r_tot_count << "%" << std::endl;
      std::cout << "---------------------------" << std::endl;
      std::cout << "NOT FULL-RANK " << std::setprecision(2) << std::setw(6)
           << static_cast<rExt>(dist_count[params.dist_max + 3]) * r_tot_count << "%" << std::endl;
   }

   return 0;
}
