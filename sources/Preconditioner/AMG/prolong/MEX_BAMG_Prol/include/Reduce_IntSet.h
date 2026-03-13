#pragma once

#include "precision.h"

void Reduce_IntSet(const rExt maxrownrm, const iReg itmax_vol, const rExt tol_vol,
                   const rExt maxcond, const iReg optimal_lwork, const iReg inod,
                   const iReg n_neigh, const iReg ntvecs, iReg &row_rank,
                   const iReg *const fcnodes, iReg *int_list, const iReg *neigh,
                   const rExt *const *const TV, rExt **TVcomp, rExt *WR, rExt *coef_P,
                   rExt &row_nrm);
