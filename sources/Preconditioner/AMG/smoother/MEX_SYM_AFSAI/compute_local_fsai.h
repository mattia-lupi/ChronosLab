//----------------------------------------------------------------------------------------

#pragma once
#include "cpt_afsai_coef.h"

//----------------------------------------------------------------------------------------

void compute_local_fsai(iReg nthread, iReg n_step, iReg step_size, rExt tau, rExt eps,
                        iReg nrows, iReg nrows_M, iExt nterm_M, const iExt * const iat_M,
                        const iReg * const ja_M, const rExt * const coef_M, iExt *nterm_G,
                        iExt *iat_G, iReg *ja_G, rExt *coef_G);



