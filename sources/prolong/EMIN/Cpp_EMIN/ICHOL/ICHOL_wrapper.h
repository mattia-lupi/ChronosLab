#include "precision.h"

// THIS WRAPPER NEEDS ALSO A DESTRUCTOR TO DEALLOCATE FORTRAN POINTERS
extern "C" {
   iReg ICHOL_wrapper(const iReg lfil, const iReg jcol_offset, const iReg nn_in,
                      const iExt nt_in, const iExt ireg_scr_size, const iExt iext_scr_size,
                      const iExt *iat_in, const iReg *ja_in, const rExt *coef_in, iExt *iU,
                      iReg *jU, rExt *mat_U, rExt *D_inv,
                      iReg *ireg_scr_in, iExt *iext_scr_in);
}
