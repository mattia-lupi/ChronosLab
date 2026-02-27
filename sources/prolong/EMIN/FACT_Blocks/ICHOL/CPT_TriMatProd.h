#include "precision.h"

// THIS WRAPPER NEEDS ALSO A DESTRUCTOR TO DEALLOCATE FORTRAN POINTERS
extern "C" {
iReg CPT_TriMatProd(const iReg lfil,const iReg nthreads,const iReg nrows_F1,
                    const iReg ncols_F1, iExt *c_iat_F1, iReg *c_ja_F1,
                    const rExt *c_coef_F1, const iReg nrows_A, const iReg ncols_A,
                    iExt *c_iat_A, iReg *c_ja_A, const rExt *c_coef_A,
                    const iReg nrows_F2, const iReg ncols_F2, iExt *c_iat_F2,
                    iReg *c_ja_F2, const rExt *c_coef_F2, iReg &nrows_S,
                    iReg &ncols_S, iExt *&c_iat_S, iReg *&c_ja_S,
                    rExt *&c_coef_S);
}
