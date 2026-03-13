#include "precision.h"

// THIS WRAPPER NEEDS ALSO A DESTRUCTOR TO DEALLOCATE FORTRAN POINTERS
extern "C" {
   iReg free_fptr_iReg( iReg *c_ptr_in , iExt len_in);
   iReg free_fptr_iExt( iExt *c_ptr_in , iExt len_in);
   iReg free_fptr_rExt( rExt *c_ptr_in , iExt len_in);
}
