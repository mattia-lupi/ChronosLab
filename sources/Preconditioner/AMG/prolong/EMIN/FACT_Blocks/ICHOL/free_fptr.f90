function dealloc_fptr_iReg(c_ptr_in,len_in) result(ierr) bind ( C, name="free_fptr_iReg" )

use class_precision
use iso_c_binding, only: c_ptr, c_f_pointer

implicit none

! Input variables
type(C_PTR), value, intent(in) :: c_ptr_in
integer(CB_IEXT), value        :: len_in

! Result
integer(CB_IREG)               :: ierr

! Locals
integer(CB_IREG), pointer      :: f_ptr(:)

! Convert C pointers into fortran pointers
call c_f_pointer(c_ptr_in,f_ptr,[len_in])

deallocate(f_ptr,stat=ierr)

end function dealloc_fptr_iReg

!-----------------------------------------------------------------------------------------

function dealloc_fptr_iExt(c_ptr_in,len_in) result(ierr) bind ( C, name="free_fptr_iExt" )

use class_precision
use iso_c_binding, only: c_ptr, c_f_pointer

implicit none

! Input variables
type(C_PTR), value, intent(in) :: c_ptr_in
integer(CB_IEXT), value        :: len_in

! Result
integer(CB_IREG)               :: ierr

! Locals
integer(CB_IEXT), pointer      :: f_ptr(:)

! Convert C pointers into fortran pointers
call c_f_pointer(c_ptr_in,f_ptr,[len_in])

deallocate(f_ptr,stat=ierr)

end function dealloc_fptr_iExt

!-----------------------------------------------------------------------------------------

function dealloc_fptr_rExt(c_ptr_in,len_in) result(ierr) bind ( C, name="free_fptr_rExt" )

use class_precision
use iso_c_binding, only: c_ptr, c_f_pointer

implicit none

! Input variables
type(C_PTR), value, intent(in) :: c_ptr_in
integer(CB_IEXT), value        :: len_in

! Result
integer(CB_IREG)               :: ierr

! Locals
real(CB_DBL), pointer          :: f_ptr(:)

! Convert C pointers into fortran pointers
call c_f_pointer(c_ptr_in,f_ptr,[len_in])

deallocate(f_ptr,stat=ierr)

end function dealloc_fptr_rExt
