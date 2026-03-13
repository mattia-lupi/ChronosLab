function ICHOL_wrapper(lfil,jcol_offset,nn_in,nt_in,ireg_scr_size,iext_scr_size,          &
                       iat_in,ja_in,coef_in,iU_in,jU_in,mat_U_in,D_inv_in,                &
                       ireg_scr_in,iext_scr_in) result(ierr)                              &
                       bind ( C, name="ICHOL_wrapper" )
!-----------------------------------------------------------------------------------------
!
!  Subroutine: ICHOL_wrapper
!
!  Coded by Carlo Janna
!  November 2020
!
!  Purpose: Computes an Incomplete Cholesky factorization of an input matrix.
!
!  NOTE: the input and output matrices are in FORTRAN style.
!        jcol_offset is used only for debug purposes currently
!
!  Variables:
!
!  lfil:        number of entries retained for each row of U in addition to those in A
!  jcol_offset: column offset of local part of A with respect to global
!
!  ierr:        error code: 0 ---> success
!                           1 ---> not enough room in the scratches ireg_scr and iext_scr
!                           2 ---> zero or negative pivot encountered during factorization
!
!-----------------------------------------------------------------------------------------

use class_precision
use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

implicit none

include 'Ichol_RootFree_shift_int.h'

! Input variables
integer(CB_IREG), value, intent(in)    :: lfil
integer(CB_IREG), value, intent(in)    :: jcol_offset
integer(CB_IREG), value, intent(in)    :: nn_in
integer(CB_IEXT), value, intent(in)    :: nt_in
integer(CB_IEXT), value, intent(in)    :: ireg_scr_size
integer(CB_IEXT), value, intent(in)    :: iext_scr_size
type(C_PTR), value, intent(in)         :: iat_in
type(C_PTR), value, intent(in)         :: ja_in
type(C_PTR), value, intent(in)         :: coef_in
type(C_PTR), value, intent(in)         :: iU_in
type(C_PTR), value, intent(in)         :: jU_in
type(C_PTR), value, intent(in)         :: mat_U_in,D_inv_in
type(C_PTR), value, intent(in)         :: ireg_scr_in
type(C_PTR), value, intent(in)         :: iext_scr_in

! Result
integer(CB_IREG)                       :: ierr

! Pointers to the output
integer(IEXT), pointer                 :: it_U(:)
integer(IREG), pointer                 :: jcol_U(:)
real(double), pointer                  :: mat_U(:),D_inv(:)

! Locals
integer(IREG)                          :: i,ineg
integer(IEXT)                          :: j,istrt,iwk_U
real(double)                           :: droptol,alpha,small
integer(IEXT), pointer                 :: iat_A(:)
integer(IREG), pointer                 :: ja_A(:)
real(double), pointer                  :: coef_A(:)
integer(IREG), pointer                 :: ireg_scr(:)
integer(IEXT), pointer                 :: iext_scr(:)
integer(IREG), pointer                 :: U_irow_lst(:)
integer(IEXT), pointer                 :: U_nxt_lst(:),U_end_lst(:),ind_U(:),JWI(:)

!-----------------------------------------------------------------------------------------

! Init error code
ierr = 0

! Compute ICHOL size
iwk_U = nt_in + nn_in * (lfil+1)

! Convert C pointers into fortran pointers
call c_f_pointer(iat_in,iat_A,[nn_in+1])
call c_f_pointer(ja_in,ja_A,[nt_in])
call c_f_pointer(coef_in,coef_A,[nt_in])
call c_f_pointer(iU_in,it_U,[nn_in+1])
call c_f_pointer(jU_in,jcol_U,[iwk_U])
call c_f_pointer(mat_U_in,mat_U,[iwk_U])
call c_f_pointer(D_inv_in,D_inv,[nn_in])
call c_f_pointer(ireg_scr_in,ireg_scr,[ireg_scr_size])
call c_f_pointer(iext_scr_in,iext_scr,[iext_scr_size])

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (jcol_offset == 1) then
!   write(*,*) 'DENTRO TAGLIA',nn_in
!   do i = 1,nn_in
!      do j = iat_A(i),iat_A(i+1)-1
!         write(777,'(2i10,e20.11)') i,ja_A(j),coef_A(j)
!      enddo
!   enddo
!endif
!call flush(777)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set some parameters for preconditioner computation
droptol = 1.d-10
droptol = ZERO
alpha = ZERO
small = ZERO

! Set internal pointers to scratches
if ( iext_scr_size < (iwk_U + 4*nn_in) .or. ireg_scr_size < (iwk_U+nn_in) ) then
   ierr = 1
   return
endif
U_irow_lst =>  ireg_scr(1:iwk_U+nn_in)
istrt = 1
U_nxt_lst => iext_scr(istrt:iwk_U+nn_in)
istrt = istrt + iwk_U+nn_in
U_end_lst => iext_scr(istrt:istrt-1+nn_in)
istrt = istrt + nn_in
ind_U => iext_scr(istrt:istrt-1+nn_in)
istrt = istrt + nn_in
JWI => iext_scr(istrt:istrt-1+nn_in)

! Compute Incomplete Cholesky Factorization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ineg = 0
!write(*,*) 'CHIAMO Ichol_RootFree_shift: ',nn_in,nt_in,lfil,iwk_U
!do i = 1,nn_in+1
!   it_U(i) = 1
!enddo
!write(*,*) 'TOCCATO it_U'
!do i = 1,iwk_U
!   jcol_U(i) = 1
!enddo
!write(*,*) 'TOCCATO jcol_U'
!do i = 1,iwk_U
!   mat_U(i) = 2.0
!enddo
!write(*,*) 'TOCCATO mat_U'
!do i = 1,nn_in
!   D_inv(i) = 2.0
!enddo
!write(*,*) 'TOCCATO D_inv'
!do i = 1,iwk_U+nn_in
!   U_irow_lst(i) = 3
!enddo
!write(*,*) 'TOCCATO D_inv'
!do i = 1,iwk_U+nn_in
!   U_nxt_lst(i) = 3
!enddo
!write(*,*) 'TOCCATO U_nxt_lst'
!do i = 1,nn_in
!   U_end_lst(i) = 3
!enddo
!write(*,*) 'TOCCATO U_end_lst'
!do i = 1,nn_in
!   ind_U(I) = 3
!enddo
!write(*,*) 'TOCCATO ind_U'
!do i = 1,nn_in
!   JWI(i) = 3
!enddo
!write(*,*) 'TOCCATO JWI'
!!if (nn_in == 139 .and. nt_in == 1622) then
!if (nn_in == 196 .and. nt_in == 2682) then
!   write(111,*) nn_in,nt_in
!   do i = 1,nn_in
!      do j = iat_A(i),iat_A(i+1)-1
!         write(111,'(2i6,e15.6)') i,ja_A(j),coef_A(j)
!      enddo
!   enddo
!   call flush(111)
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call Ichol_RootFree_shift(nn_in,nt_in,lfil,droptol,alpha,small,iat_A,ja_A,coef_A,iwk_U,   &
&                         it_U,jcol_U,mat_U,D_inv,U_irow_lst,U_nxt_lst,U_end_lst,ind_U,   &
&                         JWI,ineg,ierr)
if (ierr /= 0) then
   write(*,'(a,i2,a,i5,a,i5)') 'ERROR ',ierr,' IN CHOLESKY FACT. NEGATIVE DIAG: ',ineg,   &
                               ' OUT OF ',nn_in
   ierr = 2
   return
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!if (jcol_offset == 1) then
!   write(222,*) nn_in,it_U(nn_in+1)-1
!   do i = 1,nn_in
!      do j = it_U(i),it_U(i+1)-1
!         write(222,'(2i6,e15.6)') i,jcol_U(j),mat_U(j)
!      enddo
!   enddo
!   write(555,*) it_U(1:nn_in+1)
!   do i = 1,nn_in
!      write(333,'(2i6,e15.6)') i,i,1.d0/D_inv(i)
!   enddo
!   call flush(222)
!   call flush(333)
!endif
!write(*,*) 'USCITO DA Ichol_RootFree_shift: ',ierr,ineg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

end function ICHOL_wrapper
