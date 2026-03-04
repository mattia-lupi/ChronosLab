subroutine IcholRF_apply(nequ,nterm,iU,jU,mat_U,D_inv,vec,pvec)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: IcholRF_apply
!
!  Coded by Carlo Janna
!  February 2011
!
!  Purpose: forward and backward substitution for the application of the root-free
!           Choleski preconditioner ---> pvec:=(U^TDU)^-1*vec
!           (solves the system U^TDU*pvec=vec).
!           Note that D_inv stores the inverse of the diagonal entries
!
!  Variables:
!
!  nequ    : # of equations
!  nterm   : # of non-zero entries in U
!  iU      : topology vector of U
!  jU      : column indices of U
!  mat_U   : coefficients of U
!  D_inv   : inverse of the entries of D
!  vec     : vector to be preconditioned
!  pvec    : preconditioned vector
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer, intent(in)            :: nequ,nterm
integer, intent(in)            :: iU(nequ+1),jU(nterm)
real(kind=double), intent(in)  :: D_inv(nequ),mat_U(nterm),vec(nequ)

! Output variables
real(kind=double), intent(out) :: pvec(nequ)

! Local variables
integer                        :: k,n1,mm,i,j,m
real(kind=double)              :: a


! Initialize pvec
do i = 1,nequ
   pvec(i) = zero
enddo

! Forward substitution
do k = 1,nequ
   i = iU(k)
   j = iU(k+1) - 1
   pvec(k) = vec(k) - pvec(k)
   do m = i,j
      pvec(jU(m)) = pvec(jU(m)) + mat_U(m)*pvec(k)
   end do
end do

! Scale by D_inv
do k = 1,nequ
   pvec(k) = pvec(k)*D_inv(k)
enddo

! Backward substitution
do k = 1,nequ
   n1 = nequ-k+1
   a = zero
   i = iU(n1)
   j = iU(n1+1) - 1
   do m = i,j
      mm = j-m+i
      a = a + mat_U(mm)*pvec(jU(mm))
   end do
   pvec(n1) = pvec(n1) - a
end do

end subroutine IcholRF_apply
