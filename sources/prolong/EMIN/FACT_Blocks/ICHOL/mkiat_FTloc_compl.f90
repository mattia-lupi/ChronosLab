subroutine mkiat_FTloc_compl(nrows,nequ,nproc,firstrow,WI,iat_FT,nnz)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: mkiat_FTloc_compl
!
!  Coded by Carlo Janna
!  August 2011
!
!  Purpose: each processor compute the pointers to the beginning of each row of its part
!           of matrix FT
!
!  Variables:
!
!  nrows            : # of rows of the current processor
!  nequ             : # of equation
!  nproc            : # of computing cores
!  firstrow         : first row of the current processor
!  WI(nequ,nproc+1) : integer array s.t. WI(i,j+1) = # of F non-zeroes going in FT in
!                     i-th row and handled by j-th processor
!  iat_FT           : integer array of the pointers to the beginning of each row of
!                     matrix FT
!  nnz              : # of non-zeroes
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer(IREG), intent(in)  :: nrows,nequ,nproc,firstrow
integer(IEXT), intent(in)  :: WI(nequ,nproc+1)

! Output variables
integer(IEXT), intent(out) :: iat_FT(nrows)
integer(IEXT), intent(out) :: nnz

! Local variables
integer(IREG)              :: i,j
integer(IEXT)              :: tmp,tmp2

!-----------------------------------------------------------------------------------------

! Perform reduction
do i = 1,nrows
   iat_FT(i) = WI(firstrow-1+i,1)
enddo
do j = 2,nproc+1
   do i = 1,nrows
      iat_FT(i) = iat_FT(i) + WI(firstrow-1+i,j)
   enddo
enddo

! Compute iat_FT
tmp = iat_FT(1)
iat_FT(1) = 1
do i = 2,nrows
   tmp2 = iat_FT(i)
   iat_FT(i) = iat_FT(i-1) + tmp
   tmp = tmp2
enddo

! Compute # of non-zeroes
nnz = iat_FT(nrows) + tmp - 1

!-----------------------------------------------------------------------------------------

end subroutine mkiat_FTloc_compl
