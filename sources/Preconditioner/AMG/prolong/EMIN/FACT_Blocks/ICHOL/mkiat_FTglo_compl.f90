subroutine mkiat_FTglo_compl(myid,nrows,nequ,nproc,firstrow,WI1,nnz,iat_FT)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: mkiat_FTglo_compl
!
!  Coded by Carlo Janna
!  August 2011
!
!  Purpose: compute global iat_FT and update WI! pointers
!
!  Variables:
!
!  myid     : current processor ID
!  nrows    : # of rows of myid processor
!  nequ     : # of equations
!  nproc    : # of computing cores
!  firstrow : first row of the current processor
!  WI1      : integer scratch array
!  nnz      : integer array s.t. nnz(i) = # of FT non-zeroes for each processor
!  iat_FT   : integer array of the pointers to the beginning of each row of matrix FT
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer(IREG), intent(in)    :: myid,nrows,nequ,nproc,firstrow
integer(IEXT), intent(in)    :: nnz(nproc)

! Input/Output variables
integer(IEXT), intent(inout) :: WI1(nequ,nproc+1),iat_FT(nrows+1)

! Local variables
integer(IREG)                :: i,j
integer(IEXT)                :: ntprec

!-----------------------------------------------------------------------------------------

! Compute # of non-zeros belonging to previous processors
ntprec = 0
do i = 1,myid-1
   ntprec = ntprec + nnz(i)
enddo

! Update iat_FT
do i = 1,nrows
   iat_FT(i) = iat_FT(i) + ntprec
enddo

! Account for (nrows+1)-th component for last processor
if (myid .eq. nproc) then
   iat_FT(i) = ntprec + nnz(myid) + 1
endif

! Update WI1
do i = firstrow,firstrow+nrows-1
   WI1(i,1) = iat_FT(i+1-firstrow) + WI1(i,1)
enddo
do j = 2,nproc
   do i = firstrow,firstrow+nrows-1
      WI1(i,j) = WI1(i,j) + WI1(i,j-1)
   enddo
enddo
!$OMP barrier

!-----------------------------------------------------------------------------------------

end subroutine mkiat_FTglo_compl
