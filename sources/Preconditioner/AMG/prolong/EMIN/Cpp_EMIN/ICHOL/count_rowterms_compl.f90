subroutine count_rowterms_compl(nequ,nterm,ja,WI)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: count_rowterms_compl
!
!  Coded by Carlo Janna
!  August 2011
!
!  Purpose: Counts # of F non-zeroes assigned to myid processor in each FT row
!
!  Variables:
!
!  nequ     : # of equations
!  nterm    : # of non-zeroes of F of myid processor
!  ja       : integer array of the column indices of matrix F
!  WI(nequ) : integer array s.t. WI(i) = # of F non-zeroes of myid in FT i-th row
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer(IREG), intent(in)  :: nequ
integer(IEXT), intent(in)  :: nterm
integer(IREG), intent(in)  :: ja(nterm)

! Output variables
integer(IEXT), intent(out) :: WI(nequ)

! Local variables
integer(IEXT)              :: i

!-----------------------------------------------------------------------------------------

! Initialize WI
do i = 1,nequ
   WI(i) = 0
end do

! Count non-zeroes
do i = 1,nterm
   WI(ja(i)) = WI(ja(i)) + 1
enddo

!-----------------------------------------------------------------------------------------

end subroutine count_rowterms_compl
