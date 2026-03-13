subroutine mvcoef_compl(firstrow,nrows,nequ,nterm_F,iat_F,ja_F,ja_FT,mat_F,mat_FT,punt)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: mvcoef_compl
!
!  Coded by Carlo Janna
!  August 2011
!
!  Purpose: transpose matrix F coefficients in matrix FT
!
!  Variables:
!
!  firstrow : first row of the current processor
!  nrows    : # of rows of the current processor
!  nequ     : # of equations
!  nterm_F  : # of non-zeroes of F for the current processor
!
!  [F]
!  iat_F    : topology vector of F
!  ja_F     : column indices of F
!  mat_F    : coefficients of F
!
!  [FT]
!  ja_FT    : column indices of FT
!  mat_FT   : coefficients of FT
!
!  punt     : integer scratch array
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer(IREG), intent(in)    :: firstrow,nrows,nequ
integer(IEXT), intent(in)    :: nterm_F
integer(IREG), intent(in)    :: ja_F(nterm_F)
integer(IEXT), intent(in)    :: iat_F(nrows+1)
real(double), intent(in)     :: mat_F(nterm_F)

! Input/Output variables
integer(IEXT), intent(inout) :: punt(nequ)

! Output variables
integer(IREG), intent(out)   :: ja_FT(nterm_F)
real(double), intent(out)    :: mat_FT(nterm_F)

! Local variables
integer(IREG)                :: i,irow,shift
integer(IEXT)                :: j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer omp_get_thread_num,myid
!myid = omp_get_thread_num() + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------------------------

! Transpose local stripe of F
shift = firstrow - 1
do i = 1,nrows
   do j = iat_F(i),iat_F(i+1)-1
      irow = ja_F(j)
      ja_FT(punt(irow)) = i + shift
      mat_FT(punt(irow)) = mat_F(j)
      punt(irow) = punt(irow) + 1
   enddo
enddo

!-----------------------------------------------------------------------------------------

end subroutine mvcoef_compl
