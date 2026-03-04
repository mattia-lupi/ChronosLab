subroutine Complete_Matrix(myid,firstrow,nproc,nrows_L,nrows_C,nterm_L,nterm_C,           &
&                          iat_L,ja_L,iat_C,ja_C,WI1,WI2,coef_L,coef_C)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: Complete_Matrix
!
!  Coded by Carlo Janna
!  March 2018
!
!  Purpose: Parallel completion of an lower matrix L into C
!
!  Variables:
!
!  myid       : current processor ID
!  firstrow   : first row of current processor
!  nproc      : # of processors
!  nrows_L    : # of rows of L
!  nrows_C    : # of rows of C
!  nterm_L    : # of non-zeroes of L
!  nterm_C    : # of non-zeroes of C

!  [L]
!  iat_L      : topology vector of L
!  ja_L       : column indices of L
!  coef_L     : coefficients of L
!
!  [C]
!  iat_C      : topology vector of C
!  ja_C       : column indices of C
!  coef_C     : coefficients of C
!
!  WI1        : integer working array
!  WI2        : integer working array

!------------------------------------------------------------------------------------------

use class_precision

implicit none

! Interfaces
include './count_rowterms_compl_int.h'
include './mkiat_FTloc_compl_int.h'
include './mkiat_FTglo_compl_int.h'
include './mvcoef_compl_int.h'

! Input variables
integer(IREG), intent(in)  :: myid,firstrow,nproc
integer(IREG), intent(in)  :: nrows_L,nrows_C
integer(IEXT), intent(in)  :: nterm_L,nterm_C
integer(IEXT), intent(in)  :: iat_L(nrows_L+1)
integer(IREG), intent(in)  :: ja_L(nterm_L)
real(double), intent(in)   :: coef_L(nterm_L)

! Output variables
integer(IEXT), intent(out) :: iat_C(nrows_C+1)
integer(IREG), intent(out) :: ja_C(nterm_C)
integer(IEXT), intent(out) :: WI1(nrows_C,nproc+1),WI2(nproc)
real(double), intent(out)  :: coef_C(nterm_C)

! Local variables
integer(IREG)              :: i,shift,nnz
integer(IEXT)              :: ind,ind2
integer(IREG)              :: nt

!------------------------------------------------------------------------------------------

! Counts # of L non-zeroes assigned to each processor in each C row
call count_rowterms_compl(nrows_C,nterm_L,ja_L,WI1(1:nrows_C,myid+1))
shift = firstrow - 1
do i = 1,nrows_L
   ! Store the original number of non-zeroes of the lower part
   nnz = iat_L(i+1) - iat_L(i) - 1
   WI1(shift+i,1) = nnz
enddo
!$OMP barrier

! Each processor computes its part of iat_C
call mkiat_FTloc_compl(nrows_L,nrows_C,nproc,firstrow,WI1,iat_C(firstrow),WI2(myid))
!$OMP barrier

! Update iat_C and WI1 pointers
call mkiat_FTglo_compl(myid,nrows_L,nrows_C,nproc,firstrow,WI1,WI2,iat_C(firstrow))
!$OMP barrier

! Transpose coefficients of L into the upper part of C
call mvcoef_compl(firstrow,nrows_L,nrows_C,nterm_C,iat_L,ja_L,ja_C,coef_L,coef_C,         &
&               WI1(1:nrows_C,myid))
!$OMP barrier

! Copy the lower part of L into C
do i = 1,nrows_L
   ind = iat_L(i)
   nt = iat_L(i+1) - 1 - ind
   ind2 = iat_C(shift+i)
   if (nt > 0) then
      call SCOPY(nt,ja_L(ind),1,ja_C(ind2),1)
      call DCOPY(nt,coef_L(ind),1,coef_C(ind2),1)
   endif
enddo

!-----------------------------------------------------------------------------------------

end subroutine Complete_Matrix
