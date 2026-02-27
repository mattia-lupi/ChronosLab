subroutine RI_SortSplit(n,ncut,R_vec,I_vec)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: RI_SortSplit
!
!  Coded by Carlo Janna
!  December 2010
!
!  Purpose: Sorts the real vector R_vec such that:
!           abs(R_vec(i)) .ge. abs(R_vec(ncut)) for i .lt. ncut and
!           abs(R_vec(i)) .le. abs(R_vec(ncut)) for i .gt. ncut
!           And perfroms the same permutations on integer vector I_vec
!
!  Variables:
!
!  n     : # of R_vec
!  ncut  : # of largest R_vec entries to be sorted
!  R_vec : real vector to be sorted and splitted
!  I_vec : auxiliary vector to be sorted and splitted as R_vec
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Input variables
integer(IREG), intent(in)        :: n, ncut

! Input/Ouput variables
real(kind=double), intent(inout) :: R_vec(n)
integer(IREG), intent(inout)     :: I_vec(n)

! Local variables
integer(IREG)                    :: itmp, first, last, mid, j
real(kind=double)                :: tmp, absval

!-----------------------------------------------------------------------------------------

first = 1
last = n

! Check ncut consistency
if (ncut .lt. first .or. ncut .gt. last) return

outer_loop: do
   mid = first
   absval = abs(R_vec(mid))
   inner_loop: do j=first+1, last
      if (abs(R_vec(j)) .gt. absval) then
         mid = mid+1
!        Exchange
         tmp = R_vec(mid)
         itmp = I_vec(mid)
         R_vec(mid) = R_vec(j)
         I_vec(mid) = I_vec(j)
         R_vec(j)  = tmp
         I_vec(j) = itmp
      endif
   enddo inner_loop
!  Exchange
   tmp = R_vec(mid)
   R_vec(mid) = R_vec(first)
   R_vec(first)  = tmp
   itmp = I_vec(mid)
   I_vec(mid) = I_vec(first)
   I_vec(first) = itmp
!  Exit test
   if (mid .eq. ncut) return
   if (mid .gt. ncut) then
      last = mid-1
   else
      first = mid+1
   endif
enddo outer_loop

!-----------------------------------------------------------------------------------------

end subroutine RI_SortSplit
