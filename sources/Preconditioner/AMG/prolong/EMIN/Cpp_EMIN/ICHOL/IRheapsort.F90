subroutine IRheapsort(x1,x2,n)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: IRheapsort
!
!  Coded by Carlo Janna
!  December 2010
!
!  Purpose: Sorts an integer array x1 in such a way that x1(i) .le. x1(i+1), performing
!           the same permutation on the real double array x2
!
!  Variables:
!
!  n     : # of components of x1 and x2
!  x1(n) : array of integers
!  x2(n) : array of reals
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

! Interfaces
include 'Swap_int.h'

! Input variables
integer(IREG), intent(in)        :: n

! Input/Output variables
integer(IREG), intent(inout)     :: x1(n)
real(kind=double), intent(inout) :: x2(n)

! Local variables
integer(IREG)                    :: node,i,j,k,ik,jk
logical                          :: cont_cycle

!-----------------------------------------------------------------------------------------
do node = 2,n
   i = node
   j = i/2
   do while(x1(j).lt.x1(i))
      call swapi(x1(j),x1(i))
      call swapr(x2(j),x2(i))
      i = j
      j = i/2
      if (i.eq.1) exit
   end do
end do

do i = n,2,-1
   call swapi(x1(i),x1(1))
   call swapr(x2(i),x2(1))
   k = i-1
   ik = 1
   jk = 2
   if (k.ge.3) then
      if (x1(3).gt.x1(2))  jk = 3
   endif
   cont_cycle = .False.
   if (jk.le.k) then
      if (x1(jk).gt.x1(ik)) cont_cycle = .True.
   endif
   do while (cont_cycle)
      call swapi(x1(jk),x1(ik))
      call swapr(x2(jk),x2(ik))
      ik = jk
      jk = ik*2
      if (jk+1.le.k) then
         if (x1(jk+1).gt.x1(jk)) jk = jk+1
      endif
      cont_cycle = .False.
      if (jk.le.k) then
         if (x1(jk).gt.x1(ik)) cont_cycle = .True.
      endif
   end do
end do
!-----------------------------------------------------------------------------------------

end subroutine IRheapsort
