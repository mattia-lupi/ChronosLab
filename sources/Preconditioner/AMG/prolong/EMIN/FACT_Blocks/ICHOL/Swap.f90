subroutine swapi(i1,i2)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: swapi
!
!  Coded by Carlo Janna
!  December 2010
!
!  Purpose: Swaps 2 integer variables
!
!  Variables:
!
!  i1 : integer variable
!  i2 : integer variable
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

integer(IREG), intent(inout) :: i1,i2
integer(IREG)                :: tmp

tmp = i1
i1 = i2
i2 = tmp

!-----------------------------------------------------------------------------------------

end subroutine swapi


subroutine swapr(r1,r2)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: swapr
!
!  Coded by Carlo Janna
!  December 2010
!
!  Purpose: Swaps 2 real double variables
!
!  Variables:
!
!  r1 : real double variable
!  r2 : real double variable
!
!-----------------------------------------------------------------------------------------

use class_precision

implicit none

real(kind=double), intent(inout) :: r1,r2
real(kind=double)                :: tmp

!-----------------------------------------------------------------------------------------

tmp = r1
r1 = r2
r2 = tmp

!-----------------------------------------------------------------------------------------

end subroutine swapr
