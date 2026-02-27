module class_resize
!*****************************************************************************************
!
!  Class: resize
!
!  Coded by Carlo Janna
!  December 2011
!
!  Purpose: resize a given input array (given as a pointer), returns info .ne. 0 if some
!           error occurs
!
!*****************************************************************************************

use class_precision

implicit none

public :: resize_r,resize_ireg

contains

!-----------------------------------------------------------------------------------------

    function resize_r(n,vec) result(info)

    implicit none

!   Input variables
    integer(IEXT), intent(in) :: n
    real(double), pointer     :: vec(:)
!   Output variables
    integer(IREG)             :: info
!   Local variables
    integer(IEXT)             :: n_old
    real(double), pointer     :: tmp(:),tmp1(:)

    info = 0
    n_old = size(vec)
    if (n_old .ne. n) then
!      Allocate a temporary array
       allocate(tmp(n),stat=info)
       if (info.ne.0) return
!      Swap pointers
       tmp1 => vec
       vec  => tmp
!      Copy old entries in the new array
       call DCOPY(min(n,n_old),tmp1,1,vec,1)
!      Deallocate old array
       deallocate(tmp1,stat=info)
       tmp => null()
    endif

    end function resize_r

!-----------------------------------------------------------------------------------------

    function resize_ireg(n,vec) result(info)

    implicit none

!   Input variables
    integer(IEXT), intent(in) :: n
    integer(IREG), pointer    :: vec(:)
!   Output variables
    integer(IREG)             :: info
!   Local variables
    integer(IEXT)             :: n_old
    integer(IREG), pointer    :: tmp(:),tmp1(:)

    info = 0
    n_old = size(vec)
    if (n_old .ne. n) then
!      Allocate a temporary array
       allocate(tmp(n),stat=info)
       if (info.ne.0) return
!      Swap pointers
       tmp1 => vec
       vec  => tmp
!      Copy old entries in the new array
       call SCOPY(min(n,n_old),tmp1,1,vec,1)
!      Deallocate old array
       deallocate(tmp1,stat=info)
       tmp => null()
    endif

    end function resize_ireg

!-----------------------------------------------------------------------------------------

end module class_resize
