! Conversion of indices from C to FORTRAN and from FORTRAN to C
module C2FOR

use class_precision
implicit none

contains

   !--------------------------------------------------------------------------------------

   subroutine idx_c2for_IREG(nn,ivec)

   implicit none
   integer(IEXT), intent(in)    :: nn
   integer(IREG), intent(inout) :: ivec(nn)
   integer(IEXT)                :: i

   do i = 1,nn
      ivec(i) = ivec(i) + 1
   enddo

   end subroutine idx_c2for_IREG

   !--------------------------------------------------------------------------------------

   subroutine idx_c2for_IEXT(nn,ivec)

   implicit none
   integer(IEXT), intent(in)    :: nn
   integer(IEXT), intent(inout) :: ivec(nn)
   integer(IEXT)                :: i

   do i = 1,nn
      ivec(i) = ivec(i) + 1
   enddo

   end subroutine idx_c2for_IEXT

   !--------------------------------------------------------------------------------------

   subroutine idx_for2c_IREG(nn,ivec)

   implicit none
   integer(IEXT), intent(in)    :: nn
   integer(IREG), intent(inout) :: ivec(nn)
   integer(IEXT)                :: i

   do i = 1,nn
      ivec(i) = ivec(i) - 1
   enddo

   end subroutine idx_for2c_IREG

   !--------------------------------------------------------------------------------------

   subroutine idx_for2c_IEXT(nn,ivec)

   implicit none
   integer(IEXT), intent(in)    :: nn
   integer(IEXT), intent(inout) :: ivec(nn)
   integer(IEXT)                :: i

   do i = 1,nn
      ivec(i) = ivec(i) - 1
   enddo

   end subroutine idx_for2c_IEXT

   !--------------------------------------------------------------------------------------

end module C2FOR
