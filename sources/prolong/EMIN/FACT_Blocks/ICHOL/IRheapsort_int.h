interface IRheapsort_int
   subroutine IRheapsort(x1,x2,n)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)        :: n 

   ! Input/Output variables
   integer(IREG), intent(inout)     :: x1(n)
   real(kind=double), intent(inout) :: x2(n)

   end subroutine IRheapsort
end interface IRheapsort_int
