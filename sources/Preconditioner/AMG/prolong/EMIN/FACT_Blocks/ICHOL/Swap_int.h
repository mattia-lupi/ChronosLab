interface swapi_int
   subroutine swapi(i1,i2)

   implicit none

   integer, intent(inout) :: i1,i2

   end subroutine swapi
end interface swapi_int

interface swapr_int
   subroutine swapr(r1,r2)

   use class_precision

   implicit none

   real(kind=double), intent(inout) :: r1,r2

   end subroutine swapr
end interface swapr_int
