interface RI_SortSplit_int
   subroutine RI_SortSplit(n,ncut,R_vec,I_vec)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)        :: n, ncut

   ! Input/Ouput variables
   real(kind=double), intent(inout) :: R_vec(n)
   integer(IREG), intent(inout)     :: I_vec(n)

   end subroutine RI_SortSplit
end interface RI_SortSplit_int
