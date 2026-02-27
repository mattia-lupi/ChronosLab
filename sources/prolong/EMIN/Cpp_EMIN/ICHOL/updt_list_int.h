interface updt_list_int
   subroutine updt_list(ival,n,indices,nequ,iwk,lenlist,ind_lst,nxt_lst,end_lst)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)    :: nequ
   integer(IEXT), intent(in)    :: iwk
   integer(IREG), intent(in)    :: ival,n,indices(n)

   ! Input/Output variables
   integer(IEXT), intent(inout) :: lenlist,end_lst(nequ)

   ! Output variables
   integer(IREG), intent(out)   :: ind_lst(iwk+nequ)
   integer(IEXT), intent(out)   :: nxt_lst(iwk+nequ)

   end subroutine updt_list
end interface updt_list_int
