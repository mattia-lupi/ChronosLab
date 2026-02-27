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

!  Local variables
integer(IREG)                :: i,ind
integer(IEXT)                :: jnd

do i = 1,n
   lenlist = lenlist + 1
   ind = indices(i)
   jnd = end_lst(ind)
   ind_lst(jnd) = ival
   nxt_lst(jnd) = lenlist
   ind_lst(lenlist) = 0
   end_lst(ind) = lenlist
enddo

end subroutine updt_list
