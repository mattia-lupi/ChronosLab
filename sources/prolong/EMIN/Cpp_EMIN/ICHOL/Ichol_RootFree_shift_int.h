interface Ichol_RootFree_shift_int
   subroutine Ichol_RootFree_shift(nequ,nterm,lfil,droptol,alpha,small,iat,ja,mat_A,iwk_U,   &
&                                  it_U,jcol_U,mat_U,D,U_irow_lst,U_nxt_lst,U_end_lst,ind_U, &
&                                  JWI,ineg,ierr)

   use class_precision

   implicit none

   ! Input Variables

   ! System Matrix
   integer(IREG), intent(in)      :: nequ
   integer(IEXT), intent(in)      :: nterm
   integer(IEXT), intent(in)      :: iat(nequ+1)
   integer(IREG), intent(in)      :: ja(nterm)
   real(kind=double), intent(in)  :: mat_A(nterm)
   ! Preconditioner Parameters
   integer(IREG), intent(in)      :: lfil
   real(kind=double), intent(in)  :: droptol,alpha,small
   ! Preconditioner Dimensions
   integer(IEXT), intent(in)      :: iwk_U

   ! Output Variables
   integer(IREG), intent(out)     :: ierr,ineg
   integer(IEXT), intent(out)     :: it_U(nequ+1)
   integer(IREG), intent(out)     :: jcol_U(iwk_U)
   real(kind=double), intent(out) :: D(nequ),mat_U(iwk_U)

   ! Scratch Variables (given as Input by the calling routine)
   integer(IREG), intent(out)     :: U_irow_lst(iwk_U+nequ)
   integer(IEXT), intent(out)     :: U_nxt_lst(iwk_U+nequ)
   integer(IEXT), intent(out)     :: U_end_lst(nequ)
   integer(IEXT), intent(out)     :: ind_U(nequ)
   integer(IEXT), intent(out)     :: JWI(nequ)

   end subroutine Ichol_RootFree_shift
end interface Ichol_RootFree_shift_int
