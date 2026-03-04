interface count_rowterms_compl_int
   subroutine count_rowterms_compl(nequ,nterm,ja,WI)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)  :: nequ
   integer(IEXT), intent(in)  :: nterm
   integer(IREG), intent(in)  :: ja(nterm)

   ! Output variables
   integer(IEXT), intent(out) :: WI(nequ)

   end subroutine count_rowterms_compl
end interface count_rowterms_compl_int
