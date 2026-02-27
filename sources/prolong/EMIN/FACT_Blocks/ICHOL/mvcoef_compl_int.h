interface mvcoef_compl_int
   subroutine mvcoef_compl(firstrow,nrows,nequ,nterm_F,iat_F,ja_F,ja_FT,mat_F,mat_FT,punt)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)    :: firstrow,nrows,nequ
   integer(IEXT), intent(in)    :: nterm_F
   integer(IREG), intent(in)    :: ja_F(nterm_F)
   integer(IEXT), intent(in)    :: iat_F(nrows+1)
   real(double), intent(in)     :: mat_F(nterm_F)

   ! Input/Output variables
   integer(IEXT), intent(inout) :: punt(nequ)

   ! Output variables
   integer(IREG), intent(out)   :: ja_FT(nterm_F)
   real(double), intent(out)    :: mat_FT(nterm_F)

   end subroutine mvcoef_compl
end interface mvcoef_compl_int
