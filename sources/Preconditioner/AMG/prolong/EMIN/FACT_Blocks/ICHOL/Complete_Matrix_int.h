interface Complete_Matrix_int
   subroutine Complete_Matrix(myid,firstrow,nproc,nrows_L,nrows_C,nterm_L,nterm_C,           &
&                             iat_L,ja_L,iat_C,ja_C,WI1,WI2,coef_L,coef_C)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)  :: myid,firstrow,nproc
   integer(IREG), intent(in)  :: nrows_L,nrows_C
   integer(IEXT), intent(in)  :: nterm_L,nterm_C
   integer(IEXT), intent(in)  :: iat_L(nrows_L+1)
   integer(IREG), intent(in)  :: ja_L(nterm_L)
   real(double), intent(in)   :: coef_L(nterm_L)

   ! Output variables
   integer(IEXT), intent(out) :: iat_C(nrows_C+1)
   integer(IREG), intent(out) :: ja_C(nterm_C)
   integer(IEXT), intent(out) :: WI1(nrows_C,nproc+1),WI2(nproc)
   real(double), intent(out)  :: coef_C(nterm_C)

   end subroutine Complete_Matrix
end interface Complete_Matrix_int
