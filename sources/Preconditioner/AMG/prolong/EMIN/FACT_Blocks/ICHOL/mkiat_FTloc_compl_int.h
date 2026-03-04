interface mkiat_FTloc_compl_int
   subroutine mkiat_FTloc_compl(nrows,nequ,nproc,firstrow,WI,iat_FT,nnz)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)  :: nrows,nequ,nproc,firstrow
   integer(IEXT), intent(in)  :: WI(nequ,nproc+1)

   ! Output variables
   integer(IEXT), intent(out) :: iat_FT(nrows)
   integer(IEXT), intent(out) :: nnz

   end subroutine mkiat_FTloc_compl
end interface mkiat_FTloc_compl_int
