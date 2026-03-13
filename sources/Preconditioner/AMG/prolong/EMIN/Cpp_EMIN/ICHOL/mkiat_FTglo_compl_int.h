interface mkiat_FTglo_compl_int
   subroutine mkiat_FTglo_compl(myid,nrows,nequ,nproc,firstrow,WI1,nnz,iat_FT)

   use class_precision

   implicit none

   ! Input variables
   integer(IREG), intent(in)    :: myid,nrows,nequ,nproc,firstrow
   integer(IEXT), intent(in)    :: nnz(nproc)

   ! Input/Output variables
   integer(IEXT), intent(inout) :: WI1(nequ,nproc+1),iat_FT(nrows+1)

   end subroutine mkiat_FTglo_compl
end interface mkiat_FTglo_compl_int
