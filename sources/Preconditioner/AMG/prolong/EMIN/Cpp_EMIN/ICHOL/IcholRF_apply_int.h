interface IcholRF_apply_int
   subroutine IcholRF_apply(nequ,nterm,iU,jU,mat_U,D_inv,vec,pvec)

   use class_precision

   implicit none  

   ! Input variables
   integer, intent(in)            :: nequ,nterm
   integer, intent(in)            :: iU(nequ+1),jU(nterm)
   real(kind=double), intent(in)  :: D_inv(nequ),mat_U(nterm),vec(nequ)

   ! Output variables
   real(kind=double), intent(out) :: pvec(nequ)

   end subroutine IcholRF_apply
end interface IcholRF_apply_int
