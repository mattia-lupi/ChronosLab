interface CPT_TriMatProd_int
   function CPT_TriMatProd(lfil,nthreads,nrows_F1,ncols_F1,c_iat_F1,c_ja_F1,c_coef_F1,       &
&                          nrows_A,ncols_A,c_iat_A,c_ja_A,c_coef_A,                          &
&                          nrows_F2,ncols_F2,c_iat_F2,c_ja_F2,c_coef_F2,                     &
&                          nrows_S,ncols_S,ic_iat_S,c_ja_S,c_coef_S) result(ierr)            &
&                          bind ( C, name="CPT_TriMatProd" )

   use class_precision
   use C2FOR
   use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

   implicit none

   ! Input variables
   integer(CB_IREG), value, intent(in)    :: lfil
   integer(CB_IREG), value, intent(in)    :: nthreads
   integer(CB_IREG), value, intent(in)    :: nrows_F1,ncols_F1,nrows_F2,ncols_F2
   integer(CB_IREG), value, intent(in)    :: nrows_A,ncols_A
   type(C_PTR), value, intent(in)         :: c_iat_F1,c_ja_F1
   type(C_PTR), value, intent(in)         :: c_iat_A,c_ja_A
   type(C_PTR), value, intent(in)         :: c_iat_F2,c_ja_F2
   type(C_PTR), value, intent(in)         :: c_coef_F1,c_coef_A,c_coef_F2

   ! Output variables
   integer(CB_IREG), intent(out)          :: nrows_S,ncols_S
   type(C_PTR), intent(out)               :: c_iat_S,c_ja_S
   type(C_PTR), intent(out)               :: c_coef_S

   ! Result variables
   integer(CB_IREG)                       :: ierr

   end function CPT_TriMatProd
end interface CPT_TriMatProd_int
