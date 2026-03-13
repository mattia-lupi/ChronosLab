interface SPRS_3MP_int
   subroutine SPRS_3MP(ntmax_P,shift,lfil_P,PRE_DROP,n1,nt_F1,iat_F1,ja_F1,coef_F1,n2,m2,    &
&                      nt_F2,iat_F2,ja_F2,coef_F2,n3,m3,nt_F3,iat_F3,ja_F3,coef_F3,          &
&                      nnz,nt_P,iat_P,ja_P,coef_P,WN2,WNP,JWN2,WR2,info)

   use class_precision
   use class_resize

   implicit none

   ! Input variables
   logical, intent(in)                       :: PRE_DROP
   integer(IREG), intent(in)                 :: shift,lfil_P
   integer(IREG), intent(in)                 :: n1,n2,m2,n3,m3
   integer(IEXT), intent(in)                 :: nt_F1,nt_F2,nt_F3
   integer(IEXT), intent(in)                 :: iat_F1(n1+1),iat_F2(n2+1),iat_F3(n3+1)
   integer(IREG), intent(in)                 :: ja_F1(nt_F1),ja_F2(nt_F2),ja_F3(nt_F3)
   integer(IREG), intent(in)                 :: nnz(n1)
   real(double), intent(in)                  :: coef_F1(nt_F1)
   real(double), intent(in)                  :: coef_F2(nt_F2)
   real(double), intent(in)                  :: coef_F3(nt_F3)

   ! Input/Output variables
   integer(IEXT), intent(inout)              :: ntmax_P
   integer(IEXT), intent(inout)              :: iat_P(n1+1)
   integer(IREG), pointer, intent(inout)     :: ja_P(:)
   real(double), pointer, intent(inout)      :: coef_P(:)

   ! Output variables
   integer(IREG), intent(out)                :: info
   integer(IEXT), intent(out)                :: nt_P

   ! Scratch arrays
   integer(IEXT), intent(out)                :: WN2(m2)
   integer(IEXT), intent(out)                :: WNP(m3)
   integer(IREG), intent(out)                :: JWN2(m2)
   real(double), intent(out)                 :: WR2(m2)

   end subroutine SPRS_3MP
end interface SPRS_3MP_int
