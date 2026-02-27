subroutine SPRS_3MP(ntmax_P,shift,lfil_P,PRE_DROP,n1,nt_F1,iat_F1,ja_F1,coef_F1,n2,m2,    &
&                   nt_F2,iat_F2,ja_F2,coef_F2,n3,m3,nt_F3,iat_F3,ja_F3,coef_F3,          &
&                   nnz,nt_P,iat_P,ja_P,coef_P,WN2,WNP,JWN2,WR2,info)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: SPRS_3MP
!
!  Coded by Carlo Janna
!  March 2018
!
!  Purpose: Compute the sparse product between 3 CSR matrices [P] = [F1][F2][F3]. The
!           input matrices are given in CSR format with the consistency check of the
!           number of rows and columns left to the calling routine. The resulting matrix
!           is supposed symmetric, hence only the lower part is returned in SSR format.
!           A dropping strategy is applied during the computation in such a way that only
!           lfil entries of each row of P are retained in excess to the value specified by
!           the input array nnz.
!
!  Variables:
!
!  [F1]
!  n1        : # of rows of F1
!  nt_F1     : # of non-zeroes of F1
!  iat_F1    : topology vector of F1
!  ja_F1     : column indices of F1
!  coef_F1   : coefficients of F1
!
!  [F2]
!  n2        : # of rows of F2 (should equal the # of columns of F1)
!  nt_F2     : # of non-zeroes of F2
!  iat_F2    : topology vector of F2
!  ja_F2     : column indices of F2
!  coef_F2   : coefficients of F2
!
!  [F3]
!  n3        : # of rows of F3 (should equal the # of columns of F2)
!  m3        : # of columns of F3
!  nt_F3     : # of non-zeroes of F3
!  istart_F3 : start point of each line of F3
!  iend_F3   : start point of the first entry of each next line of F3
!  ja_F3     : column indices of F3
!  coef_F3   : coefficients of F3
!
!  shift     : difference between the local and absolute row numbering of P
!              (i.e. the absolute row index of P is given by "irow_P+shift")
!  nnz       : integer array of dimension n1 specifing the # of non-zeroes to be
!              retained for each row of P
!  lfil      : # of entries retained in excess to those specified by nnz
!  PRE_DROP  : .true.  ---> perform dropping on the partial product [F1][F2]
!              .false. ---> do not perform dropping on the partial product [F1][F2]
!
!  [P]
!  nt_P      : # of non-zeroes of P
!  ntmax_P   : maximum # of non-zeroes of P
!  iat_P     : topology vector of P (first element of iat_P is used as starting point)
!  ja_P      : column indices of P
!  coef_P    : coefficients of P
!  info      : == 0 ---> successful run
!              \= 0 ---> error reallocating P
!
!  NOTE: the working arrays WN2 and JWN2 can have a smaller dimension than n3 (i.e. the
!        # of columns of the product [F1][F2] equal to the # of rows of F3)
!-----------------------------------------------------------------------------------------

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

! Local variables
integer(IREG)                             :: i,ii,ll,i_abs,jjF1,ncut
integer(IEXT)                             :: j,k
integer(IEXT)                             :: ISTRT,ISTOP,JSTRT,JSTOP
integer(IEXT)                             :: kk,JSTRT_I64,JSTOP_I64
integer(IREG)                             :: ind1,len1,len_P
integer(IREG)                             :: ind_P,P_start,ipos
real(double)                              :: WRF1,fac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!integer omp_get_thread_num,myid
!myid = omp_get_thread_num() + 1
!write(2000+myid,*) '-----------------------------------------' 
!write(2000+myid,*) 'DENTRO',myid,n1,n2,m2,n3,m3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!-----------------------------------------------------------------------------------------

! Initialize the non-zero indicators work arrays WN2 and WNP
WN2(1:m2) = 0
WNP(1:m3) = 0

! Use the first pointer of iat_P as starting point
ind_P = iat_P(1)
info = 0

! Main cycle on the F1 (and P) rows
main_loop: do i = 1,n1

   ! Explore the i-th row of F1
   ISTRT = iat_F1(i)
   ISTOP = iat_F1(i+1)-1
   ind1 = 0
   F1_row_loop: do j = ISTRT,ISTOP

      jjF1 = ja_F1(j)
      WRF1 = coef_F1(j)

      ! Explore the jjF1 row of F2
      JSTRT_I64 = iat_F2(jjF1)
      JSTOP_I64 = iat_F2(jjF1+1)-1
      F2_row_loop: do kk = JSTRT_I64,JSTOP_I64

         ! Column index of this F2 entry
         ii = ja_F2(kk)

         if (WN2(ii) .eq. 0) then
            ! This is a new entry to add in the i-th row of [F1][F2]
            ind1 = ind1 + 1
            JWN2(ind1) = ii
            WR2(ind1) = coef_F2(kk)*WRF1
            WN2(ii) = int(ind1,IEXT)
         else
            ! Sum this entry to the corresponding entry already present in the i-th row
            ! of [F1][F2]
            WR2(WN2(ii)) = WR2(WN2(ii)) + coef_F2(kk)*WRF1
         end if

      end do F2_row_loop

   end do F1_row_loop

   ! Now WR2 and JWN2 store the i-th row of the product [F1][F2]

   if (PRE_DROP) then
      ! Apply a preliminary dropping to WR2-JWN2
      len1 = max( 0 , min(ind1,nnz(i)+lfil_P) )
      call RI_SortSplit(ind1,len1,WR2,JWN2)
   else
      ! Don't apply dropping
      len1 = ind1
   end if

   ! Matrix P = F1 * F2[1]

   ! Sparse reset of the non-zero indicator WN2
   do ll = 1,ind1
      WN2(JWN2(ll)) = 0
   end do

   ! Compute the i-th row of P
   P_start = ind_P

   ! Explore the i-th row of [F1][F2]

   ! Absolute number of the P row
   i_abs = i + shift

   F1F2_row_loop: do ll = 1,len1

      jjF1 = JWN2(ll)
      WRF1 = WR2(ll)

      ! Explore the jjF1 row of F3
      JSTRT = iat_F3(jjF1)
      JSTOP = iat_F3(jjF1+1)-1

      F3_row_loop: do k = JSTRT,JSTOP

         ii = ja_F3(k)

         ! Consider the lower part only
         if (ii .le. i_abs) then

            if (WNP(ii) .eq. 0) then
               ! This is a new entry to add in the i-th row of P
               ja_P(ind_P) = ii
               coef_P(ind_P) = coef_F3(k)*WRF1
               WNP(ii) = ind_P
               ind_P = ind_P + 1

               ! Check if the P storage is enough
               if (ind_P .gt. ntmax_P) then
                  ! Compute the new size estimating the final size with the n1/i ratio
                  fac = 1.1_double * real(n1,double) / real(i,double)
                  ntmax_P = ceiling(fac*real(ind_P,double))
                  info = resize_ireg(ntmax_P,ja_P)
                  if (info .ne. 0) return
                  info = resize_r(ntmax_P,coef_P)
                  if (info .ne. 0) return
               end if
            else
               ! Sum this entry to the corresponding entry already present in the i-th
               ! row of P
               coef_P(WNP(ii)) = coef_P(WNP(ii)) + coef_F3(k)*WRF1
            end if

         end if

      end do F3_row_loop

   end do F1F2_row_loop

   ! Move the diagonal term in first position of the current P row
   ipos = WNP(i_abs)
   if (ipos .eq. 0) then
      ! Add diagonal term if missing
      ja_P(ind_P) = i_abs
      coef_P(ind_P) = 0.0_double
      ind_P = ind_P + 1
   else
      ! Save diagonal term
      call swapi(ja_P(ipos),ja_P(ind_P-1))
      call swapr(coef_P(ipos),coef_P(ind_P-1))
   end if

   ! Set the number of non-zeroes to retain
   len_P = int(ind_P - P_start)
   ncut = max( 1 , min(len_P,nnz(i)+lfil_P) ) - 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !write(2000+myid,'(10i10)') i,n1,nnz(i),P_start,ntmax_P
   !call flush(2000+myid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Apply the dropping strategy
   call RI_SortSplit(len_P-1,ncut,coef_P(P_start:P_start+len_P-2),                        &
&                    ja_P(P_start:P_start+len_P-2))
   call IRheapsort(ja_P(P_start:P_start+ncut-1),coef_P(P_start:P_start+ncut-1),ncut)

   ! Sparse reset of the non-zero indicator WNP
   do kk = P_start,ind_P-1
      WNP(ja_P(kk)) = 0
   end do

   ! Restore diagonal term
   ja_P(P_start+ncut) = ja_P(ind_P-1)
   coef_P(P_start+ncut) = coef_P(ind_P-1)

   ! Assign the pointer of the next P row
   iat_P(i+1) = P_start + ncut + 1
   ind_P = iat_P(i+1)

end do main_loop

! Compute the # number of non-zeroes of P
nt_P = iat_P(n1+1)-iat_P(1)

!-----------------------------------------------------------------------------------------

end subroutine SPRS_3MP
