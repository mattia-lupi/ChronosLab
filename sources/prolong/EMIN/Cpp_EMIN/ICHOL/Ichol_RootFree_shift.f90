subroutine Ichol_RootFree_shift(nequ,nterm,lfil,droptol,alpha,small,iat,ja,mat_A,iwk_U,   &
&                               it_U,jcol_U,mat_U,D,U_irow_lst,U_nxt_lst,U_end_lst,ind_U, &
&                               JWI,ineg,ierr)
!-----------------------------------------------------------------------------------------
!
!  Subroutine: Ichol_RootFree_shift
!
!  Coded by Carlo Janna
!  February 2011
!
!  Purpose: Compute an Incomplete [U]^T[D][U] (Root-Free) Factorization by using the
!           symmetric Crout elimination scheme and a dual threshold dropping strategy.
!           If a negative entry arise on the main diagonal, at the end of the
!           factorization process it is converted to its absolute value
!
!-----------------------------------------------------------------------------------------
!
!  On entry:
!
!  nequ          : system dimensions
!  nterm         : # of non-zero entries of A
!  iat,ja,mat_A  : system matrix stored in SSR format
!  lfil          : number of entries retained for each row of A in addition to the
!                  original non-zeroes
!  droptol       : tolerance relative to the diagonal entry below which elements are
!                  disregarded
!  alpha         : shift value to be added to the diagonal
!  small         : value below which an entry is considered zero
!  iwk_U         : total space reserved for the upper part of the preconditioner
!                  NB: has to be at least nequ
!
!  On exit:
!
!  D             : inverse of the diagonal factor
!  it_U,jcol_U
!  mat_U         : upper factor stored in CSR format
!  ierr          : Error code
!                  0 ---> Successful run
!                  1 ---> Error, insufficient storage for U
!                  2 ---> warning, negative entries on D
!  ineg          : row index of the first diagonal negative entry encountered
!
!  Scratch Vectors:
!
!  U_irow_lst,
!  U_nxt_lst,
!  U_end_lst     : column structure of the U factor stored as a list
!  ind_U         : pointer to the first entry of each row of U that is
!                  still to be used in the elimination process
!  JWI           : non-zero indicator vector
!
!-----------------------------------------------------------------------
use class_precision

implicit none

! Interfaces
include 'RI_SortSplit_int.h'
include 'IRheapsort_int.h'
include 'updt_list_int.h'

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

! Local Variables
integer(IEXT)                  :: U_lenlst
integer(IREG)                  :: i,jjcol,ncut
integer(IEXT)                  :: j,icurr
integer(IEXT)                  :: ind,ind1
integer(IREG)                  :: ind_UT_col
integer(IEXT)                  :: pt_U,end_U,U_stop
integer(IREG)                  :: len_U0,len_U
real(kind=double)              :: pivot,abstol,scal_fac

!-----------------------------------------------------------------------

! Initialize U lists

U_lenlst = nequ
do i = 1,nequ
   U_irow_lst(i) = 0
   U_end_lst(i) = i
enddo

! Initialize JWI
do i = 1,nequ
   JWI(i) = 0
enddo

! Load the diagonal of A in D

do i = 1,nequ
   D(i) = mat_A(iat(i)) + alpha
enddo

! Start the elimination process

it_U(1) = 1

! Main loop

main_loop: do i = 1,nequ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !write(*,*) 'i: ',i
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   pt_U = it_U(i)
   end_U = pt_U
   ind_U(i) = pt_U

!  Load the i-th row of A in U
   do j = iat(i)+1,iat(i+1)-1
      ind1 = ja(j)
      jcol_U(end_U) = ind1
      JWI(ind1) = end_U
      mat_U(end_U) = mat_A(j)
      end_U = end_U + 1
   enddo
   len_U0 = end_U - pt_U

!  Use the U list to perform the elimination
   icurr = i
   ind_UT_col = U_irow_lst(icurr)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !write(*,*) 'ind_UT_col: ',ind_UT_col
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   do while (ind_UT_col .ne. 0)
!     Take the correct position in the L column
      ind = ind_U(ind_UT_col)

!     Row Elimination

!     Compute the pivotal term
      pivot = mat_U(ind)*D(ind_UT_col)
!     Take the correct position in the corresponding U row
      ind1 = ind_U(ind_UT_col)
      U_stop = it_U(ind_UT_col+1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'PRIMA = ind1,U_stop:',ind1,U_stop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      if ((ind1.lt.U_stop) .and. (jcol_U(ind1).eq.i)) then
!        Update the diagonal term
         D(i) = D(i) - pivot*mat_U(ind1)
         ind1 = ind1 + 1
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'ind1,U_stop:',ind1,U_stop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!     Update the upper part
      do j = ind1,U_stop-1
         jjcol = jcol_U(j)
         if (JWI(jjcol) .eq. 0) then
!           Fill-in term
            mat_U(end_U) = - pivot*mat_U(j)
            jcol_U(end_U) = jjcol
            JWI(jjcol) = end_U
            end_U = end_U + 1
         else
!           No fill-in term
            mat_U(JWI(jjcol)) = mat_U(JWI(jjcol)) - pivot*mat_U(j)
         endif
      enddo
!     Update the position in the UT column
      ind_U(ind_UT_col) = ind + 1

!     Consider the next UT column in the list
      icurr = U_nxt_lst(icurr)
      ind_UT_col = U_irow_lst(icurr)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'icurr: ',icurr
      !write(*,*) 'ind_UT_col: ',ind_UT_col
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   enddo

!  Check if the diagonal entry is negative
   if (D(i) .lt. small) then
      ineg = i
      ierr = 2
      return
   endif

!  Check the U storage
   if (end_U .gt. iwk_U) then
      ierr = 1
      return
   endif

!  Reset non-zero indicators

   do j = pt_U,end_U-1
      JWI(jcol_U(j)) = 0
   enddo

!  Apply dropping strategy

   abstol = droptol*D(i)**2
   scal_fac = 1.0_double / D(i)

   ind = pt_U
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !write(*,*) 'DROPPO ind,end_U',ind,end_U
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   do j = pt_U,end_U-1
      if (abs(mat_U(j)) .gt. abstol) then
         mat_U(ind) = mat_U(j)
         jcol_U(ind) = jcol_U(j)
         ind = ind + 1
      endif
   enddo
   end_U = ind

!  Set the number of non-zero terms to retain

   len_U = end_U - pt_U
   ncut = max(0,min(len_U0+lfil,len_U))
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   !write(*,*) 'len_U,ncut      ',len_U,ncut
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   if (len_U .gt. 0) then

!     Select the ncut largest entries
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'RI_SortSplit   ',pt_U
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call RI_SortSplit(len_U,ncut,mat_U(pt_U),jcol_U(pt_U))

!     Order them in increasing column index
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'IRheapsort   ',pt_U,ncut
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call IRheapsort(jcol_U(pt_U),mat_U(pt_U),ncut)

!     Scale them by the diagonal term
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'XXXXX   '
      !write(*,*) 'DSCAL   ',ncut,scal_fac,pt_U,mat_U(pt_U)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !call DSCAL(ncut,scal_fac,mat_U(pt_U),1)
      do j = pt_U,pt_U+ncut-1
         mat_U(j) = scal_fac*mat_U(j)
      enddo

!     Update the U list
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !write(*,*) 'updt_list      ',i,ncut,jcol_U(pt_U)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      call updt_list(i,ncut,jcol_U(pt_U),nequ,iwk_U,U_lenlst,U_irow_lst,                  &
&                    U_nxt_lst,U_end_lst)

   endif

!  Update the pointer of the U row
   it_U(i+1) = pt_U + ncut

enddo main_loop

! Invert the diagonal terms
do i = 1,nequ
   D(i) = 1._double / D(i)
enddo

! Successful run
ierr = 0

end subroutine Ichol_RootFree_shift
