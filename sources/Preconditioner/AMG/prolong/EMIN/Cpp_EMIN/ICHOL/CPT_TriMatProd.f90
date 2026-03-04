!_----------------------------------------------------------------------------------------
! FUNCTION: CPT_TriMatProd
!
!> @brief High level interface for the computation of the triple sparse matrix product:
!!        S = F1 * A * F2
!!
!> @author Carlo Janna
!
!> @param[in]      lfil: number of off-diagonal entries per row in the upper part of the
!!                       resulting matrix. Must be >= 0.
!> @param[in]      nthreads: number of threads for openMP
!> @param[in]      nrows_F1,iat_F1,ja_F1,coef_F1: left matrix of the product
!> @param[in]      nrows_A,iat_A,ja_A,coef_A: central matrix of the product
!> @param[in]      nrows_F2,ncols_F2,iat_F2,ja_F2,coef_F2: rigth matrix of the product
!> @param[out]     nrows_S,iat_S,ja_S,coef_S: resulting factor of the product
!> @param[out]     ierr: error code 0 ---> success
!!                                  1 ---> dimensional inconsistency between inputs
!!                                  2 ---> negative value of lfil
!!                                  3 ---> allocation/deallocation error
!
!> @version 1.0
!
!> @date November 2020
!
!-----------------------------------------------------------------------------------------

function CPT_TriMatProd(lfil,nthreads,nrows_F1,ncols_F1,c_iat_F1,c_ja_F1,c_coef_F1,       &
&                       nrows_A,ncols_A,c_iat_A,c_ja_A,c_coef_A,                          &
&                       nrows_F2,ncols_F2,c_iat_F2,c_ja_F2,c_coef_F2,                     &
&                       nrows_S,ncols_S,c_iat_S,c_ja_S,c_coef_S) result(ierr)             &
&                       bind ( C, name="CPT_TriMatProd" )

!-----------------------------------------------------------------------------------------

use omp_lib
use class_precision
use C2FOR
use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

implicit none

! Interfaces
include './Complete_Matrix_int.h'
include './SPRS_3MP_int.h'

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

! A and P handles
integer(IEXT)                          :: nterm_A,nterm_S
integer(IEXT), pointer                 :: iat_A(:),iat_S(:)
integer(IREG), pointer                 :: ja_A(:),ja_S(:)
real(double), pointer                  :: coef_A(:),coef_S(:)

! mat_F1 and mat_F2 handles
integer(IEXT)                          :: nterm_F1,nterm_F2
integer(IEXT), pointer                 :: iat_F1(:),iat_F2(:)
integer(IREG), pointer                 :: ja_F1(:),ja_F2(:)
real(double), pointer                  :: coef_F1(:),coef_F2(:)

! Local variables
integer(IOMP)                          :: myid
integer(IREG)                          :: bsize,resto,nrows_L,firstrow
integer(IREG)                          :: alinfo,ierr_L
integer(IREG)                          :: nzr_A_avg,lfil_mean
integer(IREG)                          :: ind1,ind2
integer(IEXT)                          :: ind3,ind4
integer(IEXT)                          :: nz_S,nterm_S_loc

! Local allocatable variables
integer(IREG), allocatable             :: nnz_A(:),JWN2(:)
integer(IEXT), allocatable             :: WN2(:),WNP(:)
integer(IEXT), allocatable             :: WIG1(:,:),WIG2(:)
integer(IEXT), allocatable             :: proc_iatS(:)
integer(IREG), pointer                 :: proc_jaS(:)
real(double), allocatable              :: WR2(:)
real(double), pointer                  :: proc_coefS(:)

! Parameters
logical, parameter                     :: PRE_DROP = .false.

!-----------------------------------------------------------------------------------------

! Initialize the error code
ierr = 0

! Set handles for A
call c_f_pointer(c_iat_A,iat_A,[nrows_A+1])
nterm_A =  iat_A(nrows_A+1)
call c_f_pointer(c_ja_A,ja_A,[nterm_A])
call c_f_pointer(c_coef_A,coef_A,[nterm_A])
call idx_c2for_IEXT(int(nrows_A+1,IEXT),iat_A)
call idx_c2for_IREG(nterm_A,ja_A)

! Set handles for mat_F1
call c_f_pointer(c_iat_F1,iat_F1,[nrows_F1+1])
nterm_F1 =  iat_F1(nrows_F1+1)
call c_f_pointer(c_ja_F1,ja_F1,[nterm_F1])
call c_f_pointer(c_coef_F1,coef_F1,[nterm_F1])
call idx_c2for_IEXT(int(nrows_F1+1,IEXT),iat_F1)
call idx_c2for_IREG(nterm_F1,ja_F1)

! Set handles for mat_F2
call c_f_pointer(c_iat_F2,iat_F2,[nrows_F2+1])
nterm_F2 =  iat_F2(nrows_F2+1)
call c_f_pointer(c_ja_F2,ja_F2,[nterm_F2])
call c_f_pointer(c_coef_F2,coef_F2,[nterm_F2])
call idx_c2for_IEXT(int(nrows_F2+1,IEXT),iat_F2)
call idx_c2for_IREG(nterm_F2,ja_F2)

! Check dimensional consistency between:
if ( (ncols_F1 /= nrows_A) .or. (ncols_A /= nrows_F2) ) ierr = 1

! Check fill-in parameter:
if (lfil < 0) ierr = 2

! Compute number of rows for S
nrows_S = nrows_F1

! Allocate shared work arrays for P transposition
allocate(WIG1(nrows_S,nthreads+1),WIG2(nthreads),stat=alinfo)
if (alinfo .ne. 0) ierr = 3

if (ierr /= 0) return

! Parallel computation of P and R
nterm_S = 0

nzr_A_avg = max( 1 , nterm_A / nrows_A)
lfil_mean = min( lfil+1 , 2*nzr_A_avg )

!-----------------------------------------------------------------------------------------
!$OMP parallel                                                                            &
!$OMP private(myid,bsize,resto,firstrow,nrows_L,ind1,ind2,ind3,ind4,nz_S,                 &
!$OMP         proc_iatS,proc_jaS,proc_coefS,alinfo,nnz_A,nterm_S_loc,WN2,JWN2,WNP,WR2,    &
!$OMP         ierr_L)                                                                     &
!$OMP num_threads(nthreads)
!-----------------------------------------------------------------------------------------

! Retrieve processor ID
myid = omp_get_thread_num() + 1

! Set the thread partition
bsize = nrows_F1/ nthreads;
resto = mod(nrows_F1,nthreads);
if (myid <= resto) then
   nrows_L = bsize+1
   firstrow = (myid-1)*nrows_L + 1
else
   nrows_L = bsize
   firstrow = (myid-1)*nrows_L + resto + 1;
endif

! Initialize local error code
ierr_L = 0

! Allocate local scratches
allocate(nnz_A(nrows_L),WN2(ncols_A),WNP(ncols_F2),JWN2(ncols_A),WR2(ncols_A),stat=alinfo)
if (alinfo .ne. 0) ierr_L = 3

ind1 = firstrow
ind2 = ind1 + nrows_L

ind3 = iat_A(ind1)
ind4 = iat_A(ind2)
nz_S = lfil_mean * nrows_L + (ind4-ind3)

allocate(proc_iatS(nrows_L+1),proc_jaS(nz_S),proc_coefS(nz_S),stat=alinfo)
if (alinfo .ne. 0) ierr_L = 3

if (ierr_L .eq. 0) then
   nnz_A(1:nrows_L) = 1
   proc_iatS(1) = 1

   call SPRS_3MP(nz_S,firstrow-1,lfil,PRE_DROP,nrows_L,nterm_F1,                          &
&                iat_F1(firstrow:firstrow+nrows_L),ja_F1,coef_F1,nrows_A,ncols_A,nterm_A, &
&                iat_A,ja_A,coef_A,nrows_F2,ncols_F2,nterm_F2,iat_F2,ja_F2,coef_F2,       &
&                nnz_A,nterm_S_loc,proc_iatS,proc_jaS,proc_coefS,WN2,WNP,JWN2,WR2,ierr_L)
   if (ierr_L /= 0) ierr_L = 3
endif

! Deallocate local scratches
deallocate(nnz_A,WN2,WNP,JWN2,WR2,stat=alinfo)
if (alinfo .ne. 0) ierr_L = 3

! Handle errors in computation
!$OMP atomic
ierr = max(ierr,ierr_L)

! Matrix S

! Compute the total number of term
!$OMP atomic
nterm_S = nterm_S + nterm_S_loc

!$OMP barrier
if (ierr .ne. 0) goto 10000

! Completion of S matrix
! Compute the number of term of the complete matrix
!$OMP single
nterm_S = 2 * nterm_S - nrows_S
allocate(iat_S(nrows_S+1),ja_S(nterm_S),coef_S(nterm_S),stat=ierr)
if (ierr .ne. 0) ierr = 3
!$OMP end single
if (ierr .ne. 0) goto 10000

call Complete_Matrix(myid,firstrow,nthreads,nrows_L,nrows_S,nterm_S_loc,nterm_S,          &
&                    proc_iatS,proc_jaS,iat_S,ja_S,WIG1,WIG2,proc_coefS,coef_S)

deallocate(proc_iatS,proc_jaS,proc_coefS,stat=ierr_L)
if (ierr_L .ne. 0) ierr_L = 3

! Handle errors in computation
!$OMP atomic
ierr = max(ierr,ierr_L)

! Exit line in case of errors
10000 continue

!-----------------------------------------------------------------------------------------
!$OMP end parallel
!-----------------------------------------------------------------------------------------

! Deallocate shared work arrays for S computation and completition
deallocate(WIG1,WIG2,stat=alinfo)
if (alinfo .ne. 0) ierr = 3

! Check for errors
if (ierr .ne. 0) return

! Modify indices of S to C-style
ncols_S =  nrows_S
call idx_for2c_IEXT(int(nrows_S+1,IEXT),iat_S)
call idx_for2c_IREG(nterm_S,ja_S)
c_iat_S = c_loc(iat_S)
c_ja_S = c_loc(ja_S)
c_coef_S = c_loc(coef_S)

! Return indices of mat_A, mat_F1 and mat_F2 to C-style
call idx_for2c_IEXT(int(nrows_A+1,IEXT),iat_A)
call idx_for2c_IREG(nterm_A,ja_A)
call idx_for2c_IEXT(int(nrows_F1+1,IEXT),iat_F1)
call idx_for2c_IREG(nterm_F1,ja_F1)
call idx_for2c_IEXT(int(nrows_F2+1,IEXT),iat_F2)
call idx_for2c_IREG(nterm_F2,ja_F2)

end function CPT_TriMatProd
