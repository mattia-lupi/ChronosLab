#define IREG_LONG 0
#define IEXT_LONG 0

module class_precision
!*****************************************************************************************
!
!  Class: precision
!
!  Coded by Carlo Janna
!  December 2018
!
!  Purpose: set data precision codes.
!
!*****************************************************************************************

use iso_c_binding, only: c_int, c_long, c_double

implicit none

! Precision
integer, parameter      :: single = 4
integer, parameter      :: double = 8
integer, parameter      :: CB_DBL = c_double
#if IREG_LONG
integer, parameter      :: IREG = 8
integer, parameter      :: CB_IREG = c_long
#else
integer, parameter      :: IREG = 4
integer, parameter      :: CB_IREG = c_int
#endif
#if IEXT_LONG
integer, parameter      :: IEXT = 8
integer, parameter      :: CB_IEXT = c_long
#else
integer, parameter      :: IEXT = 4
integer, parameter      :: CB_IEXT = c_int
#endif
integer, parameter      :: IOMP = 4

! Parameters
real(double), parameter :: ZERO = 0._double
real(double), parameter :: ONE  = 1._double
real(double), parameter :: HALF = 0.5_double
real(double), parameter :: MONE = -1._double
real(double), parameter :: MACH_EPS = epsilon(0._double)

end module class_precision
