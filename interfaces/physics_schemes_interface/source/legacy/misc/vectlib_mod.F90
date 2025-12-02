! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

module vectlib_mod

! Description:
! This routine acts as an interface to vector versions
! of intrinsics functions on a platform.
!
! Supported libraries:
! IBM's VMASS (compile using VMASS def)
! Documentation: http://www-01.ibm.com/support/docview.wss?uid=swg27005473
!
! Intel's MKL (compile with MKL def)
! Documentation:
! https://software.intel.com/en-us/mkl-developer-reference-fortran
!
! Default compiles to equivalent do loop over array
!
! INTERFACEs for 32/64 bit versions created as required. Implemented so far:
! -powr_v
! -exp_v
! -oneover_v
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

implicit none
private

public ::                                                                      &
  exp_v_interface,                                                             &
  powr_v_interface,                                                            &
  oneover_v_interface,                                                         &
  exp_v,                                                                       &
  powr_v,                                                                      &
  oneover_v,                                                                   &
  asin_v,                                                                      &
  rtor_v,                                                                      &
  cubrt_v,                                                                     &
  log_v,                                                                       &
  sqrt_v,                                                                      &
  atan2_v

interface exp_v_interface
  module procedure exp_v, exp_v_32b
end interface

interface powr_v_interface
  module procedure powr_v, powr_v_32b
end interface

interface oneover_v_interface
  module procedure oneover_v, oneover_v_32b
end interface

contains

subroutine exp_v(n,x,y)
use um_types, only: real_64,integer_32
implicit none
! Sets y(i) to the exponential function of x(i), for i=1,..,n
integer :: n
real (kind=real_64) :: y(n), x(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_EXP_V)
  ! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)
l_n=n
call vexp (y, x, l_n)
#elif defined(MKL_EXP_V)
call vdexp(n, x, y)
#else
do i=1, n
  y(i) = exp(x(i))
end do
#endif
return
end subroutine exp_v

subroutine exp_v_32b(n,x,y)
use um_types, only: real_32,integer_32
implicit none
! Sets y(i) to the exponential function of x(i), for i=1,..,n
integer :: n
real (kind=real_32) :: y(n), x(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_EXP_V)
  ! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)
l_n=n
call vsexp (y, x, l_n)
#elif defined(MKL_EXP_V)
call vsexp(n, x, y)
#else
do i=1, n
  y(i) = exp(x(i))
end do
#endif
return
end subroutine exp_v_32b

!-----------------------------------------------------------

subroutine powr_v(n, x, power, z)
use um_types, only: real_64,integer_32
implicit none
! Sets z(i) to x(i) raised to the power power, for i=1,..,n
integer :: n, i
real (kind=real_64) :: z(n), x(n), y(n), power
integer (kind=integer_32) :: l_n

#if defined(MKL_POWR_V)
  ! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)
l_n=n
do i=1, n
  y(i)=power
end do
call vpow (z, x, y, l_n)
#elif defined(MKL_POWR_V)
call vdpowx(n, x, power, z)
#else
do i=1, n
  z(i) = x(i)**power
end do
#endif
return
end subroutine powr_v

subroutine powr_v_32b(n, x, power, z)
use um_types, only: real_32,integer_32
implicit none
! Sets z(i) to x(i) raised to the power power, for i=1,..,n
integer :: n, i
real (kind=real_32) :: z(n), x(n), y(n), power
integer (kind=integer_32) :: l_n

#if defined(MKL_POWR_V)
  ! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)
l_n=n
do i=1, n
  y(i)=power
end do
call vspow (z, x, y, l_n)
#elif defined(MKL_POWR_V)
call vspowx(n, x, power, z)
#else
!Allowing this loop to vectorise with the Cray compiler version before
!8.4.0 and 32bits, breaks PROC comparability, even at safe
#if defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif
do i=1, n
  z(i) = x(i)**power
end do
#endif
return
end subroutine powr_v_32b

!-----------------------------------------------------------

subroutine rtor_v(n, x, y, z)
use um_types, only: real_64, integer_32
implicit none

! Sets z(i) to x(i) raised to the power y(i), for i=1,..,n

integer :: n
real (kind=real_64) :: z(n), x(n), y(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_RTOR_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vpow (z, x, y, l_n)

#elif defined(MKL_RTOR_V)

call vdpow(n, x, y, z)

#else

do i=1, n
  z(i) = x(i)**y(i)
end do

#endif

return
end subroutine rtor_v

!-----------------------------------------------------------

subroutine sqrt_v(n, x, y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the square root of x(i), for i=1,..,n

integer :: n
real (kind=real_64) :: x(n), y(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_SQRT_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vsqrt (y, x, l_n)

#elif defined(MKL_SQRT_V)

call vdsqrt(n, x, y)

#else

do i=1, n
  y(i) = sqrt(x(i))
end do

#endif

return
end subroutine sqrt_v

!-----------------------------------------------------------

subroutine oneover_v(n, x, y)
use um_types, only: real_64, integer_32
implicit none
! Sets y(i) to the reciprocal of x(i), for i=1,..,n
integer :: n
real (kind=real_64) :: x(n), y(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(VMASS)
l_n=n
call vrec (y, x, l_n)
#else
do i=1, n
  y(i) = 1/x(i)
end do
#endif
return
end subroutine oneover_v

subroutine oneover_v_32b(n, x, y)
use um_types, only: real_32, integer_32
implicit none
! Sets y(i) to the reciprocal of x(i), for i=1,..,n
integer :: n
real (kind=real_32) :: x(n), y(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(VMASS)
l_n=n
call vsrec (y, x, l_n)
#else
do i=1, n
  y(i) = 1/x(i)
end do
#endif
return
end subroutine oneover_v_32b

!-----------------------------------------------------------

subroutine log_v (n, x, y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the natural logarithm of x(i), for i=1,..,n

integer :: n
real (kind=real_64) :: x(n), y(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_LOG_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vlog (y, x, l_n)

#elif defined(MKL_LOG_V)

call vdln( n, x, y )

#else

do i=1, n
  y(i) = log(x(i))
end do

#endif

return
end subroutine log_v

!-----------------------------------------------------------

subroutine sin_v(n,x,y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the sin function of x(i), for i=1,..,n

integer :: n
real (kind=real_64) :: y(n), x(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_SIN_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vsin (y, x, l_n)

#elif defined(MKL_SIN_V)

call vdsin(n, x, y)

#else

do i=1, n
  y(i) = sin(x(i))
end do

#endif

return
end subroutine sin_v

!-----------------------------------------------------------

subroutine cos_v(n,x,y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the cos function of x(i), for i=1,..,n

integer :: n
real (kind=real_64) :: y(n), x(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_COS_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vcos (y, x, l_n)

#elif defined(MKL_COS_V)

call vdcos(n, x, y)

#else

do i=1, n
  y(i) = cos(x(i))
end do

#endif

return
end subroutine cos_v

!-----------------------------------------------------------
subroutine acos_v(n,x,y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the cos function of x(i), for i=1,..,n

integer :: n
real (kind=real_64) :: y(n), x(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_ACOS_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif



#if defined(VMASS)

l_n=n
call vacos (y, x, l_n)

#elif defined(MKL_ACOS_V)

call vdacos(n, x, y)

#else

do i=1, n
  y(i) = acos(x(i))
end do

#endif

return
end subroutine acos_v

!-----------------------------------------------------------

subroutine asin_v(n,x,y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the asin function of x(i), for i=1,..,n

integer :: n
real (kind=real_64) :: y(n), x(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_ASIN_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vasin (y, x, l_n)

#elif defined(MKL_ASIN_V)

call vdasin(n, x, y)

#else

do i=1, n
  y(i) = asin(x(i))
end do

#endif

return
end subroutine asin_v

!-----------------------------------------------------------

subroutine atan2_v(n,a,b,y)
use um_types, only: real_64, integer_32
implicit none

! Sets y(i) to the atan2 function of a(i),b(i), for i=1,..,n

integer :: n
real (kind=real_64) :: y(n), a(n),b(n)
integer (kind=integer_32) :: l_n
integer :: i

#if defined(MKL_ATAN2_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif

#if defined(VMASS)

l_n=n
call vatan2 (y, a, b, l_n)

#elif defined(MKL_ATAN2_V)

call vdatan2(n, a, b, y)

#else

do i=1, n
  y(i) = atan2(a(i),b(i))
end do

#endif

return
end subroutine atan2_v

!-----------------------------------------------------------

subroutine cubrt_v(n, x, y)

use um_types, only: real_64
#if defined(VMASS)
use um_types, only: integer_32
#endif

implicit none

! Sets y(i) to the cube root of x(i), for i=1,..,n

integer, intent(in)  :: n
real (kind=real_64), intent(in) :: x(n)
real (kind=real_64), intent(out) :: y(n)

! local variables
#if defined(VMASS)
integer (kind=integer_32) :: l_n
#endif
integer :: i

#if defined(MKL_CUBRT_V)
! Interfaces for MKL
include 'mkl_vml.f90'
#endif


#if defined(VMASS)
!! use IBM's vector maths library if available
l_n=n
call vcbrt ( y, x, l_n)
#elif defined(MKL_CUBRT_V)
!! use INTEL's vector maths library if available
call vdcbrt(n, x, y)
#else
do i=1, n
  y(i) = x(i)**(1.0/3.0)
end do
#endif

return
end subroutine cubrt_v

!-----------------------------------------------------------

end module vectlib_mod
