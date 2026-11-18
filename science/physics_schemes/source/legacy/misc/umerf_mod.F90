! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -------------------- module umErf_mod -----------------------------
!
!   Purpose: Provides a portable wrapper function for the erf intrinsic
!
!   -------------------------------------------------------------------
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Misc

module umErf_mod

use, intrinsic :: iso_c_binding, only: c_double

implicit none

private
public :: umErf

#if defined(UMERF_USE_LIBC)

interface
  function c_erf(x) bind(c,name="erf")

  import :: c_double

  implicit none

  real(kind=c_double)             :: c_erf
  real(kind=c_double), intent(in) :: x

  end function
end interface

! Use the libc standard C99 version

interface umErf
  module procedure umErf_scalar, umErf_vector
end interface

#endif

!------------------------------------------------------------------------------!
contains
!------------------------------------------------------------------------------!

#if !defined(UMERF_USE_LIBC)

elemental function umErf(x)

implicit none

real             :: umErf
real, intent(in) :: x

! Use the non-F2k3 standard version - which is elemental
umErf = erf(x)

end function umErf

!------------------------------------------------------------------------------!

#else

function umErf_scalar(x)

implicit none

real, intent(in)            :: x

real                        :: umErf_scalar

umErf_scalar = c_erf(real(x,kind=c_double))

end function umErf_scalar

!------------------------------------------------------------------------------!

function umErf_vector(x)

implicit none

real, intent(in)            :: x(:)

real                        :: umErf_vector(1:size(x))

integer                     :: i

do i=1,size(x)
  umErf_vector(i)=c_erf(real(x(i),kind=c_double))
end do

end function umErf_vector

!------------------------------------------------------------------------------!

#endif

end module umErf_mod
