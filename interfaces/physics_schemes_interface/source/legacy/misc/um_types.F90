! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

! Parameters for 32 and 64 bit kinds

module um_types

use constants_mod, only: r_bl

use, intrinsic :: iso_c_binding, only:                                         &
  c_int8_t, c_int16_t

#if defined(_OPENMP)
!$ use omp_lib, only: openmp_version
#endif

implicit none
! Precision and range for 64 bit real
integer, parameter :: prec64  = 15
integer, parameter :: range64 = 307

! Precision and range for 32 bit real
integer, parameter :: prec32  = 6
integer, parameter :: range32 = 37

! Range for integers
integer, parameter :: irange64=15
integer, parameter :: irange32=9

! Range for small logicals
integer, parameter :: lrange1=1

! Kind for 64 bit real
integer, parameter :: real_64  = selected_real_kind(prec64,range64)
! Kind for 32 bit real
integer, parameter :: real_32  = selected_real_kind(prec32,range32)
! Kind for 64 bit integer
integer, parameter :: integer_64 = selected_int_kind(irange64)
! Kind for 32 bit integer
integer, parameter :: integer_32 = selected_int_kind(irange32)
! Kind for 16 bit integer
integer, parameter :: integer_16 = c_int16_t
! Kind for 8 bit integer
integer, parameter :: integer_8 = c_int8_t

!Scheme-specific precisions

!Large scale precipitation scheme
#if defined(LSPREC_32B)
integer, parameter :: real_lsprec = real_32
#else
integer, parameter :: real_lsprec = real_64
#endif

! USSP Gravity Wave Scheme
#if defined(USSPPREC_32B)
integer, parameter :: real_usprec = real_32
#else
integer, parameter :: real_usprec = real_64
#endif

! Kind for real variables in UM Physics that don't have a switchable type
! defined above.
! NOTE : This kind is applied to all real variables in the UM physics schemes.
! Be aware that this may cause issues with stashwork arrays if the precision
! is altered from real64.
integer, parameter :: real_umphys  = real_64

! Explicit kind type for reals used in routines imported by the UM and LFRic
integer, parameter :: real_jlslsm = real_umphys


! Smallest non-zero real
real(kind=real_umphys), parameter :: real_eps = epsilon(0.0)
real(kind=r_bl), parameter :: rbl_eps = epsilon(0.0_r_bl)

! Kind for use with OpenMP functions (from omp_lib).
! Equal to kind(0) if no OpenMP, else is equal to kind(openmp_version)
#if defined(_OPENMP)
integer, parameter :: integer_omp = kind(openmp_version)
#else
integer, parameter :: integer_omp = kind(0)
#endif

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we have tested.
integer, parameter :: logical_64 = integer_64
integer, parameter :: logical_32 = integer_32
! Kind for small logicals
integer, parameter :: log_small = selected_int_kind(lrange1)


end module um_types

