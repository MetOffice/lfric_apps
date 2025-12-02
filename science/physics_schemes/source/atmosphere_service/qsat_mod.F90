! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_qsat

! Overarching module for calculating saturation vapour pressure.
! Calculates:
! Saturation Specific Humidity (Qsat): Vapour to Liquid/Ice.
! Saturation Specific Humidity (Qsat_Wat): Vapour to Liquid.

! All return a saturation mixing ratio given a temperature and pressure
! using saturation vapour pressures calculated using the Goff-Gratch
! formulae, adopted by the WMO as taken from Landolt-Bornstein, 1987
! Numerical Data and Functional Relationships in Science and Technolgy.
! Group V/vol 4B meteorology. Phyiscal and Chemical properties or air, P35

! Note regarding values stored in the lookup tables (es):
! _wat versions     : over water above and below 0 deg c.
! non-_wat versions : over water above 0 deg c
!                     over ice below 0 deg c

! Method:
! Uses lookup tables to find eSAT, calculates qSAT directly from that.

! Documentation: UMDP No.29

module qsat_mod

! Module-wide use statements
use um_types, only: real_64, real_32

implicit none

! Set everything as private, and only expose the interfaces we want to
private
public :: qsat, qsat_wat, qsat_mix, qsat_wat_mix

interface qsat
  module procedure qsat_real64_1D, qsat_real32_1D,                             &
                   qsat_real64_2D, qsat_real32_2D,                             &
                   qsat_real64_3D, qsat_real32_3D,                             &
                   qsat_real64_scalar, qsat_real32_scalar,                     &
                   qsat_real64_1D_idx
end interface

interface qsat_wat
  module procedure qsat_wat_real64_1D, qsat_wat_real32_1D,                     &
                   qsat_wat_real64_2D, qsat_wat_real32_2D,                     &
                   qsat_wat_real64_3D, qsat_wat_real32_3D,                     &
                   qsat_wat_real64_scalar, qsat_wat_real32_scalar
end interface

interface qsat_mix
  module procedure qsat_mix_real64_1D, qsat_mix_real32_1D,                     &
                   qsat_mix_real64_2D, qsat_mix_real32_2D,                     &
                   qsat_mix_real64_3D, qsat_mix_real32_3D,                     &
                   qsat_mix_real64_scalar, qsat_mix_real32_scalar
end interface

interface qsat_wat_mix
  module procedure qsat_wat_mix_real64_1D, qsat_wat_mix_real32_1D,             &
                   qsat_wat_mix_real64_2D, qsat_wat_mix_real32_2D,             &
                   qsat_wat_mix_real64_3D, qsat_wat_mix_real32_3D,             &
                   qsat_wat_mix_real64_scalar, qsat_wat_mix_real32_scalar
end interface

character(len=*), parameter, private :: ModuleName='QSAT_MOD'

contains

! Do all the 64-bit routines first, then all the 32-bit ones

! ------------------------------------------------------------------------------
! 1D subroutines.
! These are called by all the other dimensionalities.
! There are 4 routines that provide every permutation of with/without
! _wat and _mix. The _mix versions apply a slightly different algorithm,
! and _wat uses different data tables.

subroutine qsat_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real64_mod, only:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon, one_minus_epsilon
use conversions_mod,      only: zerodegc
implicit none
integer, parameter :: prec = real_64
#include "qsat_mod_qsat.h"
end subroutine qsat_real64_1D

subroutine qsat_wat_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real64_mod, only:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon, one_minus_epsilon
use conversions_mod,      only: zerodegc
implicit none
integer, parameter :: prec = real_64
#include "qsat_mod_qsat.h"
end subroutine qsat_wat_real64_1D

subroutine qsat_mix_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real64_mod, only:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon
use conversions_mod,      only: zerodegc
implicit none
integer, parameter :: prec = real_64
#include "qsat_mod_qsat_mix.h"
end subroutine qsat_mix_real64_1D

subroutine qsat_wat_mix_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real64_mod, only:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon
use conversions_mod,      only: zerodegc
implicit none
integer, parameter :: prec = real_64
#include "qsat_mod_qsat_mix.h"
return
end subroutine qsat_wat_mix_real64_1D

! ------------------------------------------------------------------------------
! 1D indirect indexing subroutines
! This routine expects that the indirect indices are unique, i.e. not
! the same point is calculate twice.

subroutine qsat_real64_1D_idx(qs, t, p, npnts, idx, ni)
! Remember to use the correct flavour of es
use qsat_data_real64_mod, only:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon, one_minus_epsilon
use conversions_mod,      only: zerodegc
implicit none
integer, parameter :: prec = real_64
#define   QSAT_WITH_IDX 1
#include "qsat_mod_qsat.h"
#undef    QSAT_WITH_IDX
end subroutine qsat_real64_1D_idx


! ------------------------------------------------------------------------------
! 2D subroutines

subroutine qsat_real64_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_real64_2D

subroutine qsat_wat_real64_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_wat_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_wat_real64_2D

subroutine qsat_mix_real64_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_mix_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_mix_real64_2D

subroutine qsat_wat_mix_real64_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_wat_mix_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_wat_mix_real64_2D

! ------------------------------------------------------------------------------
! 3D subroutines

subroutine qsat_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_real64_3D

subroutine qsat_wat_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_wat_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_wat_real64_3D

subroutine qsat_mix_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_mix_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_mix_real64_3D

subroutine qsat_wat_mix_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_64),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_64),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_wat_mix_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_wat_mix_real64_3D

! ------------------------------------------------------------------------------
! scalar subroutines

subroutine qsat_real64_scalar(qs, t, p)
implicit none
real (kind=real_64),    intent(in)  :: t, p
real (kind=real_64),    intent(out) :: qs
real (kind=real_64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_real64_scalar

subroutine qsat_wat_real64_scalar(qs, t, p)
implicit none
real (kind=real_64),    intent(in)  :: t, p
real (kind=real_64),    intent(out) :: qs
real (kind=real_64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_wat_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_wat_real64_scalar

subroutine qsat_mix_real64_scalar(qs, t, p)
implicit none
real (kind=real_64),    intent(in)  :: t, p
real (kind=real_64),    intent(out) :: qs
real (kind=real_64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_mix_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_mix_real64_scalar

subroutine qsat_wat_mix_real64_scalar(qs, t, p)
implicit none
real (kind=real_64),    intent(in)  :: t, p
real (kind=real_64),    intent(out) :: qs
real (kind=real_64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_wat_mix_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_wat_mix_real64_scalar

!Now for all the 32-bit routines

! ------------------------------------------------------------------------------
! 1D subroutines.
! These are called by all the other dimensionalities.
! There are 4 routines that provide every permutation of with/without
! _wat and _mix. The _mix versions apply a slightly different algorithm,
! and _wat uses different data tables.

subroutine qsat_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real32_mod, only:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon          => repsilon_32b,             &
                                one_minus_epsilon => one_minus_epsilon_32b
use conversions_mod,      only: zerodegc   => zerodegc_32
implicit none
integer, parameter :: prec = real_32
#include "qsat_mod_qsat.h"
end subroutine qsat_real32_1D

subroutine qsat_wat_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real32_mod, only:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon          => repsilon_32b,             &
                                one_minus_epsilon => one_minus_epsilon_32b
use conversions_mod,      only: zerodegc   => zerodegc_32
implicit none
integer, parameter :: prec = real_32
#include "qsat_mod_qsat.h"
end subroutine qsat_wat_real32_1D

subroutine qsat_mix_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real32_mod, only:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon          => repsilon_32b
use conversions_mod,      only: zerodegc   => zerodegc_32
implicit none
integer, parameter :: prec = real_32
#include "qsat_mod_qsat_mix.h"
end subroutine qsat_mix_real32_1D

subroutine qsat_wat_mix_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
use qsat_data_real32_mod, only:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right kind
use planet_constants_mod, only: repsilon          => repsilon_32b
use conversions_mod,      only: zerodegc   => zerodegc_32
implicit none
integer, parameter :: prec = real_32
#include "qsat_mod_qsat_mix.h"
return
end subroutine qsat_wat_mix_real32_1D

! ------------------------------------------------------------------------------
! 2D subroutines

subroutine qsat_real32_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_real32_2D

subroutine qsat_wat_real32_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_wat_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_wat_real32_2D

subroutine qsat_mix_real32_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_mix_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_mix_real32_2D

subroutine qsat_wat_mix_real32_2D(qs, t, p, npntsi, npntsj)
implicit none
integer,               intent(in)  :: npntsi, npntsj
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj), p(npntsi,npntsj)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj)
integer                            :: j
do j=1, npntsj
  call qsat_wat_mix_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
end do
return
end subroutine qsat_wat_mix_real32_2D

! ------------------------------------------------------------------------------
! 3D subroutines

subroutine qsat_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_real32_3D

subroutine qsat_wat_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_wat_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_wat_real32_3D

subroutine qsat_mix_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_mix_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_mix_real32_3D

subroutine qsat_wat_mix_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
implicit none
integer,               intent(in)  :: npntsi, npntsj, npntsk
real (kind=real_32),    intent(in)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
real (kind=real_32),    intent(out) :: qs(npntsi,npntsj,npntsk)
integer                            :: j, k
do k=1, npntsk
  do j=1, npntsj
    call qsat_wat_mix_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  end do
end do
return
end subroutine qsat_wat_mix_real32_3D

! ------------------------------------------------------------------------------
! scalar subroutines

subroutine qsat_real32_scalar(qs, t, p)
implicit none
real (kind=real_32),    intent(in)  :: t, p
real (kind=real_32),    intent(out) :: qs
real (kind=real_32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_real32_scalar

subroutine qsat_wat_real32_scalar(qs, t, p)
implicit none
real (kind=real_32),    intent(in)  :: t, p
real (kind=real_32),    intent(out) :: qs
real (kind=real_32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_wat_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_wat_real32_scalar

subroutine qsat_mix_real32_scalar(qs, t, p)
implicit none
real (kind=real_32),    intent(in)  :: t, p
real (kind=real_32),    intent(out) :: qs
real (kind=real_32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_mix_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_mix_real32_scalar

subroutine qsat_wat_mix_real32_scalar(qs, t, p)
implicit none
real (kind=real_32),    intent(in)  :: t, p
real (kind=real_32),    intent(out) :: qs
real (kind=real_32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
call qsat_wat_mix_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
return
end subroutine qsat_wat_mix_real32_scalar

end module qsat_mod
