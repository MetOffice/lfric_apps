! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  UM module for calculating humidity-related fields for UKCA.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
! Code Description:
!   Language:  FORTRAN 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

module atmos_ukca_humidity_mod

use parkind1,              only: jpim, jprb      ! DrHook
use yomhook,               only: lhook, dr_hook  ! DrHook

implicit none
private

character(len=*), parameter :: ModuleName='ATMOS_UKCA_HUMIDITY_MOD'

public :: atmos_ukca_humidity

contains

!------------------------------------------------------------------------------
subroutine atmos_ukca_humidity(row_length, rows, model_levels,                 &
                               l_using_rh, l_using_rh_clr, l_using_svp,        &
                               t_theta_levels, p_theta_levels,                 &
                               q_mr, qcf_mr, cloud_liq_frac, cloud_frac,       &
                               rel_humid_frac, rel_humid_frac_clr, svp)
!------------------------------------------------------------------------------
! Calculate relative humidity, clear-sky relative humidity and/or satuartion
! vapour pressure as indicated by the logical arguments given.
! Allocatable arguments need only be allocated if they will be used.
! *****************************************************************************
!!!! NOTE: Corrections are required to calculations below which currently treat
!!!! specific humidity and water vapour mixing ratio as interchangeable.
!!!! The necessary corrections are detailed in inline comments. 
!------------------------------------------------------------------------------

use qsat_mod,              only: qsat, qsat_wat_mix
use planet_constants_mod,  only: repsilon

implicit none

! Subroutine arguments

integer, intent(in) :: row_length
integer, intent(in) :: rows
integer, intent(in) :: model_levels

! Switches to request individual fields
logical, intent(in) :: l_using_rh
logical, intent(in) :: l_using_rh_clr
logical, intent(in) :: l_using_svp

real, intent(in) :: t_theta_levels(row_length, rows, model_levels)
  ! Temperature on theta levels (K)
real, intent(in) :: p_theta_levels(row_length, rows, model_levels)
  ! Pressure on theta levels (Pa)
real, intent(in) :: q_mr(row_length, rows, model_levels)
  ! Water vapour mixing ratio 
real, allocatable, intent(in) :: qcf_mr(:,:,:)
  ! Cloud ice mixing ratio (for clear-sky RH only)
real, allocatable, intent(in) :: cloud_liq_frac(:,:,:)
  ! Liquid cloud fraction (for clear-sky RH only)
real, allocatable, intent(in) :: cloud_frac(:,:,:)
  ! Bulk cloud fraction (for clear-sky RH only)
real, allocatable, intent(in out) :: rel_humid_frac(:,:,:)
  ! Relative humidity (if requested)
real, allocatable, intent(in out) :: rel_humid_frac_clr(:,:,:)
  ! Clear sky relative humidity (if requested)
real, allocatable, intent(in out) :: svp(:,:,:)
  ! Saturation vapour pressure with respect to liquid water
  ! irrespective of temperature (if requested)

! Local variables

integer :: i 
integer :: j
integer :: k

real, allocatable :: qsatmr_wat(:,:,:)
  ! sat mixing ratio with respect to liquid water irrespective of temperature
  ! Saturation mixing ratio with respect to liquid water irrespective of
  ! temperature
real, allocatable :: qsatmr_ice_wat(:,:,:)
  ! Saturation mixing ratio with respect to liquid water (T > 0C) or
  ! ice (T < 0C)

integer (kind=jpim), parameter :: zhook_in  = 0  ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1  ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle   ! DrHook tracing

character(len=*), parameter :: RoutineName='ATMOS_UKCA_HUMIDITY'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

allocate(qsatmr_wat(row_length, rows, model_levels))
if (l_using_rh_clr) allocate(qsatmr_ice_wat(row_length, rows, model_levels))

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k)                    &
!$OMP SHARED(model_levels,qsatmr_wat,t_theta_levels,                           &
!$OMP p_theta_levels,l_using_rh_clr,qsatmr_ice_wat,l_using_rh,q_mr,            &
!$OMP rel_humid_frac,l_using_svp,svp,repsilon)
do k = 1, model_levels
  ! Calculate saturation mixing ratio with respect to water.
  call qsat_wat_mix(qsatmr_wat(:,:,k),                                         &
                    t_theta_levels(:,:,k), p_theta_levels(:,:,k),              &
                    size(t_theta_levels(:,1,k)), size(t_theta_levels(1,:,k)))
  ! Calculate saturation mixing ratio with respect to
  ! ice for T < 0C or water for T > 0C.
  !!!! ************************************************************************
  !!!! CORRECTION NEEDED: qsat returns saturation specific humidity but
  !!!! saturation mixing ratio is required so qsat_mix should be
  !!!! called instead.
  !!!! ************************************************************************
  if (l_using_rh_clr)                                                          &
    call qsat(qsatmr_ice_wat(:,:,k),                                           &
              t_theta_levels(:,:,k), p_theta_levels(:,:,k),                    &
              size(t_theta_levels(:,1,k)), size(t_theta_levels(1,:,k)))
  ! Calculate relative humidity fraction
  if (l_using_rh) rel_humid_frac(:,:,k) = q_mr(:,:,k) / qsatmr_wat(:,:,k)
  ! Derive saturation vapour pressure from saturated mixing ratio
  !!!! ************************************************************************
  !!!! CORRECTION NEEDED: the formula below is that for deriving SVP from
  !!!! specific humidity but qsatmr_wat is a mixing ratio.
  !!!! ************************************************************************
  if (l_using_svp)                                                             &
    svp(:,:,k) = qsatmr_wat(:,:,k) * p_theta_levels(:,:,k) /                   &
                 (repsilon + qsatmr_wat(:,:,k) * (1.0 - repsilon))
end do
!$OMP END PARALLEL DO

! Cap relative humidity values: 0.0 <= RH < 1.0
!!!! **************************************************************************
!!!! Note that values of 1.0 or greater will be capped to 0.999 but
!!!! values between 0.999 and 1.0 will not be modified. This lack of
!!!! consistency is undesirable. It is noted but left unchanged at
!!!! present to preserve bit-comparability of results.
!!!! **************************************************************************
if (l_using_rh) then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, j, k)              &
!$OMP SHARED(row_length, rows, model_levels, rel_humid_frac)
  do k = 1, model_levels
    do j = 1, rows
      do i = 1, row_length
        if (rel_humid_frac(i,j,k) < 0.0) then
          rel_humid_frac(i,j,k) = 0.0
        else if (rel_humid_frac(i,j,k) >= 1.0) then
          rel_humid_frac(i,j,k) = 0.999
        end if
      end do
    end do
  end do
!$OMP END PARALLEL DO
end if

if (l_using_rh_clr)                                                            &
  call calc_clear_sky_rh(row_length, rows, model_levels,                       &
                         q_mr, qsatmr_ice_wat, qsatmr_wat, qcf_mr,             &
                         cloud_liq_frac, cloud_frac, rel_humid_frac_clr)

if (allocated(qsatmr_ice_wat)) deallocate(qsatmr_ice_wat)
if (allocated(qsatmr_wat)) deallocate(qsatmr_wat)

! Cap clear sky relative humidity values: 0.0 <= RH < 1.0
!!!! **************************************************************************
!!!! Note that values of 1.0 or greater will be capped to 0.999 but
!!!! values between 0.999 and 1.0 will not be modified. This lack of
!!!! consistency is undesirable. It is noted but left unchanged at
!!!! present to preserve bit-comparability of results.
!!!! **************************************************************************
if (l_using_rh_clr) then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, j, k)              &
!$OMP SHARED(row_length, rows, model_levels, rel_humid_frac_clr)
  do k = 1, model_levels
    do j = 1, rows
      do i = 1, row_length
        if (rel_humid_frac_clr(i,j,k) < 0.0) then
          rel_humid_frac_clr(i,j,k) = 0.0
        else if (rel_humid_frac_clr(i,j,k) >= 1.0) then
          rel_humid_frac_clr(i,j,k) = 0.999
        end if
      end do
    end do
  end do
!$OMP END PARALLEL DO
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine atmos_ukca_humidity

!------------------------------------------------------------------------------
subroutine calc_clear_sky_rh(row_length, rows, model_levels,                   &
                             q_mr, qsatmr_ice_wat, qsatmr_wat, qcf_mr,         &
                             cloud_liq_frac, cloud_frac, rel_humid_frac_clr)
!------------------------------------------------------------------------------

use lsp_subgrid_mod, only: lsp_qclear
use cloud_inputs_mod, only: rhcrit

implicit none

! Subroutine arguments

integer, intent(in) :: row_length
integer, intent(in) :: rows
integer, intent(in) :: model_levels

real, intent(in) :: q_mr(:,:,:)
real, intent(in) :: qsatmr_ice_wat(:,:,:)
real, intent(in) :: qsatmr_wat(:,:,:)
real, intent(in) :: qcf_mr(:,:,:)
real, intent(in) :: cloud_liq_frac(:,:,:)
real, intent(in) :: cloud_frac(:,:,:)

real, intent(in out) :: rel_humid_frac_clr(:,:,:)

! Local variables

integer :: k

real, allocatable :: q_mr_1d(:)
real, allocatable :: qsatmr_ice_wat_1d(:)
real, allocatable :: qsatmr_wat_1d(:)
real, allocatable :: qcf_mr_1d(:)
real, allocatable :: cloud_liq_frac_1d(:)
real, allocatable :: cloud_frac_1d(:)
real, allocatable :: rhcrit_1d(:)
real, allocatable :: q_mr_clr_1d(:)
real, allocatable :: rel_humid_frac_clr_1d(:)

integer (kind=jpim), parameter :: zhook_in  = 0  ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1  ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle   ! DrHook tracing

character(len=*), parameter :: RoutineName='CALC_CLEAR_SKY_RH'

! End of header

if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!!!! CALC SIZE ONCE AND PASS IN
!!!! USE AUTOMATIC STORAGE?

  allocate(q_mr_1d(row_length*rows))
  allocate(qsatmr_ice_wat_1d(row_length*rows))
  allocate(qsatmr_wat_1d(row_length*rows))
  allocate(qcf_mr_1d(row_length*rows))
  allocate(cloud_liq_frac_1d(row_length*rows))
  allocate(cloud_frac_1d(row_length*rows))
  allocate(rhcrit_1d(row_length*rows))
  allocate(q_mr_clr_1d(row_length*rows))
  allocate(rel_humid_frac_clr_1d(row_length*rows))

  do k = 1, model_levels

    q_mr_1d           = reshape(q_mr(:,:,k), [row_length * rows])
    qsatmr_ice_wat_1d = reshape(qsatmr_ice_wat(:,:,k), [row_length * rows])
    qsatmr_wat_1d     = reshape(qsatmr_wat(:,:,k), [row_length * rows])
    qcf_mr_1d         = reshape(qcf_mr(:,:,k), [row_length * rows])
    cloud_liq_frac_1d = reshape(cloud_liq_frac(:,:,k), [row_length * rows])
    cloud_frac_1d     = reshape(cloud_frac(:,:,k), [row_length * rows])

    rhcrit_1d(:)      = rhcrit(k)

    call lsp_qclear(q_mr_1d, qsatmr_ice_wat_1d, qsatmr_wat_1d,                 &
         qcf_mr_1d, cloud_liq_frac_1d, cloud_frac_1d,                          &
         rhcrit_1d, q_mr_clr_1d, size(q_mr_1d))

    rel_humid_frac_clr_1d(:) = q_mr_clr_1d(:) / qsatmr_wat_1d(:)

    rel_humid_frac_clr(:,:,k) = reshape(rel_humid_frac_clr_1d,                  &
                                        [row_length, rows])

  end do

  deallocate(rel_humid_frac_clr_1d)
  deallocate(q_mr_clr_1d)
  deallocate(rhcrit_1d)
  deallocate(cloud_frac_1d)
  deallocate(cloud_liq_frac_1d)
  deallocate(qcf_mr_1d)
  deallocate(qsatmr_wat_1d)
  deallocate(qsatmr_ice_wat_1d)
  deallocate(q_mr_1d)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_clear_sky_rh

end module atmos_ukca_humidity_mod

