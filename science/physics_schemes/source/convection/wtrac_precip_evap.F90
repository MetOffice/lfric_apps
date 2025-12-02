! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_precip_evap_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Updates water tracers for evaporation of precipitation in the convection
!   scheme (as calculated in devap for downdraughts and pevp_bcb for
!   environment or below cloud base).
!   This routine is called twice to update the water tracers for the phase
!   changes in devap.  This is because the phase change can be in both
!   directions and isotopes behave differently for rain evaporation and
!   condensation.  So the first call deals with evaporation points and the
!   second call deals with condensation points.
!
! Method:
!   Calculates the amount of water changing phase and then updates the water
!   tracers.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'WTRAC_PRECIP_EVAP_MOD'

contains

! Subroutine Interface:
subroutine wtrac_precip_evap (npnts, n_wtrac, delp, q_km1, rain, snow,         &
                              th_km1, exkm1, qs_km1,                           &
                              wtrac_conv_old, q_km1_wtrac, timestep,           &
                              wtrac_dd2, l_evap_call, wtrac_ev)

use planet_constants_mod,    only: g
use wtrac_conv_mod,          only: conv_dd2_wtrac_type, conv_ev_wtrac_type
use wtrac_conv_store_mod,    only: conv_old_wtrac_type,                        &
                                   wtrac_dealloc_conv_store2
use wtrac_all_phase_chg_mod, only: wtrac_all_phase_chg

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!
! Subroutine arguments
integer, intent(in) ::  npnts     ! Vector length
integer, intent(in) ::  n_wtrac   ! No. of water tracers

real(kind=real_umphys), intent(in) :: delp(npnts)
                                  ! If called after pevp_bcb:
                                  !  Change in pressure across layer k-1 (Pa)
                                  ! If called after devap:
                                  ! Downdraught mass flux in layer K-1 (Pa/s)

! This next field is currently unused, but will be used in the future
real(kind=real_umphys), intent(in) :: q_km1(npnts)
                                  ! Mixing ratio in layer k-1
                                  ! after phase change(kg/kg)
real(kind=real_umphys), intent(in) :: rain(npnts)
                                  ! Amount of rain after phase chg (kg/m**2/s)
real(kind=real_umphys), intent(in) :: snow(npnts)
                                  ! Amount of snow after phase chg (kg/m**2/s)

! The next 3 fields are currently unused, but will be used in the future
real(kind=real_umphys), intent(in) :: th_km1(npnts)
                                  ! Potential temperature in layer k-1
real(kind=real_umphys), intent(in) :: exkm1(npnts)
                                  ! exner ratio for layer K-1
real(kind=real_umphys), intent(in) :: qs_km1(npnts)
                                  ! Saturated specific humidity in layer k-1

real(kind=real_umphys), intent(in out) :: q_km1_wtrac(npnts,n_wtrac)
                                  ! Water tracer mixing ratio in layer k-1
                                  ! (kg/kg)

type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old
                                  ! Water values before phase change

! If called for environment calculation (i.e. after pevp_bcb),
! the following 2 fields are used:
real(kind=real_umphys), optional, intent(in) :: timestep

type(conv_ev_wtrac_type), optional, intent(in out) :: wtrac_ev(n_wtrac)
                                  ! Structure containing water tracer fields
                                  ! used for phase change and evap of precip
                                  ! in environment or below cloud base

! Or if called for downdraught calculation (i.e. after devap),
! the following field is used:
logical, optional, intent(in) :: l_evap_call
                                  ! T = Model impact of evaporation only
                                  ! F = Model impact of condensation only

type(conv_dd2_wtrac_type), optional, intent(in out) :: wtrac_dd2(n_wtrac)
                                  ! Structure containing 2nd compression
                                  ! water tracer fields used in
                                  ! downdraught calculations


! Local variables
integer :: i, i_wt       ! Loop counters

! Working arrays
real(kind=real_umphys) :: rain_wtrac(npnts,n_wtrac ) ! Water tracer rainfall
real(kind=real_umphys) :: snow_wtrac(npnts,n_wtrac)  ! Water tracer snowfall

real(kind=real_umphys) :: conv_unit(npnts) ! Unit conversion for rain/snow
real(kind=real_umphys) :: q_km1_wtrac_old(npnts,n_wtrac)
                                           ! Store input q_km1_wtrac
real(kind=real_umphys) :: q_km1_wtrac_r(npnts,n_wtrac)
                                           ! q_km1_wtrac before and after rain
                                           !  phase change
real(kind=real_umphys) :: q_km1_wtrac_s(npnts,n_wtrac)
                                           ! q_km1_wtrac before and after snow
                                           !  phase change
real(kind=real_umphys) :: rain_old(npnts)  ! Amount of rainfall prior to
                                           ! phase change (kg/kg)
real(kind=real_umphys) :: snow_old(npnts)  ! Amount of snowfall prior to
                                           ! phase change(kg/kg)
real(kind=real_umphys) :: rain_new(npnts)  ! Rainfall after phase chg (kg/kg)
real(kind=real_umphys) :: snow_new(npnts)  ! Snowfall after phase chg (kg/kg)
real(kind=real_umphys) :: qchange_r(npnts) ! Amount of rain changing phase
real(kind=real_umphys) :: qchange_s(npnts) ! Amount of snow changing phase
real(kind=real_umphys) :: qchange_wtrac_r(npnts,n_wtrac)
                                           ! Amount of water tracer
                                           !  rain changing phase
real(kind=real_umphys) :: qchange_wtrac_s(npnts,n_wtrac)
                                           ! Amount of water tracer
                                           !  snow changing phase
real(kind=real_umphys) :: q_km1_new_r(npnts)! q_km1 updated for rain evap only
real(kind=real_umphys) :: q_km1_new_s(npnts)! q_km1 updated for snow evap only

logical :: l_devap               ! = T if called after devap
logical :: l_evap                ! = T if evaporation (rain/snow -> q)

! The following is in preparation for isotope code
character(len=8) :: process_txt

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_PRECIP_EVAP'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Code is slightly different depending on whether the phase change has
! been calculated in devap or pevp_bcb.  Set up logical to control this.
if (present(wtrac_dd2)) then
  ! Downdraughts
  l_devap = .true.
  l_evap = l_evap_call        ! Is this call for evaporation?
  if (l_evap) then
    process_txt = 'ddpr_evp'
  else
    process_txt = 'ddpr_con'
  end if
else
  ! Environment
  l_devap = .false.
  l_evap = .true.             ! Always evaporation
  process_txt = 'evpr_evp'
end if

! Set up unit conversion for changing rain and snow from kg/m2/s to kg/kg
if (l_devap) then
  ! units of delp are Pa/s
  do i = 1, npnts
    conv_unit(i) = g / delp(i)
  end do
else
  ! units of delp are Pa
  do i = 1, npnts
    conv_unit(i) = g * timestep / delp(i)
  end do
end if

! Set working arrays
if (l_devap) then
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      rain_wtrac(i,i_wt) = wtrac_dd2(i_wt)%rain(i)
      snow_wtrac(i,i_wt) = wtrac_dd2(i_wt)%snow(i)
    end do
  end do
else
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      rain_wtrac(i,i_wt) = wtrac_ev(i_wt)%rain_env(i)
      snow_wtrac(i,i_wt) = wtrac_ev(i_wt)%snow_env(i)
    end do
  end do
end if

! Store original q_km1_wtrac values
do i_wt = 1, n_wtrac
  do i = 1, npnts
    q_km1_wtrac_old(i,i_wt) = q_km1_wtrac(i,i_wt)
    q_km1_wtrac_r(i,i_wt)   = q_km1_wtrac(i,i_wt)
    q_km1_wtrac_s(i,i_wt)   = q_km1_wtrac(i,i_wt)
  end do
end do

do i = 1, npnts
  ! Temporarily convert units of rain and snow from kg/m2/s to kg/kg
  rain_old(i) = wtrac_conv_old%rain(i) * conv_unit(i)  ! Convert to kg/kg
  snow_old(i) = wtrac_conv_old%snow(i) * conv_unit(i)
  rain_new(i) = rain(i) * conv_unit(i)
  snow_new(i) = snow(i) * conv_unit(i)

  do i_wt = 1, n_wtrac
    rain_wtrac(i,i_wt) = rain_wtrac(i,i_wt) * conv_unit(i)
    snow_wtrac(i,i_wt) = snow_wtrac(i,i_wt) * conv_unit(i)
  end do

  ! Calculate amount of water changing phase
  ! (Note, wtrac_all_phase_chg only updates points where qchange is non zero)

  if (l_evap) then
    ! Evaporation
    qchange_r(i) = max(rain_old(i) - rain_new(i), 0.0)
    qchange_s(i) = max(snow_old(i) - snow_new(i), 0.0)
  else
    ! Condensation
    ! (Change direction so flux is vap -> rain/snow)
    qchange_r(i) = - min(rain_old(i) - rain_new(i), 0.0)
    qchange_s(i) = - min(snow_old(i) - snow_new(i), 0.0)
  end if

  ! Create temporary q_km1 values updated for each phase change
  q_km1_new_r(i) = wtrac_conv_old%q_km1(i) + qchange_r(i)
  q_km1_new_s(i) = wtrac_conv_old%q_km1(i) + qchange_s(i)

end do

! -----------------------------------------------------------------
! Update water tracers for phase change
! -----------------------------------------------------------------

if (l_evap) then

  ! Snow sublimation (snow -> vap)

  call wtrac_all_phase_chg(npnts, n_wtrac,                                     &
                           snow_old, wtrac_conv_old%q_km1, qchange_s,          &
                           snow_new, q_km1_new_s,                              &
                           'sno', 'vap', 'one_way', snow_wtrac, q_km1_wtrac_s, &
                           qchange_wtrac = qchange_wtrac_s)

  ! Rain evaporation (rain -> vap)

  call wtrac_all_phase_chg(npnts, n_wtrac,                                     &
                           rain_old, wtrac_conv_old%q_km1, qchange_r,          &
                           rain_new, q_km1_new_r,                              &
                           'rai', 'vap', 'one_way', rain_wtrac, q_km1_wtrac_r, &
                           qchange_wtrac = qchange_wtrac_r)

  ! Update water tracers q_km1_wtrac - do snow and rain processes in parallel
  ! (Note, rain_wtrac and snow_wtrac are correctly updated in
  !  wtrac_all_phase_chg)
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      q_km1_wtrac(i,i_wt) = q_km1_wtrac_old(i,i_wt)                            &
                          + qchange_wtrac_s(i,i_wt) + qchange_wtrac_r(i,i_wt)
    end do
  end do

else

  ! Deposition directly to snow (vap -> snow)
  call wtrac_all_phase_chg(npnts, n_wtrac,                                     &
                           wtrac_conv_old%q_km1, snow_old, qchange_s,          &
                           q_km1_new_s, snow_new,                              &
                           'vap', 'ice', 'one_way', q_km1_wtrac_s, snow_wtrac, &
                           qchange_wtrac = qchange_wtrac_s)

  ! Condensation directly to rain (vap -> rain)
  call wtrac_all_phase_chg(npnts, n_wtrac,                                     &
                           wtrac_conv_old%q_km1, rain_old, qchange_r,          &
                           q_km1_new_r, rain_new,                              &
                           'vap', 'liq', 'one_way', q_km1_wtrac_r, rain_wtrac, &
                           qchange_wtrac = qchange_wtrac_r)

  ! Update water tracers q_km1_wtrac - do snow and rain processes in parallel
  ! (although this makes no difference as the phase change can only involve
  !  rain or snow at each point)
  ! (Note, rain_wtrac and snow_wtrac are correctly updated in
  !  wtrac_all_phase_chg)
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      q_km1_wtrac(i,i_wt) = q_km1_wtrac_old(i,i_wt)                            &
                          - qchange_wtrac_s(i,i_wt) - qchange_wtrac_r(i,i_wt)
    end do
  end do

end if

! Conversion of water tracer rain/snow back to kg/m2/s (from kg/kg)
do i_wt = 1, n_wtrac
  do i = 1, npnts
    rain_wtrac(i,i_wt) = rain_wtrac(i,i_wt) / conv_unit(i)
    snow_wtrac(i,i_wt) = snow_wtrac(i,i_wt) / conv_unit(i)
  end do
end do

if (l_evap) then
  ! If all precip has been evaporated, ensure that water tracer is consistent
  do i = 1, npnts
    if (rain(i) <= 0.0) then
      do i_wt = 1, n_wtrac
        rain_wtrac(i,i_wt) = 0.0
      end do
    end if
    if (snow(i) <= 0.0) then
      do i_wt = 1, n_wtrac
        snow_wtrac(i,i_wt) = 0.0
      end do
    end if
  end do
end if

! Update fields in water tracer stucture
if (l_devap) then
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      wtrac_dd2(i_wt)%rain(i)    = rain_wtrac(i,i_wt)
      wtrac_dd2(i_wt)%snow(i)    = snow_wtrac(i,i_wt)
    end do
  end do
else
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      wtrac_ev(i_wt)%rain_env(i)    = rain_wtrac(i,i_wt)
      wtrac_ev(i_wt)%snow_env(i)    = snow_wtrac(i,i_wt)
      wtrac_ev(i_wt)%dqbydt_km1(i)  = wtrac_ev(i_wt)%dqbydt_km1(i)             &
               + (q_km1_wtrac(i,i_wt) - q_km1_wtrac_old(i,i_wt)) / timestep
    end do
  end do
end if

! Deallocate temporary arrays (but don't do this on first call from ddraught)
if (l_evap) then
  call wtrac_dealloc_conv_store2(wtrac_conv_old)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_precip_evap

end module wtrac_precip_evap_mod
