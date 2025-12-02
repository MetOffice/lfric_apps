! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module llcs

! Description:
!  Module containing all subroutines for the Lambert-Lewis
!  Convection Scheme (LLCS)
!  Further description of the scheme can be found below

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

use cv_run_mod,           only: llcs_cloud_precip,                             &
    llcs_opt_all_rain, llcs_opt_all_cloud,                                     &
    llcs_opt_crit_condens, llcs_opt_const_frac,                                &
    llcs_detrain_coef, llcs_rhcrit, llcs_timescale,                            &
    fac_qsat, mparwtr, qlmin,                                                  &
    llcs_rain_frac

use gen_phys_inputs_mod,  only: l_mr_physics

use nlsizes_namelist_mod, only: model_levels

use planet_constants_mod, only: cp, r, kappa, g, repsilon, p_zero

use timestep_mod,         only: timestep

use water_constants_mod,  only: lc, hcapv

use um_types, only: real_umphys

implicit none
private

public llcs_control

character(len=*), parameter :: ModuleName = 'LLCS'

contains

! ----------------------------------------------------------------------
! -- UM WRAPPER -- !
! ----------------------------------------------------------------------

subroutine llcs_control(                                                       &
    ! Inputs
    row_length, rows, t_s_arr, q_s_arr, qcl_s_arr, p_half_arr, p_full_arr,     &
    ! Outputs
    t_f_arr, q_f_arr, qcl_inc_arr, cfl_arr, rain_arr)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Inputs
integer, intent(in) :: row_length
integer, intent(in) :: rows

real(kind=real_umphys), intent(in) :: t_s_arr   (row_length, rows, model_levels)
real(kind=real_umphys), intent(in) :: q_s_arr   (row_length, rows, model_levels)
real(kind=real_umphys), intent(in) :: qcl_s_arr (row_length, rows, model_levels)
real(kind=real_umphys), intent(in) ::                                          &
                    p_half_arr(row_length, rows, model_levels + 1)
real(kind=real_umphys), intent(in) :: p_full_arr(row_length, rows, model_levels)

! Outputs
real(kind=real_umphys), intent(out) ::                                         &
                     t_f_arr    (row_length, rows, model_levels)
real(kind=real_umphys), intent(out) ::                                         &
                     q_f_arr    (row_length, rows, model_levels)
real(kind=real_umphys), intent(out) ::                                         &
                     qcl_inc_arr(row_length, rows, model_levels)
real(kind=real_umphys), intent(out) ::                                         &
                     cfl_arr    (row_length, rows, model_levels)
real(kind=real_umphys), intent(out) :: rain_arr   (row_length, rows)

! Locals
integer :: i
integer :: j
integer :: nlevels

integer(kind=jpim), parameter :: zhook_in = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName = 'LLCS_CONTROL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nlevels = model_levels

do j = 1, rows
  do i = 1, row_length
    call convection_scheme(                                                    &
        ! Inputs
        nlevels, t_s_arr(i, j, :), q_s_arr(i, j, :), qcl_s_arr(i, j, :),       &
        p_half_arr(i, j, :), p_full_arr(i, j, :),                              &
        ! Outputs
        t_f_arr(i, j, :), q_f_arr(i, j, :), qcl_inc_arr(i, j, :),              &
        cfl_arr(i, j, :), rain_arr(i, j))
  end do
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine llcs_control

! ----------------------------------------------------------------------
! -- THE MAIN LLCS subroutine -- !
! ----------------------------------------------------------------------
subroutine convection_scheme(                                                  &
    ! Inputs
    nlevels, theta_start, q_start, qcl_start, p_half, p_full,                  &
    ! Outputs
    theta_final, q_final, qcl_inc, cloud_frac, rain_sec)
!------------------------------------------------------------------------------
! Run Lambert-Lewis convection scheme.
!
! Required inputs are:
! potential temperature, specific humidity, pressure
! on full model levels, pressure on half model levels
! Potential temperature and specific humidity have (nlevels) levels.
! The arrays are set up such that, for example, theta_start(1) is
! potential temperature at the lowest model level and theta_start(nlevels) is
! potential temperature at the top of atmosphere.
!
! The convection scheme will output final potential temperature and
! final specific humidity, along with rainfall / second (kgm^-2s^-1)
! within the column.
!------------------------------------------------------------------------------
implicit none

! Inputs
integer, intent(in) :: nlevels ! number of full model levels
real(kind=real_umphys), intent(in) :: theta_start(:)
                                   ! initial pot. temp. at full model levels.
real(kind=real_umphys), intent(in) :: q_start(:) ! initial specific humidity
real(kind=real_umphys), intent(in) :: qcl_start(:)
                                 ! initial cloud liquid content
real(kind=real_umphys), intent(in) :: p_full(:) ! pressure at full model levels
real(kind=real_umphys), intent(in) :: p_half(:) ! pressure at half model levels
! Outputs
real(kind=real_umphys), intent(out) :: theta_final(:) ! final pot. temp. [K]
real(kind=real_umphys), intent(out) :: q_final(:)
                                ! final specific humidity [kg kg-1]
real(kind=real_umphys), intent(out) :: qcl_inc(:)
                                ! increment of cloud liquid water [kg kg-1]
real(kind=real_umphys), intent(out) :: cloud_frac(:) ! cloud fraction [1]
real(kind=real_umphys), intent(out) :: rain_sec
                              ! convective precipitation [kg m-2 s-1]
! Local variables
integer :: j_lev ! level loop index
integer :: convbase ! the lowest level from which any form of mixing occurs.
integer :: cloudbase ! level from which adjustment to a moist adiabat can occur
integer :: final_level ! last level of a convective event.
                       ! Set to -1 whilst the event is in progress
real(kind=real_umphys) :: theta_adjust     (1:nlevels)
                                     ! theta after correction to reference
                                     ! profile
real(kind=real_umphys) :: theta_noq_adjust (1:nlevels)
                                     ! dummy theta profile where only
                                     ! dry corr. happens, used for enthalpy
                                     ! conservation.
real(kind=real_umphys) :: theta_adj_rlx    (1:nlevels)
                                     ! theta_adjust corrected to account for
                                     ! timestep over timescale
real(kind=real_umphys) :: theta_noq_adj_rlx(1:nlevels)
                                     ! dry correction theta adjusted for
                                     ! timestep over timescale
real(kind=real_umphys) :: q_adjust         (1:nlevels)
                                     ! spec. hum. after correction to q_sat
real(kind=real_umphys) :: q_adj_rlx        (1:nlevels)
                                     ! q_adjust adjusted for timestep/timescale
real(kind=real_umphys) :: q_conserve       (1:nlevels)
                                     ! conserved sp. hum. in the column, other
                                     ! than moisture lost through rainfall
real(kind=real_umphys) :: q_dummy          (1:nlevels)
                                     ! dummy q profile to ensure moisture is
                                     ! not mixed twice
real(kind=real_umphys) :: delta_p          (1:nlevels)
                                     ! Change in pressure between the top and
                                     ! the bottom of a layer
real(kind=real_umphys) :: xmin             (1:nlevels)
                                     ! critical cloud condensate content
real(kind=real_umphys) :: rh_use
               ! RH adjustment limit, currently is constant and equal to 1.0

logical :: convection_flag ! indicator that convection takes place
logical :: moist_has_happened ! flag indicating whether moist convection occurs
logical :: moist_trigger ! flag for triggering moist convection

delta_p = p_half(1:nlevels) - p_half(2:nlevels+1)

! Initialise variables
call initialise(                                                               &
    ! Inputs
    nlevels, theta_start, q_start,                                             &
    ! Outputs
    final_level,                                                               &
    theta_adjust, theta_noq_adjust,                                            &
    theta_adj_rlx, theta_noq_adj_rlx, theta_final,                             &
    q_adjust, q_adj_rlx, q_conserve, q_dummy, q_final, qcl_inc,                &
    cloud_frac, xmin, rain_sec,                                                &
    moist_has_happened, convection_flag)

! Calculate critical cloud condensate
if (llcs_cloud_precip == llcs_opt_crit_condens) then
  call critical_condensate(                                                    &
      ! Inputs
      nlevels, theta_start,p_full,                                             &
      ! Outputs
      xmin)
end if

! Diagnose and then perform convection.
do j_lev = 1, nlevels-1
  if (final_level == j_lev) then
    call diagnose(                                                             &
        ! Inputs
        nlevels, j_lev, theta_adjust, q_adjust, p_full,                        &
        ! In/out
        convection_flag,                                                       &
        ! Outputs
        convbase, cloudbase, final_level, moist_trigger, rh_use)
    if (convbase /=  -1) then
      ! Convection: dry or moist
      if (cloudbase == -1) then
        call dry_convection(                                                   &
            ! Inputs
            nlevels, convbase,                                                 &
            ! Inputs/outputs
            theta_adjust, theta_noq_adjust, q_adjust,                          &
            ! Outputs
            final_level)
      else
        call moist_convection(                                                 &
            ! Inputs
            nlevels, convbase, cloudbase, rh_use, p_full,                      &
            ! Inputs/outputs
            theta_adjust, theta_noq_adjust, q_adjust,                          &
            moist_has_happened, moist_trigger,                                 &
            ! Outputs
            final_level)
      end if
    end if
  end if
end do

if (convection_flag) then
  ! Convection has happened so we need to perform the following steps.

  ! Timescale relaxation
  ! In case of moist convection this step is performed twice:
  !   - for the actual theta profile,
  !   - for the 'dry convection only' profile used for enthalpy conservation
  call time_relax(                                                             &
      ! Inputs
      theta_start, theta_adjust, q_start, q_adjust,                            &
      ! Inputs/outputs
      theta_adj_rlx, q_adj_rlx)
  if (moist_has_happened) then
    ! Second call for enthalpy conservation calculation
    call time_relax(                                                           &
        ! Inputs
        theta_start, theta_noq_adjust, q_start, q_adjust,                      &
        ! Inputs/outputs
        theta_noq_adj_rlx, q_dummy)
  end if

  ! Conserve Moisture
  call conserve_moisture(                                                      &
      ! Inputs
      nlevels, delta_p, q_start, q_adj_rlx,                                    &
      ! Outputs
      q_conserve)

  if (moist_has_happened) then
    ! Latent heating & equiv. amount of instantaneous rainfall
    call rain(                                                                 &
        ! Inputs
        nlevels, theta_start, theta_noq_adj_rlx, q_conserve, qcl_start, xmin,  &
        p_full, delta_p,                                                       &
        ! Inputs/outputs
        theta_adj_rlx,                                                         &
        ! Outputs
        rain_sec, q_final, qcl_inc, cloud_frac)
  end if

  ! Conserve enthalpy
  if (moist_has_happened) then
    call conserve_enthalpy(                                                    &
        ! Inputs
        nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx, p_full,        &
        delta_p,                                                               &
        ! Outputs
        theta_final)
  else
    call conserve_enthalpy(                                                    &
        ! Inputs
        nlevels, theta_start, theta_adj_rlx, theta_adj_rlx, p_full, delta_p,   &
        ! Outputs
        theta_final)
  end if

end if

! Note that if no convection happens, then theta_final, q_final and
! rainsec are left as they are initialised in 'initialise', so in this
! case theta_final = theta_start, q_final = q_start and rainsec = 0.
end subroutine convection_scheme
! ----------------------------------------------------------------------
! -- end OF MAIN LLCS subroutine -- !
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! -- ROUTINES CALLED BY THE CONVECTION SCHEME -- !
! ----------------------------------------------------------------------
subroutine initialise(                                                         &
    ! Inputs
    nlevels, theta_start, q_start,                                             &
    ! Outputs
    final_level,                                                               &
    theta_adjust, theta_noq_adjust,                                            &
    theta_adj_rlx, theta_noq_adj_rlx, theta_final,                             &
    q_adjust, q_adj_rlx, q_conserve, q_dummy, q_final, qcl_inc,                &
    cloud_frac, xmin, rain_sec,                                                &
    moist_has_happened, convection_flag)
!------------------------------------------------------------------------------
! Initialise working arrays by setting them equal to input arrays.
!
! Further, the convbase and cloudbase flags are set to -1, as no
! convection has happened yet.
!
! Finally, rainsec is set to 0.
!------------------------------------------------------------------------------
implicit none

! Inputs
integer, intent(in) :: nlevels ! number of full model levels

real(kind=real_umphys), intent(in) :: theta_start(1:nlevels)
                                           ! initial pot. temp. at full levels.
real(kind=real_umphys), intent(in) :: q_start(1:nlevels)
                                       ! initial specific humidity

! Outputs
integer, intent(out) :: final_level

real(kind=real_umphys), intent(out) :: theta_adjust     (1:nlevels)
real(kind=real_umphys), intent(out) :: theta_noq_adjust (1:nlevels)
real(kind=real_umphys), intent(out) :: theta_adj_rlx    (1:nlevels)
real(kind=real_umphys), intent(out) :: theta_noq_adj_rlx(1:nlevels)
real(kind=real_umphys), intent(out) :: theta_final      (1:nlevels)
real(kind=real_umphys), intent(out) :: q_adjust         (1:nlevels)
real(kind=real_umphys), intent(out) :: q_adj_rlx        (1:nlevels)
real(kind=real_umphys), intent(out) :: q_conserve       (1:nlevels)
real(kind=real_umphys), intent(out) :: q_dummy          (1:nlevels)
real(kind=real_umphys), intent(out) :: q_final          (1:nlevels)
real(kind=real_umphys), intent(out) :: qcl_inc          (1:nlevels)
real(kind=real_umphys), intent(out) :: cloud_frac       (1:nlevels)
real(kind=real_umphys), intent(out) :: xmin             (1:nlevels)
real(kind=real_umphys), intent(out) :: rain_sec

logical, intent(out) :: convection_flag
logical, intent(out) :: moist_has_happened

theta_adjust = theta_start
theta_noq_adjust = theta_start
theta_adj_rlx = theta_start
theta_noq_adj_rlx = theta_start
theta_final = theta_start

q_adjust = q_start
q_adj_rlx = q_start
q_conserve = q_start
q_final = q_start
q_dummy = q_start

qcl_inc = 0.0
cloud_frac = 0.0
xmin = 0.0
rain_sec = 0.0

final_level = 1
moist_has_happened = .false.
convection_flag = .false.

end subroutine initialise


subroutine diagnose(                                                           &
     ! Inputs
     nlevels, current_level, theta_start, q_start, p_full,                     &
     ! In/out
     convection_flag,                                                          &
     ! Outputs
     convbase, cloudbase, final_level, moist_trigger, rh_use)
!------------------------------------------------------------------------------
! Diagnose on which levels dry or moist convection occurs.
!
! This subroutine decides whether it is possible for convection to
! originate on the model level under consideration, and if so, whether
! at any point above this level moist convection is possible (i.e
! whether a parcel rising from this level would saturate).
!
! The criterion for whether mixing is possible is the dry instability:
! where d_theta / d_z < 0.
!
! Further the requirement for a rising parcel to saturate is
! q(convbase) > rh_crit * q_sat(level where it saturates).
!------------------------------------------------------------------------------

use qsat_mod, only: qsat_wat, qsat_wat_mix

implicit none

! Inputs
integer, intent(in) :: nlevels
integer, intent(in) :: current_level

real(kind=real_umphys), intent(in) :: theta_start(1:nlevels)
real(kind=real_umphys), intent(in) :: q_start    (1:nlevels)
real(kind=real_umphys), intent(in) :: p_full     (1:nlevels)

! Input/Outputs
logical, intent(in out) :: convection_flag

! Outputs
integer, intent(out) :: convbase
integer, intent(out) :: cloudbase
integer, intent(out) :: final_level

real(kind=real_umphys), intent(out) :: rh_use

logical, intent(out) :: moist_trigger

! Local variables
integer :: j_lev
real(kind=real_umphys) :: temp
real(kind=real_umphys) :: q_sat
real(kind=real_umphys) :: q_crit
!real :: rh

! initialise cloudbase
cloudbase = -1
convbase = -1
final_level = -1
rh_use = 1.0
moist_trigger = .false.

temp = real_temp(theta_start(current_level), p_full(current_level))

if (l_mr_physics) then
  call qsat_wat_mix(q_sat, temp, p_full(current_level))
else
  call qsat_wat(q_sat, temp, p_full(current_level))
end if

q_crit = llcs_rhcrit * q_sat
moist_trigger = (q_start(current_level) >= q_crit)
if ( (theta_start(current_level) >= theta_start(current_level+1)) .or.         &
     moist_trigger ) then
  ! Thermal instability or moist instability

  if (q_start(current_level) >= 0.0) then

    convbase = current_level

    ! Convection happens, change flag
    convection_flag = .true. ! TODO: shouldn't it be outside this if-statement?

    do j_lev = current_level, nlevels - 1
      temp = real_temp(theta_start(j_lev), p_full(j_lev))

      if ( l_mr_physics ) then
        call qsat_wat_mix(q_sat, temp, p_full(j_lev))
      else
        call qsat_wat(q_sat, temp, p_full(j_lev))
      end if

      q_crit = llcs_rhcrit * q_sat
      if (q_start(j_lev) >= q_crit) then
        cloudbase = max(j_lev, convbase + 1) ! convbase cannot be cloudbase
        ! if (cloudbase == convbase) then
        !   cloudbase = convbase + 1 ! convbase cannot be cloudbase
        ! end if
        exit
      end if
    end do

    ! Allows reference profile RH to adopt RH of layer where moist
    ! convection begins.
    ! This is commented out for now pending further testing, but the code
    ! is left here for ease of future development
    ! if (cloudbase /= -1) then
    !   temp = real_temp(theta_start(cloudbase), p_full(cloudbase))
    !   call qsat_wat_mix(q_sat, temp, p_full(cloudbase), 1, .false.)
    !   rh = q_start(cloudbase) / q_sat
    !   if ((rh > rh_use) .and. (rh <= 1.0)) then
    !     rh_use = rh
    !   else if (rh > 1.0) then
    !     rh_use = 1.0
    !   end if
    ! end if

  else
    final_level = current_level + 1
  end if
else
  final_level = current_level + 1
end if

end subroutine diagnose


subroutine dry_convection(                                                     &
    ! Inputs
    nlevels, convbase,                                                         &
    ! Inputs/outputs
    theta_adjust, theta_noq_adjust, q_adjust,                                  &
    ! Outputs
    final_level)
!------------------------------------------------------------------------------
! Perform dry convection.
!
! This routine increases the potential temperature on model levels to
! those given by the dry adiabat, where there is dry instability.
!
! We also perform the same adjustment for theta_noq_adjust, a dummy 'dry
! convection only' profile that we use for enthalpy conservation.
!
! Specific humidity is also mixed in a similar way, by setting
! q_adjust(level_under_consideration) to q_start(convbase).
!
! Both of these adjustments will then be reduced to account for
! timestep / timescale.
!
! Returns
! final_level: the level where the convective event ends.
!   It is from this level that the routine 'diagnose' will be re-called, in
!   an effort to find additional thermal instabilities that would permit
!   convection.
!------------------------------------------------------------------------------

implicit none

! Inputs
integer, intent(in) :: nlevels
integer, intent(in) :: convbase

! Inputs/outputs
real(kind=real_umphys), intent(in out) :: theta_adjust    (1:nlevels)
real(kind=real_umphys), intent(in out) :: theta_noq_adjust(1:nlevels)
real(kind=real_umphys), intent(in out) :: q_adjust        (1:nlevels)

! Outputs
integer, intent(out) :: final_level

! Local variables
integer :: j_lev

final_level = -1 ! if this remains unchanged then, in circumstances
                 ! where dry_convection has been called by moist_
                 ! convection, it signals to moist_convection that
                 ! the event is able to proceed. I.e the dry mixing
                 ! reached the level cloudbase.

do j_lev = convbase + 1, nlevels
  if ( (theta_adjust(j_lev-1) >= theta_adjust(j_lev)) .and.                    &
       (q_adjust(j_lev-1) >=  0.0) ) then
    theta_adjust(j_lev) = theta_adjust(j_lev-1)
    theta_noq_adjust(j_lev) = theta_noq_adjust(j_lev-1)
    if ( (q_adjust(j_lev-1) >= 0.0) .and.                                      &
         (q_adjust(j_lev-1) > q_adjust(j_lev)) ) then
      q_adjust(j_lev) = q_adjust(j_lev-1)
    end if
  else
    final_level = j_lev
    exit
  end if
end do

end subroutine dry_convection

! ----------------------------------------------------------------------

subroutine moist_convection(                                                   &
    ! Inputs
    nlevels, convbase, cloudbase, rh_use, p_full,                              &
    ! Inputs/outputs
    theta_adjust, theta_noq_adjust, q_adjust,                                  &
    moist_has_happened, moist_trigger,                                         &
    ! Outputs
    final_level)
!------------------------------------------------------------------------------
! Perform moist convection.
!
! This routine is called if we know that a level cloudbase exists
! such that a parcel of air rising from convbase would saturate.
!
! In this situation the theta profile is adjustted beyond that of a
! dry reference profile to the reference profile for the moist
! pseudo-adiabat. This is as a result of the latent heating provided
! by saturation.
!
! The specific humidity within the columns is mixed as before, unless
! this would lead to a level having specific humidity greater than
! saturation specific humidity, in which case we set specific humidity
! equal to saturation specific humidity.
!
! If cloudbase-convbase >= 2 then we perform dry convection until
! we reach the level cloudbase.
!
! If dry convection ceases to be possible before we reach the level
! cloudbase then we claim that a parcel rising from convbase does not
! reach a level at which it saturates, and so we end the convective event
! accordingly without actually performing moist convection. Moist
! covection may still be possible in such a scenario if, when we re-call
! diagnose from 'final_level', it finds another level convbase with
! accompanying level cloudbase.
!
! Finally it returns an integer, moist_has_happened, which signals whether
! moist convection has infact taken place. This is initialised in the
! routine 'initialise' to be moist_has_happened = -1. If at any point
! moist convection takes place this is changed to 1.
!------------------------------------------------------------------------------

use qsat_mod, only: qsat_wat, qsat_wat_mix

implicit none

! Inputs
integer, intent(in) :: nlevels
integer, intent(in) :: convbase
integer, intent(in) :: cloudbase

real(kind=real_umphys), intent(in) :: rh_use
real(kind=real_umphys), intent(in) :: p_full(1:nlevels)

! Input/outputs
real(kind=real_umphys), intent(in out) :: theta_adjust    (1:nlevels)
real(kind=real_umphys), intent(in out) :: theta_noq_adjust(1:nlevels)
real(kind=real_umphys), intent(in out) :: q_adjust        (1:nlevels)

logical, intent(in out) :: moist_has_happened
logical, intent(in out) :: moist_trigger

! Outputs
integer, intent(out) :: final_level

! Locals
integer :: j_lev
real(kind=real_umphys) :: temp(1:nlevels)
real(kind=real_umphys) :: q_sat
real(kind=real_umphys) :: r_sat
real(kind=real_umphys) :: temp_v
real(kind=real_umphys) :: gamma_s
real(kind=real_umphys) :: gamma_p
real(kind=real_umphys) :: rho
real(kind=real_umphys) :: q_adjust_save

final_level = -1 ! final_level = -1 allows moist convection to happen.

temp = real_temp(theta_adjust, p_full)

if ( (cloudbase-convbase) >= 2 ) then
  ! Note that in this situation we only want to perform dry convection
  ! up until the level below cloudbase, so we use cloudbase-1 as the
  ! 'nlevels' input.
  call dry_convection(                                                         &
      ! Inputs
      cloudbase-1, convbase,                                                   &
      ! Inputs/outputs
      theta_adjust, theta_noq_adjust, q_adjust,                                &
      ! Outputs
      final_level)
  temp = real_temp(theta_adjust, p_full)
end if

if (final_level == -1) then
  ! Dry convection reached level below cloudbase, so we do moist convection.
  do j_lev = cloudbase, nlevels
    if ( (theta_adjust(j_lev-1) >= theta_adjust(j_lev))                        &
        .or. moist_trigger ) then

      moist_trigger = .false. ! reset moist trigger

      if (q_adjust(j_lev-1) >= 0.0) then
        ! Adjust specific humidity    CHANGES HERE REMOVE RH_CRITS
        if ( l_mr_physics ) then
          call qsat_wat_mix(q_sat, temp(j_lev), p_full(j_lev))
        else
          call qsat_wat(q_sat, temp(j_lev), p_full(j_lev))
        end if
        q_sat = rh_use * q_sat

        ! Save pre-adjusted spec. hum. as local scalar
        q_adjust_save = q_adjust(j_lev)

        ! Adjust spec. hum. at the current level
        q_adjust(j_lev) = max(q_adjust(j_lev), q_sat)

        ! Next, adjust potential temperature

        call qsat_wat_mix(r_sat, temp(j_lev), p_full(j_lev))
        r_sat = rh_use * r_sat

        ! Saturation dT / dz
        gamma_s = -g *                                                         &
            ( (1+r_sat) * (1+(lc*r_sat) / (r * temp(j_lev-1))) )               &
            / (cp + r_sat*hcapv +                                              &
            (repsilon + r_sat) * ((r_sat*lc**2)/(r*temp(j_lev-1)**2)))

        call qsat_wat_mix(r_sat, temp(j_lev-1), p_full(j_lev-1))
        r_sat = rh_use * r_sat

        temp_v = temp(j_lev-1) * (1 + r_sat / repsilon) / (1 + r_sat)
        ! Virtual temperature temp_v used for density of moist air
        rho = p_full(j_lev-1) / (r * temp_v) ! moist air density
        gamma_p = (1.0 / (rho * g)) * gamma_s ! saturated dT/dp
        if ( (temp(j_lev-1) +                                                  &
             (gamma_p * (p_full(j_lev-1)-p_full(j_lev)))) <= temp(j_lev) ) then
          final_level = j_lev
          q_adjust(j_lev) = q_adjust_save
          exit
        end if

        ! Reference profile temperature
        temp(j_lev) = temp(j_lev-1) +                                          &
            gamma_p * (p_full(j_lev-1) - p_full(j_lev))
        ! adjust potential temperature accordingly
        theta_adjust(j_lev) = potential_temp(temp(j_lev), p_full(j_lev))

        moist_has_happened = .true.

        ! Adjust dummy dry potential temperature profile
        theta_noq_adjust(j_lev) = max(theta_noq_adjust(j_lev-1),               &
                                      theta_noq_adjust(j_lev))
      else
        final_level = j_lev
        exit
      end if

    else ! record final level
      final_level = j_lev
      exit
    end if
  end do
end if

end subroutine moist_convection


subroutine time_relax(                                                         &
    ! Inputs
    theta_start, theta_adjust, q_start, q_adjust,                              &
    ! Inputs/outputs
    theta_adj_rlx, q_adj_rlx)
!------------------------------------------------------------------------------
! Relax temperature and humidity increments according to the conv. time scale.
!
! This routine modifies the amount by which potential temperature and
! specific humidity are increased to account for convective mixing
! timescale and model timestep.
!
! In cases where moist convection has happened, this routine is
! called twice; once for theta_adjust, and once for theta_noq_adjust.
!------------------------------------------------------------------------------

implicit none

! Inputs
real(kind=real_umphys), intent(in) :: theta_start(:)
real(kind=real_umphys), intent(in) :: theta_adjust(:)
real(kind=real_umphys), intent(in) :: q_start(:)
real(kind=real_umphys), intent(in) :: q_adjust(:)

! Outputs
real(kind=real_umphys), intent(in out) :: theta_adj_rlx(:)
real(kind=real_umphys), intent(in out) :: q_adj_rlx(:)

! Locals
real(kind=real_umphys) :: exp_fac

! Can be simplified to
! A_rlx = A_adj - exp() * (A_adj - A_start)

exp_fac = 1.0 - exp(-timestep / llcs_timescale)

theta_adj_rlx = theta_start + (theta_adjust - theta_start) * exp_fac
q_adj_rlx = q_start + (q_adjust - q_start) * exp_fac

end subroutine time_relax


subroutine conserve_moisture(                                                  &
    ! Inputs
    nlevels, delta_p, q_start, q_adj_rlx,                                      &
    ! Outputs
    q_conserve)

!------------------------------------------------------------------------------
! Conserve moisture after convection adjustment.
!
! This routine ensures that at this point in the scheme there is no net
! change in moisture within the column (this will change later, as some
! is lost through rainfall).
!
! To do this, the routine find_bounds is used to find the top and
! bottom of convective events that have occured within the column.
!
! For each event, we find the total mass of water contained on the levels
! involved in this event both before and after convection. We then take
! the ratio of water : before/after and multiply the specfic humidity on
! each level involved in the event by this ratio.
!------------------------------------------------------------------------------

implicit none

! Inputs
integer, intent(in) :: nlevels

real(kind=real_umphys), intent(in) :: delta_p  (1:nlevels)
real(kind=real_umphys), intent(in) :: q_start  (1:nlevels)
real(kind=real_umphys), intent(in) :: q_adj_rlx(1:nlevels)

! Outputs
real(kind=real_umphys), intent(out) :: q_conserve(1:nlevels)

! Locals
integer :: bottom ! used by find_bounds()
integer :: top ! used by find_bounds()
integer :: j_lev
integer :: k
integer :: n

real(kind=real_umphys) :: q_start_total_mass
real(kind=real_umphys) :: q_adj_rlx_total_mass
real(kind=real_umphys) :: q_inc_ratio
real(kind=real_umphys) :: q_start_mass  (1:nlevels)
real(kind=real_umphys) :: q_adj_rlx_mass(1:nlevels)
real(kind=real_umphys) :: q_inc         (1:nlevels)

logical :: stopnow ! used by find_bounds()

stopnow = .false.
q_conserve = q_adj_rlx

! Find increase in specific humidity after convection
! find_bounds() uses the difference between the start and adjusted profile
! to determine where a convective event begins and ends. For logistical
! reasons it requires a positive difference. For more information, see
! description of find_bounds().
q_inc = abs(q_adj_rlx - q_start)

n = 2 ! part of find_bounds set up.

do j_lev = 1, nlevels
  q_start_mass(:) = 0.0
  q_adj_rlx_mass(:) = 0.0

  if (j_lev == n) then
    call find_bounds(                                                          &
        ! Inputs
        nlevels, q_inc,                                                        &
        ! Inputs/outputs
        n, stopnow,                                                            &
        ! Outputs
        bottom, top)

    if (.not. stopnow) then
      do k = bottom, top ! find masses on each level within event
        q_start_mass(k) = abs(q_start(k)) * delta_p(k) / g
        q_adj_rlx_mass(k) = abs(q_adj_rlx(k)) * delta_p(k) / g
      end do

      ! find total masses for each event
      q_start_total_mass = sum(q_start_mass(bottom:top))
      q_adj_rlx_total_mass = sum(q_adj_rlx_mass(bottom:top))

      ! find ratio before / after
      q_inc_ratio = q_start_total_mass / q_adj_rlx_total_mass

      ! set q_conserve accordingly
      do k = bottom, top
        q_conserve(k) = abs(q_adj_rlx(k)) * q_inc_ratio
      end do

    else
      exit
    end if
  end if

end do

end subroutine conserve_moisture


subroutine rain(                                                               &
    ! Inputs
    nlevels, theta_start, theta_noq_adj_rlx, q_conserve, qcl_start, xmin,      &
    p_full, delta_p,                                                           &
    ! Inputs/outputs
    theta_adj_rlx,                                                             &
    ! Outputs
    rain_sec, q_final, qcl_inc, cloud_frac)
!------------------------------------------------------------------------------
! Calculate latent heating and rain out the equivalent amount of moisture.
!
! This routine calculates the total amount of latent heating within a
! convective event that was required for the moist convection to take
! place.
!
! It then removes the mass of water necessary to do this in the form of
! instantaneous precipitation and adjusts the specific humidity of the
! levels under consideration accordingly.
!
! Total rainfall in the column is found via Q = mass * L where Q is
! latent heating.
!------------------------------------------------------------------------------

implicit none

! Inputs
integer, intent(in) :: nlevels

real(kind=real_umphys), intent(in) :: theta_start      (1:nlevels)
real(kind=real_umphys), intent(in) :: theta_noq_adj_rlx(1:nlevels)
real(kind=real_umphys), intent(in) :: q_conserve       (1:nlevels)
real(kind=real_umphys), intent(in) :: qcl_start        (1:nlevels)
real(kind=real_umphys), intent(in) :: xmin             (1:nlevels)
real(kind=real_umphys), intent(in) :: p_full           (1:nlevels)
real(kind=real_umphys), intent(in) :: delta_p          (1:nlevels)

! Inputs/Outputs
real(kind=real_umphys), intent(in out) :: theta_adj_rlx (1:nlevels)

! Outputs
real(kind=real_umphys), intent(out) :: rain_sec
real(kind=real_umphys), intent(out) :: q_final         (1:nlevels)
real(kind=real_umphys), intent(out) :: qcl_inc         (1:nlevels)
real(kind=real_umphys), intent(out) :: cloud_frac      (1:nlevels)

! Locals
integer :: j_lev
integer :: k
integer :: n
integer :: bottom ! used by find_bounds()
integer :: top ! used by find_bounds()

real(kind=real_umphys) :: event_q_total
real(kind=real_umphys) :: total_latent
real(kind=real_umphys) :: total_latent_new
real(kind=real_umphys) :: rain_out
real(kind=real_umphys) :: rain_out_new
real(kind=real_umphys) :: latent_ratio
real(kind=real_umphys) :: temp_new
real(kind=real_umphys) :: q_level(1:nlevels)
real(kind=real_umphys) :: latent(1:nlevels)
real(kind=real_umphys) :: temp(1:nlevels)
real(kind=real_umphys) :: temp_noq(1:nlevels)
real(kind=real_umphys) :: theta_inc(1:nlevels)
real(kind=real_umphys) :: qcl_excess(1:nlevels)
real(kind=real_umphys) :: qcl_tot

logical :: stopnow ! used by find_bounds()
logical :: is_conv(1:nlevels)

stopnow = .false. ! initialise stop now
is_conv(:) = .false.
rain_sec = 0.0

q_final = q_conserve

! Find increase in potential temp after convection
! find_bounds uses the difference between the start and adjusted profile
! (as before, see description in moisture_conserve)
theta_inc = abs(theta_adj_rlx - theta_start)
qcl_inc(:) = 0.0
qcl_excess(:) = 0.0
cloud_frac(:) = 0.0
qcl_tot = 0.0

n = 2 ! part of find_bounds set up.

! calculate temperature for dummy 'only dry convection' theta profile
temp_noq = real_temp(theta_noq_adj_rlx, p_full)

do j_lev = 1, nlevels
  if (j_lev == n) then
    call find_bounds(                                                          &
        ! Inputs
        nlevels, theta_inc,                                                    &
        ! Inputs/outputs
        n, stopnow,                                                            &
        ! Outputs
        bottom, top)

    if (stopnow) then
      exit
    else
      ! calculate temperature array for actual theta profile
      temp = real_temp(theta_adj_rlx, p_full)

      ! reset all arrays and totals to 0 for new event
      latent(:) = 0.0
      q_level(:) = 0.0
      event_q_total = 0.0
      total_latent = 0.0
      total_latent_new = 0.0
      latent_ratio = 0.0
      rain_out = 0.0
      rain_out_new = 0.0

      ! Rainfall / latent heating adjustment
      do k = bottom, top
        latent(k) = delta_p(k) / g * cp * (temp(k) - temp_noq(k))
        is_conv(k) = .true.
      end do

      ! Calculate rainfall within this column from the amount
      ! of latent heating that has taken place.
      total_latent = sum(latent(bottom:top))
      if (total_latent > 0.0) then
        do k = bottom, top
          qcl_inc(k) = latent(k) * (g / delta_p(k)) / lc
          ! CHANGE HERE WAS LC*TIMESTEP
          if (qcl_inc(k) > 0.0) then
            ! Cloud fraction representing detrainment
            ! Needs to be high to emulate the effect of anvils. At the same
            ! time, because the current scheme implies that rain falls "below"
            ! cloud fraction, large values of cloud_frac lead to
            ! precipitation cores being spread out across the area of anvils,
            ! thus making precipitation too slow and evaporating too much.
            cloud_frac(k) = llcs_detrain_coef * timestep / llcs_timescale
          end if
        end do

        rain_out = total_latent / lc

        ! find total water on levels involved in convective event.
        do k = bottom, top
          q_level(k) = (delta_p(k) / g) * q_conserve(k)
        end do
        event_q_total = sum(q_level(bottom:top))

        ! The following code handles situations where the amount of rain
        ! demanded by a convective event is greater than the amount of
        ! rain stored in the levels participating in the convective event.
        ! In this scenario we limit rain_out to the amount of moisture
        ! present in the event and revise our convective potential
        ! temperature adjustment accordingly.
        if (rain_out > event_q_total) then
          rain_out_new = event_q_total
          total_latent_new = rain_out_new * lc
          latent_ratio  = total_latent_new / total_latent
          latent = latent * latent_ratio
          do k = bottom, top
            temp_new = temp_noq(k) + latent(k) / ((delta_p(k) / g) * cp)
            theta_adj_rlx(k) = potential_temp(temp_new, p_full(k))
          end do
          qcl_inc = qcl_inc * event_q_total / rain_out
          event_q_total = rain_out_new
          rain_out = rain_out_new
        end if

        ! Rescale specific humidity on levels within this event to
        ! account for rainfall
        do k = bottom, top
          q_level(k) = q_level(k) * (event_q_total - rain_out) / event_q_total
          q_final(k) = (g/delta_p(k)) * q_level(k)
        end do
      end if ! total_latent > 0
    end if ! stopnow
  end if ! j_lev == n
end do ! j_lev-loop

select case (llcs_cloud_precip)
case (llcs_opt_all_rain)
  ! produce rainfall
  rain_sec = sum(qcl_inc * delta_p / g) / timestep
  ! zero the cloud increments
  qcl_inc(:) = 0.0
  cloud_frac(:) = 0.0

case (llcs_opt_all_cloud)
  ! zero the rain-rate and hand all condensate to LS cloud
  rain_sec = 0.0

case (llcs_opt_crit_condens)
  ! Rain out excess of moisture and pass the rest of condensate to LS cloud
  do k = 1, nlevels
    qcl_tot = qcl_start(k) + qcl_inc(k)
    if ( (qcl_tot > xmin(k)) .and. is_conv(k) ) then
      ! rain_out = (qcl_tot - xmin(k)) * delta_p(k) / g
      ! Excess of qcl to be rained out
      qcl_excess(k) = qcl_tot - xmin(k)
      ! The next line is equivalent to setting qcl_tot = xmin
      qcl_inc(k) = xmin(k) - qcl_start(k)
    end if
  end do
  rain_sec = sum(qcl_excess * delta_p / g) / timestep

case (llcs_opt_const_frac)
  ! Fixed fraction of rain vs cloud
  rain_sec = sum(qcl_inc * llcs_rain_frac * delta_p / g) / timestep
  qcl_inc(:) = qcl_inc(:) * (1 - llcs_rain_frac)
  cloud_frac(:) = cloud_frac(:) * (1 - llcs_rain_frac)
end select

end subroutine rain


subroutine conserve_enthalpy(                                                  &
    ! Inputs
    nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx, p_full, delta_p,   &
    ! Outputs
    theta_final)

!------------------------------------------------------------------------------
! Conserve energy that was added to the column in the form of dry heating.
!
! Within a given convective event between the levels bottom and top
! (found by find_bounds), we find out how much energy was required to
! allow the increase in potential temperature from theta_start to
! theta_noq_adj_rlx (our dummy dry profile).
!
! In order to calculate how much energy we have added during the event
! under consideration we use dQ = cp * dT * mass.
!------------------------------------------------------------------------------

implicit none

! Inputs
integer, intent(in) :: nlevels

real(kind=real_umphys), intent(in) :: theta_start      (1:nlevels)
real(kind=real_umphys), intent(in) :: theta_adj_rlx    (1:nlevels)
real(kind=real_umphys), intent(in) :: theta_noq_adj_rlx(1:nlevels)
real(kind=real_umphys), intent(in) :: p_full           (1:nlevels)
real(kind=real_umphys), intent(in) :: delta_p          (1:nlevels)

! Inputs/Outputs
real(kind=real_umphys), intent(in out) :: theta_final  (1:nlevels)

! Locals
integer :: bottom ! used by find_bounds()
integer :: top ! used by find_bounds()
integer :: n
integer :: j_lev
integer :: k

real(kind=real_umphys) :: total_heating
real(kind=real_umphys) :: total_mass
real(kind=real_umphys) :: temp_correction
real(kind=real_umphys) :: temp_start      (1:nlevels)
real(kind=real_umphys) :: temp_adj_rlx    (1:nlevels)
real(kind=real_umphys) :: temp_noq_adj_rlx(1:nlevels)
real(kind=real_umphys) :: temp_final      (1:nlevels)
real(kind=real_umphys) :: dry_heating     (1:nlevels)
real(kind=real_umphys) :: mass            (1:nlevels)
real(kind=real_umphys) :: theta_dry_inc   (1:nlevels)

logical :: stopnow ! used by find_bounds()

! Find temperatures
temp_start = real_temp(theta_start, p_full)
temp_adj_rlx = real_temp(theta_adj_rlx, p_full)
temp_noq_adj_rlx = real_temp(theta_noq_adj_rlx, p_full)
temp_final = real_temp(theta_final, p_full)

stopnow = .false.

theta_dry_inc = theta_noq_adj_rlx - theta_start
theta_final = theta_adj_rlx

n = 2 ! required for find_bounds

do j_lev = 1, nlevels
  if (j_lev == n) then
    call find_bounds(                                                          &
        ! Inputs
        nlevels, theta_dry_inc,                                                &
        ! Inputs/outputs
        n, stopnow,                                                            &
        ! Outputs
        bottom, top)
    if (.not. stopnow) then
      do k = bottom, top
        mass(k) = 0.0
        dry_heating(k) = 0.0
        mass(k) = delta_p(k) / g
        dry_heating(k) = cp * (temp_noq_adj_rlx(k) - temp_start(k)) * mass(k)
      end do
      total_mass = sum(mass(bottom:top))
      total_heating = sum(dry_heating(bottom:top))
      temp_correction = total_heating / (cp * total_mass)
      do k = bottom, top
        temp_final(k) = temp_adj_rlx(k) - temp_correction
        theta_final(k) = potential_temp(temp_final(k), p_full(k))
      end do
    else
      exit
    end if
  end if
end do

end subroutine conserve_enthalpy


subroutine critical_condensate(                                                &
      ! Inputs
      nlevels, theta_start, p_full,                                            &
      ! Inputs/outputs
      xmin)
!------------------------------------------------------------------------------
! Calculate critical cloud condensate from initial temperature and humidity,
! given the free parameters: minimum, maximum and qsat scaling.
!
! The result is used as a threshold to separate convective rainfall and cloud
! condensed water (the latter is then passed to the cloud scheme)
! The subroutine is analogous to cloud_w_6a() and is called if
! llcs_cloud_precip is equal to llcs_opt_crit_condens.
!------------------------------------------------------------------------------
use qsat_mod, only: qsat_wat, qsat_wat_mix

implicit none

! Inputs
integer, intent(in) :: nlevels

real(kind=real_umphys), intent(in) :: theta_start(1:nlevels)
real(kind=real_umphys), intent(in) :: p_full     (1:nlevels)

! Inputs/outputs
real(kind=real_umphys), intent(in out) :: xmin    (1:nlevels)

! Locals
integer :: j_lev

real(kind=real_umphys) :: temp
real(kind=real_umphys) :: q_sat

do j_lev = 1, nlevels
  temp = real_temp(theta_start(j_lev), p_full(j_lev))
  if (l_mr_physics) then
    call qsat_wat_mix(q_sat, temp, p_full(j_lev))
  else
    call qsat_wat(q_sat, temp, p_full(j_lev))
  end if

  ! Equivalent to opt. 4 in cloud_w_6a()
  xmin(j_lev) = max(min(mparwtr, fac_qsat * q_sat), qlmin)
end do

end subroutine critical_condensate

! ----------------------------------------------------------------------
! -- end OF ROUTINES CALLED BY CONVECTION SCHEME -- !
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! -- EXTRA ROUTINES -- !
! ----------------------------------------------------------------------

subroutine find_bounds(                                                        &
    ! Inputs
    nlevels, inc,                                                              &
    ! Inputs/outputs
    n, stopnow,                                                                &
    ! Outputs
    bottom, top)
!------------------------------------------------------------------------------
! FIND BOUNDS
! This routine finds the bounds for a convective event from either the
! potential temperature or specific humidity increment for that event.
!
! Note that this routine also outputs the integer n which is the level
! above the top level for the convective event under consideration. In
! other words, this is the first level above a convective event where a
! different convective event could be found.
!------------------------------------------------------------------------------

implicit none

! Inputs
integer, intent(in) :: nlevels

real(kind=real_umphys), intent(in) :: inc(1:nlevels)

! Inputs/Outputs
integer, intent(in out) :: n

logical, intent(in out) :: stopnow

! Outputs
integer, intent(out) :: bottom, top

! Locals
integer :: j_lev

bottom = -1

do j_lev = n, nlevels
  ! find bottom
  if (bottom == -1) then
    if (inc(j_lev) > 0) then
      bottom = j_lev - 1
    else if (j_lev == nlevels) then
      stopnow = .true.
      exit
    end if
  end if

  ! find top
  if (j_lev == nlevels) then
    top = j_lev - 1
    exit
  else if ( (j_lev - 1 /= bottom) .and. (inc(j_lev) == 0) ) then
    if (bottom /= -1) then
      top = j_lev - 1
      exit
    end if
  end if

end do

n = top + 1

end subroutine find_bounds


pure elemental function real_temp(theta, p)
! Calculate real temperature.
! This routine calculates temperature from scalar values or an array of
! potential temperature and air pressure

implicit none

real(kind=real_umphys), intent(in) :: theta
real(kind=real_umphys), intent(in) :: p

real(kind=real_umphys) :: real_temp

real_temp = theta * (p / p_zero) ** kappa

end function real_temp


pure elemental function potential_temp(t, p)
! Calculate potential temperature.
! This routine calculates theta from scalar values or an array of
! temperature and air pressure

implicit none

real(kind=real_umphys), intent(in) :: t
real(kind=real_umphys), intent(in) :: p

real(kind=real_umphys) :: potential_temp

potential_temp = t * (p_zero / p) ** kappa

end function potential_temp

end module llcs
