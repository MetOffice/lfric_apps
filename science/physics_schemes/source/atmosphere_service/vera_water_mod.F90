! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module contains the subroutines required to compute
! the saturated specific humidity from input (temperature, pressure)
!

module vera_water_mod

use vera_kind_mod, only: wp => vera_real, um_r, um_i

implicit none

! Description:
!   This module computes the saturated specific humidity from the
!   input tempertaure and pressure.
!
! Method:
!   This module comprises four subroutines, and uses the
!   UM subroutine qsat_wat.
!
!   The four subroutines are:
!
!     vera_e_sat
!       Computes the saturated vapour pressure as a function of temperature.
!       This subroutine is only used for the standalone version of Vera. For
!       implementation in the UM, the subroutine qsat_wat is used instead.
!
!     vera_e_sat_to_q_sat
!       Computes the saturated specific humidity from the atmospheric pressure
!       and the saturated water vapour pressure.
!
!     vera_temperature_pressure_to_q_sat
!       This subroutine computes the saturated specific humidity given
!       the atmospheric pressure and temperature.
!
!     vera_perturb_q_total
!       This subroutine applies a perturbation to the total water, to mimic
!       the scheme in VISBTY - this has the effect of increasing the
!       total atmospheric water available in near saturation conditions.
!
!      This subroutine uses the VISBTY parameters PROB and RHCRIT, found in
!      the global data module as vera_visbty%prob, vera_visbty%rhcrit.
!
!      Whether or not this subroutine is used is controlled by a switch in
!      the Vera global data module, vera_visbty%switch_q_total.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_WATER_MOD'

private

! these two subroutines are called from vera_scheme_mod.F90
public :: vera_temperature_pressure_to_q_sat,                                  &
          vera_perturb_q_total

contains

  !=============================================================================
  !
  ! vera_e_sat
  !   This subroutine computes the saturated vapour pressure as a function
  !   of the temperature.
  !
  !=============================================================================

subroutine vera_e_sat( temperature, e_sat )

  !
  ! This subroutine calculates the saturated vapour pressure for
  ! a given temperature. It is assumed that if the temperature is
  ! below 0 Celsius then the saturated vapour pressure over ice is
  ! required.
  !
  ! The computation uses the parametrization of Murray, F.W. 1967
  ! "On the Computation of Saturation Vapor Pressure" J. Appl. Meteor. 6,
  ! p 203-204. This gives the saturation vapour pressure, es, as
  !
  ! es = 6.1078EXP[ a(T-273.16) / (T-b) ]
  !
  ! where T is the temperature in Kelvin and constants a and b are:
  !
  ! a = 17.2693882, b = 35.86 over water
  ! a = 21.8745584, b = 7.66  over ice
  !
  !
  ! INPUT:
  !    temperature - the temperature, in Kelvin, at which to
  !                  calculate the saturation vapour pressure.
  !
  ! OUTPUT:
  !    e_sat - the saturation vapour pressure, in Pa
  !
  ! EXAMPLE use:
  !              call vera_e_sat( -10.0+273.16, e_sat )
  !
  !              then e_sat => 259.457
  !
  ! RESTRICTIONS ON USAGE: This parametrization is designed for a temperature
  !                        range of -50 to + 50 Celsius.
  !

use vera_global_mod, only: vera_constants, vera_physics

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_e_sat
!===========================================================================

! input temperature, in Kelvin
real (wp), intent(in)  :: temperature

! output saturation vapour pressure, in Pa
real (wp), intent(out) :: e_sat

!===========================================================================
! local variables for vera_e_sat
!===========================================================================

! set the parametrization constants for water
real (wp), parameter   :: aw = 17.2693882_wp
real (wp), parameter   :: bw = 35.86_wp

! set the parametrization constants for ice
real (wp), parameter   :: ai = 21.8745584_wp
real (wp), parameter   :: bi = 7.66_wp

! set the scaling parametrization constant
real (wp), parameter   :: param_scaling = 6.1078_wp

! parametrization constants used, depending on the temperature
real (wp)              :: a
real (wp)              :: b

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_E_SAT'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_e_sat
!===========================================================================

! decide whether to use the ice, or liquid water parametrization
if ( temperature < vera_physics%water_freeze ) then

  ! use the ice values
  a = ai
  b = bi

else

  ! use the liquid water values
  a = aw
  b = bw

end if

! compute the saturation vapour pressure, in hPa
e_sat = param_scaling                                      *                   &
        exp( a * (temperature - vera_physics%water_freeze) /                   &
             ( temperature - b )                             )

! convert the saturation vapour pressure to Pa
e_sat = e_sat * vera_constants%hundred

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_e_sat

!=============================================================================
!
! vera_e_sat_to_q_sat
!   This subroutine computes the saturated specific humidity given
!   the atmospheric pressure and the saturated water vapour pressure.
!
!
!=============================================================================

subroutine vera_e_sat_to_q_sat( pressure, e_sat, q_sat )

use vera_global_mod, only: vera_physics

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_e_sat_to_q_sat
!===========================================================================

! input atmospheric pressure in Pa
real (wp), intent(in)  :: pressure

! input saturated vapour pressure in Pa
real (wp), intent(in)  :: e_sat

! output saturation specific humidity in kg / kg
real (wp), intent(out) :: q_sat

!===========================================================================
! local variables for vera_e_sat_to_q_sat
!===========================================================================

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_E_SAT_TO_Q_SAT'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_e_sat_to_q_sat
!===========================================================================

q_sat = ( ( vera_physics%mw_water / vera_physics%mw_dry_air ) *                &
           e_sat / ( pressure - e_sat)                          )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_e_sat_to_q_sat

!=============================================================================
!
! vera_temperature_pressure_to_q_sat
!   This subroutine computes the saturated specific humidity, q_sat, given
!   the atmospheric pressure and temperature.
!
!=============================================================================

subroutine vera_temperature_pressure_to_q_sat( temperature, pressure, q_sat )

use vera_global_mod, only: vera_water

! use the the new UM qsat routines
use qsat_mod,        only: qsat_wat_um => qsat_wat

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_temperature_pressure_to_q_sat
!===========================================================================

! input atmospheric temperature in Kelvin
real (wp), intent(in)  :: temperature

! input atmospheric pressure in Pa
real (wp), intent(in)  :: pressure

! output saturation specific humidity in kg / kg
real (wp), intent(out) :: q_sat

!===========================================================================
! local variables for vera_temperature_pressure_to_q_sat
!===========================================================================

! saturated vapour pressure in Pa
real    (wp)                   :: e_sat

! inputs for the UM subroutine qsat_wat_um
real    (um_r)                 :: qsat_wat_t, qsat_wat_p

! output for the UM subroutine qsat_wat_um
real    (um_r)                 :: qsat_wat_esat

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName =                                &
                                  'VERA_TEMPERATURE_PRESSURE_TO_Q_SAT'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_temperature_pressure_to_q_sat
!===========================================================================

! decide whether to use the UM routines to compute q_sat, or the equivalent
! routine in this module
if ( vera_water%switch_qsat_um == vera_water%switch_qsat_um_on ) then

  ! using the UM routines to compute q_sat, so
  ! pop the temperature and pressure into variables for input into
  ! the UM subroutine qsat_wat_um
  qsat_wat_t = temperature
  qsat_wat_p = pressure

  ! call the UM subroutine qsat_wat_um to compute the saturation
  ! specific humidity, q_sat
  call qsat_wat_um( qsat_wat_esat, qsat_wat_t, qsat_wat_p )

  ! capture the output q_sat from subroutine qsat_wat
  q_sat = qsat_wat_esat

else

  ! using the local routines to compute q_sat, so start with
  ! computing the saturation vapour pressure
  call vera_e_sat( temperature, e_sat )

  ! compute the saturation specific humidity q_sat
  call vera_e_sat_to_q_sat( pressure, e_sat, q_sat )

end if

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_temperature_pressure_to_q_sat

!=============================================================================
!
! vera_perturb_q_total
!   This subroutine applies a perturbation to the total water, to mimic
!   the scheme in VISBTY - this has the effect of increasing the
!   total atmospheric water available in near saturation conditions.
!
!   This subroutine uses the VISBTY parameters PROB and RHCRIT, found in
!   the global data module as vera_visbty%prob, vera_visbty%rhcrit.
!
!   Whether or not this subroutine is used is controlled by a switch in the
!   Vera global data module, vera_visbty%switch_q_total.
!
!=============================================================================

subroutine vera_perturb_q_total( q_total, q_sat, q_total_perturbed )

use vera_global_mod, only: vera_visbty, vera_constants

! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_perturb_q_total
!===========================================================================

! input atmospheric total water, ie q + ql
real (wp), intent(in)  :: q_total

! input atmospheric saturation specific humidity
real (wp), intent(in)  :: q_sat

! output perturbed total water in kg / kg
real (wp), intent(out) :: q_total_perturbed

!===========================================================================
! local variables for vera_perturb_q_total
!===========================================================================

! a function of the PROB global variable vera_visbty%prob
real (wp)              :: f_prob

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_PERTURB_Q_TOTAL'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_perturb_q_total
!===========================================================================

if ( vera_visbty%prob < vera_constants%half ) then

  f_prob = vera_constants%one                           -                      &
           sqrt( vera_constants%two * vera_visbty%prob )

else

  f_prob = -( vera_constants%one                             -                 &
              sqrt( vera_constants%two                       *                 &
                  ( vera_constants%one - vera_visbty%prob ) )  )

end if

q_total_perturbed = q_total +                                                  &
                    ( (vera_constants%one - vera_visbty%rhcrit) *              &
                      f_prob * q_sat                              )

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_perturb_q_total

end module vera_water_mod
