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
! This module uses the old MURK power law scaling to compute the aerosol
! dry radius, Rd, and number concentration, Nc, from the aerosol mass mixing
! ratio am.
!

module vera_murk_scaling_mod

use vera_kind_mod,   only: wp => vera_real

use vera_global_mod, only: vera_scaling   ,                                    &
                           vera_constants ,                                    &
                           vera_scaled

implicit none

! Description:
!   This module computes the aerosol particle dry radius, Rd, and number
!   concentration, Nc, from the MURK field estimate of the aerosol mass
!   mixing ratio.
!
! Method:
!   The power law scaling scheme used in the visibility diagnostic VISBTY
!   is used to compute Rd and Nc. This module comprises just the single
!   subroutine, vera_murk_scaling.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

! name of this module
character (len=*), parameter, private :: ModuleName='VERA_MURK_SCALING_MOD'

private

public :: vera_murk_scaling  ! subroutine that produces the monodisperse
                             ! aerosol scaled off the MURK field

contains

  !=============================================================================
  !
  ! vera_murk_scaling -Power law scaling from the MURK aerosol field to
  ! estimate the aerosol particle dry radius and the aerosol number
  ! concentration, ie (Nc, Rd).
  !
  !=============================================================================

subroutine vera_murk_scaling( am )

  ! use the DrHook stuff
use parkind1,   only: jpim, jprb
use yomhook,    only: lhook, dr_hook

implicit none

!===========================================================================
! arguments for vera_murk_scaling
!===========================================================================

! aerosol mass mixing ratio, [micrograms / kilograms]
real (wp), intent(in)  :: am

!===========================================================================
! local variables for vera_murk_scaling
!===========================================================================

! computed aerosol dry radius, [m]
real (wp) :: rd

! computed aerosol number concentration, per [m^3]
real (wp) :: nc

! ratio of aerosol mass mixing ratio, am,  to the base state, am0
real (wp) :: am_over_am0

! DrHook variables
integer (kind=jpim), parameter :: zhook_in  = 0 ! DrHook tracing entry
integer (kind=jpim), parameter :: zhook_out = 1 ! DrHook tracing exit
real    (kind=jprb)            :: zhook_handle  ! DrHook tracing

! name of this routine
character (len=*)  , parameter :: RoutineName = 'VERA_MURK_SCALING'

!end of header, effect entry for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===========================================================================
! executable code for vera_murk_scaling
!===========================================================================

! compute m, the aerosol mass mixing ratio, from the base state aerosol
! mass mixing ratio m0
! NOTE: the input am is in micorgams / kg and the base state am0
! is in kg /kg so need to convert am to kg / kg using a factor
! of 1.0e-9, ie vera_constants%nano
am_over_am0 = (am / vera_scaling%am0) * vera_constants%nano

!compute the number concentration nc
nc = vera_scaling%n0                                                *          &
     ( am_over_am0**( vera_constants%one                            -          &
                      (vera_constants%three * vera_scaling%power) )   )

! compute dry radius rd
rd = vera_scaling%radius0 * ( am_over_am0**vera_scaling%power )

! update the values of the scaled rd and nc held in the vera_scaled workspace
vera_scaled%rd          = rd
vera_scaled%nc          = nc
vera_scaled%am_over_am0 = am_over_am0
vera_scaled%nc_rd_cubed = ( rd ** vera_constants%three ) * nc

!end of routine, effect exit for DrHook
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine vera_murk_scaling

end module vera_murk_scaling_mod
